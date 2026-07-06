"""
Cluster the best kept channel from each structure based on R-value shape profile.

By default the script reads results/tunnel_selection.csv (written by run_batch.py) and
picks the first 'kept' tunnel per structure — i.e. the highest-throughput tunnel that
passed all active filters (pocket distance, curvature, length).  Use --tunnel N to
override and always cluster a fixed tunnel ID instead.

Each tunnel is resampled to N points along a normalized distance axis [0, 1],
making profiles of different lengths directly comparable. Clustering is done
with hierarchical/agglomerative clustering; a dendrogram helps choose the
number of clusters without committing upfront.
"""

import argparse
import csv
import math
import os
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import interp1d
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.spatial.distance import pdist
import matplotlib.pyplot as plt
import matplotlib.cm as cm

sys.path.insert(0, os.path.dirname(__file__))
from extract_tunnel_R import get_tunnel_axis_values
from run_batch import load_pocket_selection

ROOT = Path(__file__).parent.parent


def _shade_region(ax, regions):
    """Grey-highlight the normalized-distance window(s) used for region-restricted clustering."""
    for lo, hi in regions or []:
        ax.axvspan(lo, hi, color='gray', alpha=0.15, zorder=0)


def load_best_kept_tunnels(csv_path):
    """
    Read tunnel_selection.csv and return {stem: tunnel_id} for the first 'kept'
    tunnel per structure (already sorted by throughput descending by run_batch.py).
    """
    best = {}
    if not csv_path.exists():
        return best
    with open(csv_path, newline='') as f:
        for row in csv.DictReader(f):
            stem = row.get('stem', '').strip()
            if stem and row.get('filter_status', '').strip() == 'kept' and stem not in best:
                try:
                    best[stem] = int(row['tunnel'])
                except (KeyError, ValueError):
                    pass
    return best


def load_tunnel(results_dir, stem, tunnel_id):
    """Return (r_array, dist_array, length_angstrom) or None if not found."""
    csv_path = results_dir / stem / 'out' / 'analysis' / 'tunnel_profiles.csv'
    if not csv_path.exists():
        return None
    r    = get_tunnel_axis_values(str(csv_path), tunnel_id, 'R')
    dist = get_tunnel_axis_values(str(csv_path), tunnel_id, 'distance')
    if not r or not dist:
        return None
    n = min(len(r), len(dist))
    return np.array(r[:n]), np.array(dist[:n]), dist[n - 1]


def build_feature_vector(r, dist, n_points):
    """
    Resample R onto n_points uniform positions in normalized distance [0, 1],
    then append the exact entry, exit, and bottleneck radii so they are never
    lost to grid discretization.  All values share the same unit (Å).
    """
    norm = dist / dist[-1]
    f = interp1d(norm, r, kind='linear', bounds_error=False, fill_value=(r[0], r[-1]))
    resampled = f(np.linspace(0, 1, n_points))
    # Pinned key features: entry (near heme), exit, bottleneck (narrowest point)
    pinned = np.array([r[0], r[-1], r.min()])
    return np.concatenate([resampled, pinned])


def main():
    parser = argparse.ArgumentParser(
        description='Cluster tunnel-1 R profiles across all structures in results/.'
    )
    parser.add_argument('--results-dir', type=Path, default=ROOT / 'results',
                        help='Root results folder (default: results/)')
    parser.add_argument('--fpocket-dir', type=Path, default=None,
                        help='fpocket outputs folder — if given, only PASS structures are clustered')
    parser.add_argument('--tunnel',      type=int, default=None,
                        help='Override: always cluster this tunnel ID (default: use best kept tunnel from tunnel_selection.csv)')
    parser.add_argument('--n-points',    type=int, default=100,
                        help='Points to resample each profile to (default: 100)')
    parser.add_argument('--n-clusters',  type=int, default=None,
                        help='Number of clusters to colour (default: read from dendrogram)')
    parser.add_argument('--out',         type=Path, default=ROOT / 'results' / 'clustering',
                        help='Output directory for plots (default: results/clustering/)')
    parser.add_argument('--ylim',        type=float, default=5,
                        help='Y-axis half-range for overlay plot (default: 5)')
    parser.add_argument('--region',      type=float, nargs=2, metavar=('LO', 'HI'),
                        action='append', default=None,
                        help='Cluster using only the normalized-distance window [LO, HI] '
                             '(0 = heme, 1 = exit). Repeat the flag to combine several windows — '
                             'e.g. --region 0.08 0.24 --region 0.69 1.00 clusters on their union. '
                             'Plots still show the full profiles; the pinned entry/exit/bottleneck '
                             'features are dropped so clustering depends solely on the chosen region(s).')
    args = parser.parse_args()

    # ── fpocket allow-list ────────────────────────────────────────────────────
    pocket_selection = {}
    if args.fpocket_dir:
        pocket_selection = load_pocket_selection(args.fpocket_dir)
        print(f'fpocket PASS structures: {", ".join(pocket_selection) or "none"}')

    # ── tunnel selection: best kept per structure, or fixed override ──────────
    tunnel_csv = args.results_dir / 'tunnel_selection.csv'
    best_kept = load_best_kept_tunnels(tunnel_csv)
    if best_kept and args.tunnel is None:
        print(f'Using best kept tunnel per structure from {tunnel_csv.name}')
    elif args.tunnel is not None:
        print(f'Using fixed tunnel ID {args.tunnel} for all structures (--tunnel override)')
    else:
        print(f'tunnel_selection.csv not found — falling back to tunnel 1 for all structures')

    # ── collect profiles ──────────────────────────────────────────────────────
    stems, profiles, lengths, tunnel_ids_used = [], [], [], []
    for stem in sorted(p.name for p in args.results_dir.iterdir() if p.is_dir()
                       and p.name != 'clustering'):
        if args.fpocket_dir and stem not in pocket_selection:
            print(f'  [{stem}] skipped — not in fpocket PASS list')
            continue

        if args.tunnel is not None:
            tunnel_id = args.tunnel
        elif stem in best_kept:
            tunnel_id = best_kept[stem]
        else:
            tunnel_id = 1

        result = load_tunnel(args.results_dir, stem, tunnel_id)
        if result is None:
            print(f'  [{stem}] tunnel {tunnel_id} not found — skipping')
            continue
        r, dist, length = result
        resampled = build_feature_vector(r, dist, args.n_points)
        stems.append(stem)
        profiles.append(resampled)
        lengths.append(length)
        tunnel_ids_used.append(tunnel_id)
        print(f'  [{stem}] tunnel {tunnel_id} — length {length:.1f} Å, {len(r)} points')

    if len(stems) < 2:
        print('Need at least 2 structures with tunnel data to cluster.')
        sys.exit(1)

    X = np.vstack(profiles)          # shape (n_structures, n_points + 3 pinned)
    args.out.mkdir(parents=True, exist_ok=True)

    # ── select features for clustering (optionally normalized-distance region(s)) ─
    region = [sorted(pair) for pair in args.region] if args.region else None
    if region:
        grid = np.linspace(0, 1, args.n_points)
        mask = np.zeros(args.n_points, dtype=bool)
        for lo, hi in region:
            mask |= (grid >= lo) & (grid <= hi)
        regions_str = ', '.join(f'[{lo:.3f}, {hi:.3f}]' for lo, hi in region)
        if mask.sum() < 2:
            sys.exit(f'Region(s) {regions_str} cover fewer than 2 resampled points '
                     f'(n_points={args.n_points}); widen the window(s) or raise --n-points.')
        X_cluster = X[:, :args.n_points][:, mask]   # in-region points only, no pinned features
        print(f'Region-restricted clustering: {regions_str} '
              f'→ {int(mask.sum())} of {args.n_points} points (union)')
    else:
        X_cluster = X

    # ── hierarchical clustering ───────────────────────────────────────────────
    Z = linkage(X_cluster, method='ward', metric='euclidean')

    # ── choose cluster count (also used to colour the dendrogram) ─────────────
    n_clusters = args.n_clusters if args.n_clusters else max(2, len(stems) // 3)

    # ── dendrogram ───────────────────────────────────────────────────────────
    # Colour the tree to match the n_clusters cut: threshold just below the
    # (n_clusters-1)th largest merge so exactly n_clusters clades are coloured.
    color_thresh = Z[-(n_clusters - 1), 2] if n_clusters > 1 else 0
    dend_labels = [f'{s}  (t{t})' for s, t in zip(stems, tunnel_ids_used)]
    dend_w = max(10, len(stems) * 0.5)
    dend_h = max(8, len(stems) * 0.45)
    fig, ax = plt.subplots(figsize=(dend_w, dend_h))
    dendrogram(Z, labels=dend_labels, ax=ax, orientation='left', leaf_font_size=9,
               color_threshold=color_thresh)
    region_note = ('  [regions ' + ', '.join(f'{lo:.2f}–{hi:.2f}' for lo, hi in region) + ']'
                   ) if region else ''
    ax.set_title('Best kept tunnel per structure — hierarchical clustering (Ward linkage)'
                 + region_note)
    ax.set_xlabel('Distance')
    fig.tight_layout()
    dend_path = args.out / 'best_tunnel_dendrogram.png'
    fig.savefig(dend_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'\nDendrogram saved to {dend_path}')

    # ── assign cluster labels ─────────────────────────────────────────────────
    labels = fcluster(Z, n_clusters, criterion='maxclust')
    colors = [cm.tab10(i / 10) for i in range(n_clusters)]

    print(f'\nCluster assignments (n={n_clusters}):')
    for stem, label, length, tid in zip(stems, labels, lengths, tunnel_ids_used):
        print(f'  cluster {label}  {stem}  tunnel {tid}  (length {length:.1f} Å)')

    # ── persist cluster assignments (stem → cluster) for downstream joins ──────
    assign_path = args.out / 'cluster_assignments.csv'
    with open(assign_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['stem', 'tunnel', 'cluster', 'length'])
        for stem, label, length, tid in zip(stems, labels, lengths, tunnel_ids_used):
            writer.writerow([stem, tid, int(label), f'{length:.2f}'])
    print(f'Cluster assignments saved to {assign_path}')

    # ── overlay plot: all profiles coloured by cluster ────────────────────────
    x_norm = np.linspace(0, 1, args.n_points)

    fig, ax = plt.subplots(figsize=(10, 4))
    for i, (stem, profile, label) in enumerate(zip(stems, profiles, labels)):
        shape = profile[:args.n_points]   # strip pinned features for plotting
        c = colors[label - 1]
        ax.plot(x_norm,  shape, color=c, alpha=0.7, linewidth=1.5,
                label=f'Cluster {label}' if stem == stems[labels.tolist().index(label)] else '_')
        ax.plot(x_norm, -shape, color=c, alpha=0.7, linewidth=1.5)

    ax.axhline(0, color='black', linewidth=0.5, linestyle='--')
    ax.set_ylim(-args.ylim, args.ylim)
    ax.set_xlabel('Normalized distance along tunnel (0 = heme, 1 = exit)')
    ax.set_ylabel('R (Å)')
    ax.set_title(f'Best kept tunnel per structure — overlay by cluster (n={n_clusters})')

    handles = [plt.Line2D([0], [0], color=colors[i], linewidth=2, label=f'Cluster {i+1}')
               for i in range(n_clusters)]
    ax.legend(handles=handles, loc='upper right')
    _shade_region(ax, region)
    fig.tight_layout()
    overlay_path = args.out / 'best_tunnel_overlay.png'
    fig.savefig(overlay_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'Overlay plot saved to {overlay_path}')

    # ── overlay plot: all profiles, single colour ("before" — pre-clustering) ──
    fig, ax = plt.subplots(figsize=(4, 3))
    for profile in profiles:
        shape = profile[:args.n_points]   # strip pinned features for plotting
        ax.plot(x_norm,  shape, color='steelblue', alpha=0.6, linewidth=1)
        ax.plot(x_norm, -shape, color='steelblue', alpha=0.6, linewidth=1)
    ax.axhline(0, color='black', linewidth=0.5, linestyle='--')
    ax.set_ylim(-args.ylim, args.ylim)
    ax.set_xlabel('Normalized distance along tunnel (0 = heme, 1 = exit)')
    ax.set_ylabel('R (Å)')
    _shade_region(ax, region)
    fig.tight_layout()
    plain_path = args.out / 'best_tunnel_overlay_uncolored.png'
    fig.savefig(plain_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'Uncoloured overlay plot saved to {plain_path}')

    # ── overlay plot: coloured by cluster, no legend/labels ("after", clean) ───
    fig, ax = plt.subplots(figsize=(4, 3))
    for profile, label in zip(profiles, labels):
        shape = profile[:args.n_points]   # strip pinned features for plotting
        c = colors[label - 1]
        ax.plot(x_norm,  shape, color=c, alpha=0.7, linewidth=1)
        ax.plot(x_norm, -shape, color=c, alpha=0.7, linewidth=1)
    ax.axhline(0, color='black', linewidth=0.5, linestyle='--')
    ax.set_ylim(-args.ylim, args.ylim)
    ax.set_xlabel('Normalized distance along tunnel (0 = heme, 1 = exit)')
    ax.set_ylabel('R (Å)')
    _shade_region(ax, region)
    fig.tight_layout()
    nolabel_path = args.out / 'best_tunnel_overlay_nolabels.png'
    fig.savefig(nolabel_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'Unlabelled cluster overlay saved to {nolabel_path}')

    # ── per-cluster mean shape (grid) ────────────────────────────────────────
    ncols = min(n_clusters, max(2, math.ceil(math.sqrt(n_clusters))))
    nrows = math.ceil(n_clusters / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows), sharey=True)
    axes_flat = np.array(axes).flatten() if n_clusters > 1 else [axes]
    for ci in range(n_clusters):
        ax = axes_flat[ci]
        members = [profiles[i][:args.n_points] for i, lbl in enumerate(labels) if lbl == ci + 1]
        member_names = [stems[i] for i, lbl in enumerate(labels) if lbl == ci + 1]
        mean_r = np.mean(members, axis=0)
        std_r  = np.std(members, axis=0)
        c = colors[ci]
        for m in members:
            ax.plot(x_norm,  m, color=c, alpha=0.3, linewidth=1)
            ax.plot(x_norm, -m, color=c, alpha=0.3, linewidth=1)
        ax.plot(x_norm,  mean_r, color=c, linewidth=2)
        ax.plot(x_norm, -mean_r, color=c, linewidth=2)
        ax.fill_between(x_norm, -(mean_r + std_r), mean_r + std_r, color=c, alpha=0.15)
        ax.axhline(0, color='black', linewidth=0.5, linestyle='--')
        ax.set_ylim(-args.ylim, args.ylim)
        ax.set_title(f'Cluster {ci+1}  (n={len(members)})\n' +
                     '\n'.join(member_names), fontsize=8)
        ax.set_xlabel('Normalized distance')
        if ci % ncols == 0:
            ax.set_ylabel('R (Å)')
    # hide any unused grid cells
    for ci in range(n_clusters, len(axes_flat)):
        axes_flat[ci].set_visible(False)
    for a in axes_flat[:n_clusters]:
        _shade_region(a, region)
    fig.suptitle('Best kept tunnel per structure — per-cluster mean ± std', y=1.01)
    fig.tight_layout()
    mean_path = args.out / 'best_tunnel_cluster_means.png'
    fig.savefig(mean_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'Cluster means plot saved to {mean_path}')

    # ── per-cluster true overlay (grid, individual channels only, no mean/std) ─
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows), sharey=True)
    axes_flat = np.array(axes).flatten() if n_clusters > 1 else [axes]
    for ci in range(n_clusters):
        ax = axes_flat[ci]
        members = [profiles[i][:args.n_points] for i, lbl in enumerate(labels) if lbl == ci + 1]
        member_names = [stems[i] for i, lbl in enumerate(labels) if lbl == ci + 1]
        c = colors[ci]
        for m in members:
            ax.plot(x_norm,  m, color=c, alpha=0.7, linewidth=1.2)
            ax.plot(x_norm, -m, color=c, alpha=0.7, linewidth=1.2)
        ax.axhline(0, color='black', linewidth=0.5, linestyle='--')
        ax.set_ylim(-args.ylim, args.ylim)
        ax.set_title(f'Cluster {ci+1}  (n={len(members)})\n' +
                     '\n'.join(member_names), fontsize=8)
        ax.set_xlabel('Normalized distance')
        if ci % ncols == 0:
            ax.set_ylabel('R (Å)')
    # hide any unused grid cells
    for ci in range(n_clusters, len(axes_flat)):
        axes_flat[ci].set_visible(False)
    for a in axes_flat[:n_clusters]:
        _shade_region(a, region)
    fig.suptitle('Best kept tunnel per structure — per-cluster channel overlay (no mean/std)', y=1.01)
    fig.tight_layout()
    overlay_grid_path = args.out / 'best_tunnel_cluster_overlays.png'
    fig.savefig(overlay_grid_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'Per-cluster overlay grid saved to {overlay_grid_path}')

    # ── individual cluster figures ────────────────────────────────────────────
    for ci in range(n_clusters):
        members = [profiles[i][:args.n_points] for i, lbl in enumerate(labels) if lbl == ci + 1]
        member_names = [stems[i] for i, lbl in enumerate(labels) if lbl == ci + 1]
        mean_r = np.mean(members, axis=0)
        std_r  = np.std(members, axis=0)
        c = colors[ci]
        fig, ax = plt.subplots(figsize=(7, 4))
        for m in members:
            ax.plot(x_norm,  m, color=c, alpha=0.3, linewidth=1)
            ax.plot(x_norm, -m, color=c, alpha=0.3, linewidth=1)
        ax.plot(x_norm,  mean_r, color=c, linewidth=2)
        ax.plot(x_norm, -mean_r, color=c, linewidth=2)
        ax.fill_between(x_norm, -(mean_r + std_r), mean_r + std_r, color=c, alpha=0.15)
        ax.axhline(0, color='black', linewidth=0.5, linestyle='--')
        ax.set_ylim(-args.ylim, args.ylim)
        ax.set_xlabel('Normalized distance along tunnel (0 = heme, 1 = exit)')
        ax.set_ylabel('R (Å)')
        ax.set_title(f'Cluster {ci+1}  (n={len(members)}): ' + ', '.join(member_names),
                     fontsize=9)
        _shade_region(ax, region)
        fig.tight_layout()
        single_path = args.out / f'cluster_{ci+1}.png'
        fig.savefig(single_path, dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f'  Cluster {ci+1} figure saved to {single_path}')

        # overlay-only variant: individual channels, no mean/std
        fig, ax = plt.subplots(figsize=(7, 4))
        for m in members:
            ax.plot(x_norm,  m, color=c, alpha=0.7, linewidth=1.2)
            ax.plot(x_norm, -m, color=c, alpha=0.7, linewidth=1.2)
        ax.axhline(0, color='black', linewidth=0.5, linestyle='--')
        ax.set_ylim(-args.ylim, args.ylim)
        ax.set_xlabel('Normalized distance along tunnel (0 = heme, 1 = exit)')
        ax.set_ylabel('R (Å)')
        ax.set_title(f'Cluster {ci+1}  (n={len(members)}) — channel overlay: '
                     + ', '.join(member_names), fontsize=9)
        _shade_region(ax, region)
        fig.tight_layout()
        overlay_single_path = args.out / f'cluster_{ci+1}_overlay.png'
        fig.savefig(overlay_single_path, dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f'  Cluster {ci+1} overlay figure saved to {overlay_single_path}')


if __name__ == '__main__':
    main()
