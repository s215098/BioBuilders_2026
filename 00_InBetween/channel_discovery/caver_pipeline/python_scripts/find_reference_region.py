"""
Find the channel region best suited for region-restricted clustering, using a small
set of experimentally characterized reference structures (default: 5OXU_heme,
5FUJ_heme, 9HE9_heme).

Two views are produced:

  * SIMPLE  (reference_region_simple.png) — the region where the profiles are most
    far apart, measured by per-position spread (std or range). Easy to explain, but
    "far apart" is not the same as "informative".

  * SMOOTHNESS (reference_region.png) — the region where the lines are locally smooth.
    The idea: a region is low-information when the lines independently jitter up and down
    (no clean pattern), and high-information when each line follows a smooth local trend —
    flat OR steadily rising — regardless of how far apart they are (magnitude is ignored).
    Each line is split into a Gaussian-smoothed trend plus a residual 'wiggle'; roughness
    is the local RMS of that residual (averaged over lines) and smoothness = 1 - normalized
    roughness, in [0, 1]. A flat-but-separated region scores high; a jittery region scores
    low.

The pairwise absolute-difference panel is kept in both as a check that the chosen window
separates ALL references, not just one outlier.

Because there are only a few references (n ≈ 3), this is deliberately descriptive.
Reuses the exact loading/normalization/resampling of cluster_tunnels.py.
"""

import argparse
import csv
import os
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy.ndimage import gaussian_filter1d
import matplotlib.pyplot as plt
import matplotlib.cm as cm

sys.path.insert(0, os.path.dirname(__file__))
from cluster_tunnels import load_tunnel, build_feature_vector, load_best_kept_tunnels

ROOT = Path(__file__).parent.parent


def _norm(a):
    """Min-max normalize to [0, 1]; return all-ones if the array is ~constant."""
    a = np.asarray(a, float)
    lo, hi = np.nanmin(a), np.nanmax(a)
    if hi - lo < 1e-9:
        return np.ones_like(a)
    return (a - lo) / (hi - lo)


def local_smoothness(M, sigma):
    """
    Per-position smoothness of the *pairwise differences* between lines, in [0, 1]
    (1 = smoothest, 0 = roughest). For every pair of lines we take their signed
    difference (a - b), split it into a Gaussian-smoothed trend plus a residual 'wiggle',
    and measure the local RMS of that residual (averaged over pairs). smoothness =
    1 - normalized roughness.

    This captures "do the lines move together or separately": if the lines move together
    (correlated wiggle) their differences stay smooth → high; if they independently jitter
    up and down their differences are rough → low. A flat-but-separated region has a
    constant difference → zero roughness → high, so magnitude is ignored.
    Returns (smoothness, raw_roughness).
    """
    pairs = list(combinations(range(M.shape[0]), 2))
    diffs = np.array([M[a] - M[b] for a, b in pairs])          # (n_pairs, n)
    trend = gaussian_filter1d(diffs, sigma=sigma, axis=1)
    resid = diffs - trend
    roughness = np.sqrt(gaussian_filter1d((resid ** 2).mean(axis=0), sigma))
    return 1.0 - _norm(roughness), roughness


def peak_window(score, frac):
    """Contiguous window around the peak where score >= frac * peak."""
    pk = int(np.argmax(score))
    th = frac * score[pk]
    lo, hi = pk, pk
    while lo > 0 and score[lo - 1] >= th:
        lo -= 1
    while hi < len(score) - 1 and score[hi + 1] >= th:
        hi += 1
    return pk, lo, hi, th


def main():
    parser = argparse.ArgumentParser(
        description='Locate the best region for region-restricted clustering from references.'
    )
    parser.add_argument('--results-dir', type=Path, default=ROOT / 'results',
                        help='Root results folder (default: results/)')
    parser.add_argument('--references', default='5OXU_heme,5FUJ_heme,9HE9_heme',
                        help='Comma-separated reference stems')
    parser.add_argument('--tunnel', type=int, default=None,
                        help='Force this tunnel ID for all references (default: best kept)')
    parser.add_argument('--n-points', type=int, default=100,
                        help='Points to resample each profile to (default: 100)')
    parser.add_argument('--spread-metric', choices=['std', 'range'], default='std',
                        help='Magnitude-of-separation measure for the simple view (default: std)')
    parser.add_argument('--smooth-sigma', type=float, default=3.0,
                        help='Gaussian sigma (in points) for the trend/wiggle split used by '
                             'the smoothness score (default: 3)')
    parser.add_argument('--window-frac', type=float, default=0.5,
                        help='Window = contiguous region around the peak where the score >= this '
                             'fraction of the peak (default: 0.5)')
    parser.add_argument('--ylim', type=float, default=5,
                        help='Y-axis half-range for the overlay panel (default: 5)')
    parser.add_argument('--out', type=Path, default=ROOT / 'results' / 'clustering',
                        help='Output directory (default: results/clustering/)')
    args = parser.parse_args()

    refs = [s.strip() for s in args.references.split(',') if s.strip()]
    best_kept = load_best_kept_tunnels(args.results_dir / 'tunnel_selection.csv')

    # ── collect resampled reference profiles ──────────────────────────────────
    x_norm = np.linspace(0, 1, args.n_points)
    profiles, used_ids = {}, {}
    for stem in refs:
        tid = args.tunnel if args.tunnel is not None else best_kept.get(stem, 1)
        result = load_tunnel(args.results_dir, stem, tid)
        if result is None:
            print(f'  [{stem}] tunnel {tid} not found — skipping')
            continue
        r, dist, length = result
        profiles[stem] = build_feature_vector(r, dist, args.n_points)[:args.n_points]
        used_ids[stem] = tid
        print(f'  [{stem}] tunnel {tid} — length {length:.1f} Å')

    present = list(profiles)
    if len(present) < 2:
        sys.exit('Need at least 2 reference structures with tunnel data.')
    M = np.vstack([profiles[s] for s in present])           # (k, n_points)

    # ── two per-position scores ───────────────────────────────────────────────
    spread = M.std(axis=0) if args.spread_metric == 'std' else (M.max(axis=0) - M.min(axis=0))
    smoothness, roughness = local_smoothness(M, args.smooth_sigma)

    s_pk, s_lo, s_hi, s_th = peak_window(spread, args.window_frac)       # simple (magnitude)
    t_pk, t_lo, t_hi, t_th = peak_window(smoothness, args.window_frac)   # smoothness
    swin = (x_norm[s_lo], x_norm[s_hi])
    twin = (x_norm[t_lo], x_norm[t_hi])

    print(f'\nReferences used: {", ".join(present)} (n={len(present)})')
    print(f'Smoothness window (magnitude-free, wiggle → low): peak {x_norm[t_pk]:.3f}, '
          f'[{twin[0]:.3f}, {twin[1]:.3f}]  → --region {twin[0]:.2f} {twin[1]:.2f}')
    print(f'Simple window     (spread only):                  peak {x_norm[s_pk]:.3f}, '
          f'[{swin[0]:.3f}, {swin[1]:.3f}]  → --region {swin[0]:.2f} {swin[1]:.2f}')

    # ── pairwise |difference| curves (one-vs-rest / crossing check) ───────────
    pair_diffs = {}
    for a, b in combinations(present, 2):
        d = np.abs(profiles[a] - profiles[b])
        pair_diffs[(a, b)] = d
        print(f'  |{a} − {b}| max at {x_norm[int(np.argmax(d))]:.3f} (Δ={d.max():.3f} Å)')

    # ── per-position table ────────────────────────────────────────────────────
    args.out.mkdir(parents=True, exist_ok=True)
    table_path = args.out / 'reference_region_scores.csv'
    with open(table_path, 'w', newline='') as f:
        wr = csv.writer(f)
        wr.writerow(['norm_distance', 'spread', 'roughness', 'smoothness']
                    + [f'|{a}-{b}|' for a, b in pair_diffs])
        for i in range(args.n_points):
            wr.writerow([f'{x_norm[i]:.4f}', f'{spread[i]:.4f}', f'{roughness[i]:.4f}',
                         f'{smoothness[i]:.4f}']
                        + [f'{pair_diffs[k][i]:.4f}' for k in pair_diffs])
    print(f'\nPer-position scores saved to {table_path}')

    # ── figures ───────────────────────────────────────────────────────────────
    colors = {s: cm.tab10(i / 10) for i, s in enumerate(present)}

    def draw_figure(out_path, window, panel2):
        """3-panel figure: overlay (top), a custom middle panel, pairwise diffs (bottom)."""
        wlo, whi = window

        def shade(ax):
            ax.axvspan(wlo, whi, color='gray', alpha=0.15, zorder=0)

        fig, axes = plt.subplots(3, 1, figsize=(8, 11), sharex=True)

        ax = axes[0]
        for s in present:
            c = colors[s]
            ax.plot(x_norm,  profiles[s], color=c, linewidth=1.8, label=f'{s} (t{used_ids[s]})')
            ax.plot(x_norm, -profiles[s], color=c, linewidth=1.8)
        ax.axhline(0, color='black', linewidth=0.5, linestyle='--')
        shade(ax)
        ax.set_ylim(-args.ylim, args.ylim)
        ax.set_ylabel('R (Å)')
        ax.set_title('Reference channel profiles (shaded = chosen clustering window)')
        ax.legend(loc='upper left', fontsize=8)

        panel2(axes[1], shade)

        ax = axes[2]
        for i, (k, d) in enumerate(pair_diffs.items()):
            ax.plot(x_norm, d, linewidth=1.6, label=f'|{k[0]} − {k[1]}|', color=cm.Dark2(i / 8))
        shade(ax)
        ax.set_ylabel('|ΔR| (Å)')
        ax.set_xlabel('Normalized distance along channel (0 = heme, 1 = exit)')
        ax.set_title('Pairwise differences (confirm the window separates all references)')
        ax.legend(loc='upper left', fontsize=8)

        fig.tight_layout()
        fig.savefig(out_path, dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f'Figure saved to {out_path}')

    # simple figure — spread only
    def panel2_simple(ax, shade):
        ax.plot(x_norm, spread, color='black', linewidth=1.8)
        ax.axvline(x_norm[s_pk], color='red', linewidth=1, linestyle=':')
        ax.axhline(s_th, color='gray', linewidth=0.8, linestyle='--')
        shade(ax)
        ax.set_ylabel(f'spread across refs ({args.spread_metric}, Å)')
        ax.set_title(f'Where the references differ — peak {x_norm[s_pk]:.2f}, '
                     f'window [{swin[0]:.2f}, {swin[1]:.2f}]')

    # smoothness figure — local smoothness (magnitude ignored; wiggle = low)
    def panel2_smoothness(ax, shade):
        ax.plot(x_norm, smoothness, color='black', linewidth=2.2, label='smoothness (1 = smooth)')
        ax.plot(x_norm, _norm(spread), color='tab:blue', linewidth=1.0, linestyle=':',
                alpha=0.7, label='spread (context only, not scored)')
        ax.axvline(x_norm[t_pk], color='red', linewidth=1, linestyle=':')
        ax.axhline(t_th, color='gray', linewidth=0.8, linestyle='--')
        shade(ax)
        ax.set_ylim(-0.02, 1.05)
        ax.set_ylabel('score (0–1)')
        ax.set_title(f'Local smoothness (independent wiggle → low) — peak {x_norm[t_pk]:.2f}, '
                     f'window [{twin[0]:.2f}, {twin[1]:.2f}]')
        ax.legend(loc='lower left', fontsize=8)

    draw_figure(args.out / 'reference_region_simple.png', swin, panel2_simple)
    draw_figure(args.out / 'reference_region.png', twin, panel2_smoothness)


if __name__ == '__main__':
    main()
