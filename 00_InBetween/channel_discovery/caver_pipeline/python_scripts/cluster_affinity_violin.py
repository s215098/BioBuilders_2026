"""
Violin plots of predicted affinity across tunnel-shape clusters.

Joins a Boltz-2 summary table (one row per accession, with a ligand-confidence
column and a predicted-affinity column) to the tunnel-shape cluster assignments
written by cluster_tunnels.py, then draws a violin plot of the affinity values
grouped by cluster.

Pipeline:
  1. Read the Boltz summary CSV (--summary-csv).
  2. Keep only rows whose confidence column (--iptm-col, default ligand_iptm)
     is >= the cutoff (--iptm-cutoff, default 0.8).
  3. Join on accession to cluster_assignments.csv (--clusters-csv) produced by
     cluster_tunnels.py.
  4. Draw one violin per cluster of the affinity column (--affinity-col, default
     affinity_pred_value), overlaying the individual candidate points.

Run cluster_tunnels.py first so cluster_assignments.csv exists.
"""

import argparse
import csv
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

ROOT = Path(__file__).parent.parent


def load_clusters(csv_path):
    """Return {stem: cluster_int} from cluster_assignments.csv."""
    clusters = {}
    with open(csv_path, newline='') as f:
        for row in csv.DictReader(f):
            stem = row.get('stem', '').strip()
            try:
                clusters[stem] = int(row['cluster'])
            except (KeyError, ValueError):
                pass
    return clusters


def load_summary(csv_path, acc_col, iptm_col, aff_col, cutoff):
    """
    Read the Boltz summary CSV and return a list of dicts
    {accession, iptm, affinity} for rows that pass the confidence cutoff.
    Also returns (total_rows, dropped_by_cutoff, dropped_missing).
    """
    kept, total, dropped_cut, dropped_missing = [], 0, 0, 0
    with open(csv_path, newline='') as f:
        reader = csv.DictReader(f)
        for col in (acc_col, iptm_col, aff_col):
            if col not in reader.fieldnames:
                sys.exit(f"ERROR: column '{col}' not found in {csv_path}. "
                         f"Available: {', '.join(reader.fieldnames)}")
        for row in reader:
            total += 1
            acc = row.get(acc_col, '').strip()
            try:
                iptm = float(row[iptm_col])
                aff  = float(row[aff_col])
            except (ValueError, TypeError):
                dropped_missing += 1
                continue
            if iptm < cutoff:
                dropped_cut += 1
                continue
            kept.append({'accession': acc, 'iptm': iptm, 'affinity': aff})
    return kept, total, dropped_cut, dropped_missing


def main():
    parser = argparse.ArgumentParser(
        description='Violin plots of predicted affinity across tunnel-shape clusters.'
    )
    parser.add_argument('--summary-csv', type=Path, required=True,
                        help='Boltz summary CSV (e.g. boltz2_summary/boltz_nnbt_summary.csv)')
    parser.add_argument('--clusters-csv', type=Path,
                        default=ROOT / 'results' / 'clustering' / 'cluster_assignments.csv',
                        help='cluster_assignments.csv from cluster_tunnels.py '
                             '(default: results/clustering/cluster_assignments.csv)')
    parser.add_argument('--accession-col', default='accession',
                        help='Column holding the accession / structure stem (default: accession)')
    parser.add_argument('--iptm-col', default='ligand_iptm',
                        help='Confidence column used for filtering (default: ligand_iptm)')
    parser.add_argument('--iptm-cutoff', type=float, default=0.8,
                        help='Keep rows with iptm-col >= this value (default: 0.8)')
    parser.add_argument('--affinity-col', default='affinity_pred_value',
                        help='Affinity column plotted in the violins (default: affinity_pred_value)')
    parser.add_argument('--out', type=Path, default=ROOT / 'boltz2_summary',
                        help='Output directory (default: boltz2_summary/)')
    args = parser.parse_args()

    if not args.summary_csv.exists():
        sys.exit(f'ERROR: summary CSV not found: {args.summary_csv}')
    if not args.clusters_csv.exists():
        sys.exit(f'ERROR: cluster assignments not found: {args.clusters_csv}\n'
                 f'Run cluster_tunnels.py first to generate it.')

    clusters = load_clusters(args.clusters_csv)
    print(f'Loaded {len(clusters)} cluster assignments from {args.clusters_csv.name}')

    kept, total, dropped_cut, dropped_missing = load_summary(
        args.summary_csv, args.accession_col, args.iptm_col, args.affinity_col, args.iptm_cutoff)
    print(f'Summary: {total} rows | dropped {dropped_missing} (missing values), '
          f'{dropped_cut} ({args.iptm_col} < {args.iptm_cutoff}) | {len(kept)} pass cutoff')

    # ── join on accession ─────────────────────────────────────────────────────
    joined = []
    for rec in kept:
        cl = clusters.get(rec['accession'])
        if cl is not None:
            joined.append({**rec, 'cluster': cl})
    print(f'Joined to clusters: {len(joined)} candidates have both a cluster and '
          f'{args.affinity_col}')
    if not joined:
        sys.exit('No accessions overlap between the summary and the cluster assignments.')

    # ── group affinity by cluster ─────────────────────────────────────────────
    by_cluster = {}
    for rec in joined:
        by_cluster.setdefault(rec['cluster'], []).append(rec['affinity'])
    ordered = sorted(by_cluster)

    print('\nAffinity by cluster:')
    for cl in ordered:
        vals = by_cluster[cl]
        print(f'  cluster {cl}: n={len(vals)}  '
              f'mean={np.mean(vals):.3f}  median={np.median(vals):.3f}  '
              f'min={min(vals):.3f}  max={max(vals):.3f}')

    # ── ranked candidate table (written + printed head) ───────────────────────
    args.out.mkdir(parents=True, exist_ok=True)
    ranked = sorted(joined, key=lambda r: r['affinity'])
    rank_path = args.out / 'affinity_by_cluster.csv'
    with open(rank_path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['accession', 'cluster', args.iptm_col, args.affinity_col])
        for r in ranked:
            w.writerow([r['accession'], r['cluster'], f"{r['iptm']:.4f}", f"{r['affinity']:.4f}"])
    print(f'\nRanked candidate table saved to {rank_path} '
          f'(sorted by {args.affinity_col} ascending)')

    # ── violin plot ───────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(max(6, len(ordered) * 1.3), 5))
    positions = list(range(1, len(ordered) + 1))

    # violins only for clusters with >= 2 points (kde needs variance)
    vln_data = [by_cluster[cl] for cl in ordered if len(by_cluster[cl]) >= 2]
    vln_pos  = [p for p, cl in zip(positions, ordered) if len(by_cluster[cl]) >= 2]
    if vln_data:
        parts = ax.violinplot(vln_data, positions=vln_pos, showmeans=False,
                              showmedians=True, showextrema=False)
        for pos, body in zip(vln_pos, parts['bodies']):
            cl = ordered[positions.index(pos)]
            body.set_facecolor(cm.tab10((cl - 1) / 10))
            body.set_edgecolor('black')
            body.set_alpha(0.35)
        parts['cmedians'].set_color('black')

    # overlay individual candidate points (jittered) for every cluster
    rng = np.random.default_rng(0)
    for pos, cl in zip(positions, ordered):
        vals = by_cluster[cl]
        jitter = rng.uniform(-0.08, 0.08, size=len(vals))
        ax.scatter(np.full(len(vals), pos) + jitter, vals,
                   color=cm.tab10((cl - 1) / 10), edgecolor='black',
                   linewidth=0.4, s=30, zorder=3, alpha=0.9)

    ax.set_xticks(positions)
    ax.set_xticklabels([f'Cluster {cl}\n(n={len(by_cluster[cl])})' for cl in ordered])
    ax.set_ylabel(args.affinity_col)
    ax.set_xlabel('Tunnel-shape cluster')
    ax.set_title(f'{args.affinity_col} by tunnel cluster '
                 f'({args.iptm_col} ≥ {args.iptm_cutoff}, n={len(joined)})')
    ax.grid(axis='y', linestyle='--', alpha=0.3)
    fig.tight_layout()
    violin_path = args.out / 'affinity_violin_by_cluster.png'
    fig.savefig(violin_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'Violin plot saved to {violin_path}')


if __name__ == '__main__':
    main()
