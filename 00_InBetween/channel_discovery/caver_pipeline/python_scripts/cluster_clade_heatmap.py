"""
Cross-tabulate tunnel-shape clusters against CUPP clades.

Answers: do structures that share a tunnel shape (cluster from cluster_tunnels.py)
also share a CUPP sequence clade?  Joins the two labellings on accession and draws
a heatmap / contingency table (rows = tunnel clusters, columns = CUPP clades, cell =
number of structures with that combination).

Inputs:
  --clusters-csv  cluster_assignments.csv from cluster_tunnels.py (stem, cluster)
  --cupp-file     cupp_clades.txt: tab-separated, accession in the first field and the
                  clade as the last 3 characters of the line (e.g. '2.2')

Run cluster_tunnels.py first so cluster_assignments.csv exists.
"""

import argparse
import csv
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

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


def load_clades(cupp_path):
    """
    Return {accession: clade} from the CUPP file. Accession is the first
    tab-separated field; the clade is the last 3 characters of the line
    (after stripping trailing whitespace / carriage returns).
    """
    clades = {}
    with open(cupp_path) as f:
        for line in f:
            line = line.rstrip('\n').rstrip('\r').rstrip()
            if not line:
                continue
            acc = line.split('\t')[0].strip()
            if acc:
                clades[acc] = line[-3:]
    return clades


def main():
    parser = argparse.ArgumentParser(
        description='Heatmap of tunnel-shape clusters vs CUPP clades.'
    )
    parser.add_argument('--clusters-csv', type=Path,
                        default=ROOT / 'results' / 'clustering' / 'cluster_assignments.csv',
                        help='cluster_assignments.csv from cluster_tunnels.py')
    parser.add_argument('--cupp-file', type=Path,
                        default=ROOT / 'cupp_summary' / 'cupp_clades.txt',
                        help='CUPP clades file (accession + clade as last 3 chars)')
    parser.add_argument('--out', type=Path, default=ROOT / 'cupp_summary',
                        help='Output directory (default: cupp_summary/)')
    parser.add_argument('--cmap', default='Blues',
                        help='Matplotlib colormap for the heatmap (default: Blues)')
    args = parser.parse_args()

    for p in (args.clusters_csv, args.cupp_file):
        if not p.exists():
            sys.exit(f'ERROR: file not found: {p}')

    clusters = load_clusters(args.clusters_csv)
    clades   = load_clades(args.cupp_file)
    print(f'Loaded {len(clusters)} cluster assignments and {len(clades)} CUPP clades')

    # ── join on accession / stem ──────────────────────────────────────────────
    joined = [(clusters[s], clades[s]) for s in clusters if s in clades]
    missing = [s for s in clusters if s not in clades]
    print(f'Joined {len(joined)} structures have both a cluster and a clade'
          + (f' ({len(missing)} clustered structures absent from CUPP: '
             f'{", ".join(missing)})' if missing else ''))
    if not joined:
        sys.exit('No structures overlap between clusters and CUPP clades.')

    row_labels = sorted({c for c, _ in joined})          # tunnel clusters
    col_labels = sorted({cl for _, cl in joined})        # CUPP clades
    row_idx = {c: i for i, c in enumerate(row_labels)}
    col_idx = {cl: j for j, cl in enumerate(col_labels)}

    mat = np.zeros((len(row_labels), len(col_labels)), dtype=int)
    for c, cl in joined:
        mat[row_idx[c], col_idx[cl]] += 1

    # ── print + save contingency table ────────────────────────────────────────
    args.out.mkdir(parents=True, exist_ok=True)
    table_path = args.out / 'cluster_vs_clade.csv'
    with open(table_path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['cluster\\clade'] + col_labels + ['total'])
        for i, c in enumerate(row_labels):
            w.writerow([c] + list(mat[i]) + [int(mat[i].sum())])
        w.writerow(['total'] + [int(mat[:, j].sum()) for j in range(len(col_labels))]
                   + [int(mat.sum())])
    print(f'\nContingency table (rows=cluster, cols=clade):')
    header = 'clu\\cld  ' + '  '.join(f'{cl:>4}' for cl in col_labels)
    print(header)
    for i, c in enumerate(row_labels):
        print(f'  {c:>4}   ' + '  '.join(f'{v:>4}' for v in mat[i]))
    print(f'Table saved to {table_path}')

    # ── heatmap ───────────────────────────────────────────────────────────────
    fig_w = max(4, 1.2 + len(col_labels) * 0.7)
    fig_h = max(3.5, 1.2 + len(row_labels) * 0.6)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(mat, cmap=args.cmap, aspect='auto')

    ax.set_xticks(range(len(col_labels)))
    ax.set_xticklabels(col_labels)
    ax.set_yticks(range(len(row_labels)))
    ax.set_yticklabels([f'Cluster {c}' for c in row_labels])
    ax.set_xlabel('CUPP clade')
    ax.set_ylabel('Tunnel-shape cluster')
    ax.set_title('Tunnel-shape cluster vs CUPP clade '
                 f'(n={len(joined)} structures)')

    # annotate counts; white text on dark cells for contrast
    thresh = mat.max() / 2 if mat.max() else 0
    for i in range(len(row_labels)):
        for j in range(len(col_labels)):
            v = mat[i, j]
            if v:
                ax.text(j, i, str(v), ha='center', va='center',
                        color='white' if v > thresh else 'black', fontsize=9)

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('number of structures')
    fig.tight_layout()
    heatmap_path = args.out / 'cluster_vs_clade_heatmap.png'
    fig.savefig(heatmap_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'Heatmap saved to {heatmap_path}')


if __name__ == '__main__':
    main()
