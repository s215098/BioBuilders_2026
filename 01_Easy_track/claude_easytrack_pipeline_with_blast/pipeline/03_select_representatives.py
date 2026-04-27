#!/usr/bin/env python3
"""
03_select_representatives.py
=============================
Step 3 of the EasyTrack pipeline: cluster the phylogenetic tree and select
one representative sequence per cluster.

Two modes (set in config under selection.mode):
  auto   → automatically picks the first sequence in each cluster
  manual → generates a cluster report, pauses for user to edit, then filters

How clustering works:
  - Pairwise branch-length distances are computed from the Newick tree
  - Hierarchical clustering (average linkage) groups sequences by distance
  - The cluster_threshold sets how different sequences must be to end up
    in separate clusters (smaller = more clusters = more diverse final set)

Input
-----
  <output_dir>/sequences/raw_sequences.fasta
  <output_dir>/phylogeny/tree.nwk

Output
------
  <output_dir>/phylogeny/clusters.tsv        — all sequences with cluster IDs
  <output_dir>/phylogeny/selections.txt      — editable list of chosen IDs
  <output_dir>/sequences/representatives.fasta
"""

import argparse
import sys
from pathlib import Path

import yaml

try:
    from Bio import Phylo, SeqIO
    BIOPYTHON_OK = True
except ImportError:
    BIOPYTHON_OK = False

try:
    import numpy as np
    from scipy.cluster.hierarchy import linkage, fcluster
    from scipy.spatial.distance import squareform
    SCIPY_OK = True
except ImportError:
    SCIPY_OK = False


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def load_config(config_path: str) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def sanitise_name(name: str) -> str:
    """Biopython sometimes replaces spaces with underscores in Newick names."""
    return name.strip().replace(" ", "_")


# ─────────────────────────────────────────────────────────────────────────────
# Clustering
# ─────────────────────────────────────────────────────────────────────────────

def build_distance_matrix(tree_nwk: str):
    """
    Compute pairwise branch-length distances between all leaf nodes.
    Returns (names_list, numpy distance matrix).
    """
    if not BIOPYTHON_OK:
        print("[ERROR] biopython required. pip install biopython")
        sys.exit(1)
    if not SCIPY_OK:
        print("[ERROR] scipy and numpy required. pip install scipy numpy")
        sys.exit(1)

    tree = Phylo.read(tree_nwk, "newick")
    terminals = tree.get_terminals()
    names = [sanitise_name(t.name) for t in terminals]
    n = len(names)

    print(f"[Cluster] Computing pairwise distances for {n} sequences …")
    dist = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            d = tree.distance(terminals[i], terminals[j])
            dist[i, j] = d
            dist[j, i] = d

    return names, dist


def cluster_sequences(names, dist_matrix, threshold: float):
    """
    Hierarchical clustering (average linkage) on distance matrix.
    Returns dict: {cluster_id: [seq_name, seq_name, ...]}
    """
    condensed = squareform(dist_matrix, checks=False)
    Z = linkage(condensed, method="average")
    cluster_ids = fcluster(Z, threshold, criterion="distance")

    clusters = {}
    for name, cid in zip(names, cluster_ids):
        clusters.setdefault(int(cid), []).append(name)

    return clusters


# ─────────────────────────────────────────────────────────────────────────────
# Writing / Reading selection files
# ─────────────────────────────────────────────────────────────────────────────

def write_clusters_tsv(clusters: dict, output_tsv: str):
    """Write a TSV with columns: cluster_id, rank_in_cluster, sequence_name."""
    with open(output_tsv, "w") as f:
        f.write("cluster_id\trank\tsequence_name\n")
        for cid in sorted(clusters):
            for rank, name in enumerate(clusters[cid], 1):
                f.write(f"{cid}\t{rank}\t{name}\n")
    print(f"[Cluster] Cluster assignments → {output_tsv}")


def write_selections_txt(clusters: dict, output_txt: str):
    """
    Write an editable text file with one representative per cluster (the first member).
    Users can edit this file to change which sequences are kept.
    Lines starting with '#' are comments and are ignored.
    """
    with open(output_txt, "w") as f:
        f.write("# ============================================================\n")
        f.write("# EasyTrack — Representative sequence selection\n")
        f.write("# ============================================================\n")
        f.write("# Each line below is one sequence ID that will be kept.\n")
        f.write("# To change the representative for a cluster:\n")
        f.write("#   - Comment out the current line with '#'\n")
        f.write("#   - Add a new line with the accession you prefer\n")
        f.write("#     (find alternatives in phylogeny/clusters.tsv)\n")
        f.write("# To skip a cluster entirely, comment out its line.\n")
        f.write("# ============================================================\n\n")

        for cid in sorted(clusters):
            members = clusters[cid]
            f.write(f"# Cluster {cid} — {len(members)} member(s)\n")
            # Default: first member as representative
            f.write(f"{members[0]}\n")
            # Show alternatives as comments
            for alt in members[1:]:
                f.write(f"#   {alt}\n")
            f.write("\n")

    print(f"[Cluster] Editable selection file → {output_txt}")


def read_selections(selections_txt: str) -> list:
    """Read the (possibly edited) selections file. Returns list of selected IDs."""
    selected = []
    with open(selections_txt) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                selected.append(line)
    return selected


# ─────────────────────────────────────────────────────────────────────────────
# Filter FASTA
# ─────────────────────────────────────────────────────────────────────────────

def filter_fasta(input_fasta: str, selected_ids: list, output_fasta: str):
    """
    Extract sequences from input_fasta whose IDs appear in selected_ids.
    Handles both plain (>ACCESSION) and pipe-delimited (>db|ACCESSION|name) headers.
    """
    if not BIOPYTHON_OK:
        print("[ERROR] biopython required. pip install biopython")
        sys.exit(1)

    # Build a set of normalised IDs for fast lookup
    selected_set = set(s.strip().replace(" ", "_") for s in selected_ids)

    def extract_id(record_id: str) -> str:
        parts = record_id.split("|")
        if len(parts) >= 2:
            return parts[1]
        return record_id.split()[0]

    kept = []
    with open(input_fasta) as f:
        for record in SeqIO.parse(f, "fasta"):
            rid = sanitise_name(extract_id(record.id))
            if rid in selected_set or sanitise_name(record.id) in selected_set:
                kept.append(record)

    if not kept:
        print(f"[WARNING] No sequences matched the selected IDs.")
        print(f"          Check that IDs in selections.txt match those in {input_fasta}")

    with open(output_fasta, "w") as f:
        SeqIO.write(kept, f, "fasta")

    print(f"[Filter] Wrote {len(kept)} representative sequences → {output_fasta}")
    return len(kept)


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Step 3: Cluster tree and select representative sequences."
    )
    parser.add_argument("--config", required=True,
                        help="Path to pipeline YAML config file")
    parser.add_argument("--apply-selections", action="store_true",
                        help="Re-run only the filtering step using an already-edited "
                             "selections.txt (skip re-clustering)")
    args = parser.parse_args()

    cfg = load_config(args.config)
    config_dir = Path(args.config).parent.parent
    out_root = config_dir / cfg["output_dir"]

    seq_dir = out_root / "sequences"
    phylo_dir = out_root / "phylogeny"

    input_fasta = str(seq_dir / "raw_sequences.fasta")
    tree_nwk = str(phylo_dir / "tree.nwk")
    clusters_tsv = str(phylo_dir / "clusters.tsv")
    selections_txt = str(phylo_dir / "selections.txt")
    output_fasta = str(seq_dir / "representatives.fasta")

    sel_cfg = cfg.get("selection", {})
    threshold = float(sel_cfg.get("cluster_threshold", 0.5))
    mode = sel_cfg.get("mode", "auto").lower()

    print("=" * 60)
    print(f"  Step 3 — Representative Selection  [mode: {mode}]")
    print("=" * 60)

    # ── If --apply-selections flag is passed, skip clustering ─────────────
    if args.apply_selections:
        if not Path(selections_txt).exists():
            print(f"[ERROR] selections.txt not found: {selections_txt}")
            sys.exit(1)
        selected = read_selections(selections_txt)
        print(f"[Apply] Read {len(selected)} selected IDs from {selections_txt}")
        filter_fasta(input_fasta, selected, output_fasta)
        print(f"\n[✓] Representatives saved → {output_fasta}")
        return

    # ── Normal flow: cluster ──────────────────────────────────────────────
    if not Path(tree_nwk).exists():
        print(f"[ERROR] Tree file not found: {tree_nwk}")
        print("        Run step 2 first:  python pipeline/02_build_phylogeny.py --config <config>")
        sys.exit(1)

    names, dist_matrix = build_distance_matrix(tree_nwk)
    clusters = cluster_sequences(names, dist_matrix, threshold)

    n_clusters = len(clusters)
    n_seqs = len(names)
    print(f"[Cluster] {n_seqs} sequences → {n_clusters} clusters "
          f"(threshold: {threshold})")

    write_clusters_tsv(clusters, clusters_tsv)
    write_selections_txt(clusters, selections_txt)

    # ── Auto mode: just filter immediately ───────────────────────────────
    if mode == "auto":
        selected = read_selections(selections_txt)
        print(f"[Auto] Using first member of each cluster ({len(selected)} sequences).")
        filter_fasta(input_fasta, selected, output_fasta)

    # ── Manual mode: pause for user review ───────────────────────────────
    elif mode == "manual":
        print("\n" + "=" * 60)
        print("  MANUAL REVIEW REQUIRED")
        print("=" * 60)
        print(f"\n  A cluster report has been saved:")
        print(f"    Cluster assignments : {clusters_tsv}")
        print(f"    Selection file      : {selections_txt}")
        print()
        print("  Please open the selection file in a text editor:")
        print(f"    nano {selections_txt}")
        print(f"    # or: code {selections_txt}")
        print()
        print("  The file lists one representative per cluster.")
        print("  Edit it to keep/change/remove any representatives.")
        print("  Also review the tree image:")
        print(f"    {str(phylo_dir / 'tree.png')}")
        print()
        print("  When you're done editing, re-run with:")
        print(f"    python pipeline/03_select_representatives.py --config {args.config} --apply-selections")
        print()
        print("  Or to continue automatically with the defaults (no editing):")
        print(f"    python pipeline/03_select_representatives.py --config {args.config} --apply-selections")
        print("=" * 60)
        return   # stop here — user must re-run

    print(f"\n[✓] Representative sequences → {output_fasta}")
    print(f"    Next step: python pipeline/04_prepare_receptor.py --config {args.config}")


if __name__ == "__main__":
    main()
