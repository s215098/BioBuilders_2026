#!/usr/bin/env python3
"""
02_build_phylogeny.py
=====================
Step 2 of the EasyTrack pipeline: align sequences and build a phylogenetic tree.

Steps performed:
  1. Multiple Sequence Alignment (MSA) using MAFFT
  2. Phylogenetic tree inference using FastTree (ML) or neighbor-joining (fallback)
  3. Tree visualisation saved as PNG

Requires:
  External binaries: mafft, FastTree (or fasttree)
  Python packages: biopython, matplotlib

Input
-----
  <output_dir>/sequences/raw_sequences.fasta

Output
------
  <output_dir>/phylogeny/aligned.fasta
  <output_dir>/phylogeny/tree.nwk
  <output_dir>/phylogeny/tree.png
"""

import argparse
import os
import subprocess
import sys
import tempfile
from io import StringIO
from pathlib import Path
from shutil import which

import yaml

try:
    from Bio import AlignIO, Phylo, SeqIO
    from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
    BIOPYTHON_OK = True
except ImportError:
    BIOPYTHON_OK = False

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    MATPLOTLIB_OK = True
except ImportError:
    MATPLOTLIB_OK = False


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def load_config(config_path: str) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def check_binary(name: str) -> bool:
    return which(name) is not None


def count_seqs(fasta_path: str) -> int:
    n = 0
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                n += 1
    return n


# ─────────────────────────────────────────────────────────────────────────────
# Step A: MAFFT alignment
# ─────────────────────────────────────────────────────────────────────────────

def run_mafft(input_fasta: str, output_fasta: str, mafft_args: str = "--auto"):
    """
    Run MAFFT multiple sequence alignment.

    Falls back to a warning and exits if MAFFT is not installed.
    Install: conda install -c bioconda mafft   OR   brew install mafft
    """
    if not check_binary("mafft"):
        print("[ERROR] 'mafft' binary not found on PATH.")
        print("        Install it with:  conda install -c bioconda mafft")
        print("        Or on macOS:      brew install mafft")
        sys.exit(1)

    n = count_seqs(input_fasta)
    print(f"[MAFFT] Aligning {n} sequences with args: {mafft_args}")

    cmd = ["mafft"] + mafft_args.split() + [input_fasta]
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"[ERROR] MAFFT failed (exit code {result.returncode}):")
        print(result.stderr[:500])
        sys.exit(1)

    aligned = result.stdout
    # Normalise: replace '.' gaps (some aligners) with '-', strip blank lines
    aligned = "\n".join(line for line in aligned.splitlines() if line.strip())
    aligned = aligned.replace(".", "-")

    with open(output_fasta, "w") as f:
        f.write(aligned)

    print(f"[MAFFT] Alignment done → {output_fasta}")
    return output_fasta


# ─────────────────────────────────────────────────────────────────────────────
# Step B: Tree inference
# ─────────────────────────────────────────────────────────────────────────────

def run_fasttree(aligned_fasta: str, output_nwk: str):
    """
    Build a maximum-likelihood tree with FastTree.
    Uses protein model by default (WAG + CAT).

    Install: conda install -c bioconda fasttree
    """
    binary = which("FastTree") or which("fasttree") or which("FastTreeMP")
    if not binary:
        print("[WARNING] FastTree binary not found. Falling back to neighbor-joining.")
        return False

    print(f"[FastTree] Building ML tree with: {binary}")
    cmd = [binary, "-wag", "-quiet", aligned_fasta]

    with open(output_nwk, "w") as f:
        result = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, text=True)

    if result.returncode != 0:
        print(f"[ERROR] FastTree failed (exit code {result.returncode}):")
        print(result.stderr[:500])
        sys.exit(1)

    print(f"[FastTree] Tree saved → {output_nwk}")
    return True


def run_neighbor_joining(aligned_fasta: str, output_nwk: str):
    """
    Fallback tree method using Biopython's neighbor-joining algorithm.
    No external binary required. Slower and less accurate than FastTree.
    """
    if not BIOPYTHON_OK:
        print("[ERROR] biopython required for NJ fallback. pip install biopython")
        sys.exit(1)

    print("[NJ] Building neighbor-joining tree (FastTree not available).")
    print("     Note: FastTree is recommended for larger datasets.")

    try:
        alignment = AlignIO.read(aligned_fasta, "fasta")
    except Exception as e:
        print(f"[ERROR] Could not read alignment: {e}")
        sys.exit(1)

    calculator = DistanceCalculator("blosum62")
    dm = calculator.get_distance(alignment)

    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)

    # Root at midpoint for cleaner visualisation
    tree.root_at_midpoint()

    with open(output_nwk, "w") as f:
        Phylo.write(tree, f, "newick")

    print(f"[NJ] Tree saved → {output_nwk}")
    return True


# ─────────────────────────────────────────────────────────────────────────────
# Step C: Tree visualisation
# ─────────────────────────────────────────────────────────────────────────────

def draw_tree(tree_nwk: str, output_png: str, max_labels: int = 60):
    """
    Draw and save the phylogenetic tree as a PNG.

    Labels are hidden for trees with more than max_labels leaves to avoid clutter.
    Users should still examine the .nwk file in FigTree or iTOL for large trees.
    """
    if not BIOPYTHON_OK:
        print("[WARNING] biopython not available; skipping tree visualisation.")
        return
    if not MATPLOTLIB_OK:
        print("[WARNING] matplotlib not available; skipping tree visualisation.")
        return

    tree = Phylo.read(tree_nwk, "newick")
    n_leaves = len(tree.get_terminals())

    # Scale figure height to number of leaves
    height = max(8, min(n_leaves * 0.25, 80))
    fig = plt.figure(figsize=(16, height))
    axes = fig.add_subplot(1, 1, 1)

    show_labels = n_leaves <= max_labels

    Phylo.draw(
        tree,
        axes=axes,
        do_show=False,
        label_func=(lambda x: x.name if show_labels else ""),
    )

    axes.set_title(
        f"Phylogenetic Tree — {n_leaves} sequences"
        + ("" if show_labels else "\n(Labels hidden: >60 sequences — open tree.nwk in FigTree or iTOL)"),
        fontsize=12
    )
    axes.set_xlabel("Branch length")
    plt.tight_layout()
    plt.savefig(output_png, dpi=150, bbox_inches="tight")
    plt.close()

    print(f"[Tree] Visualisation saved → {output_png}")
    if not show_labels:
        print(f"       Tip: load tree.nwk into https://itol.embl.de or FigTree for interactive labels.")


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Step 2: Multiple Sequence Alignment + Phylogenetic Tree."
    )
    parser.add_argument("--config", required=True,
                        help="Path to pipeline YAML config file")
    parser.add_argument("--input-fasta", default=None,
                        help="Override input FASTA path (default: from config output_dir)")
    args = parser.parse_args()

    cfg = load_config(args.config)
    config_dir = Path(args.config).parent.parent
    out_root = config_dir / cfg["output_dir"]

    # Resolve input
    if args.input_fasta:
        input_fasta = args.input_fasta
    else:
        input_fasta = str(out_root / "sequences" / "raw_sequences.fasta")

    if not Path(input_fasta).exists():
        print(f"[ERROR] Input FASTA not found: {input_fasta}")
        print("        Run step 1 first:  python pipeline/01_fetch_sequences.py --config <config>")
        sys.exit(1)

    # Output directory
    phylo_dir = out_root / "phylogeny"
    phylo_dir.mkdir(parents=True, exist_ok=True)

    aligned_fasta = str(phylo_dir / "aligned.fasta")
    tree_nwk = str(phylo_dir / "tree.nwk")
    tree_png = str(phylo_dir / "tree.png")

    phylo_cfg = cfg.get("phylogeny", {})
    mafft_args = phylo_cfg.get("mafft_args", "--auto")
    tree_method = phylo_cfg.get("tree_method", "fasttree").lower()

    print("=" * 60)
    print("  Step 2 — MSA + Phylogenetic Tree")
    print("=" * 60)

    # A — Alignment
    run_mafft(input_fasta, aligned_fasta, mafft_args)

    # B — Tree
    built = False
    if tree_method == "fasttree":
        built = run_fasttree(aligned_fasta, tree_nwk)
        if not built:
            print("[INFO] Falling back to neighbor-joining...")
            run_neighbor_joining(aligned_fasta, tree_nwk)
    else:
        run_neighbor_joining(aligned_fasta, tree_nwk)

    # C — Visualise
    draw_tree(tree_nwk, tree_png)

    print(f"\n[✓] Phylogeny outputs in: {phylo_dir}")
    print(f"    Aligned FASTA : {aligned_fasta}")
    print(f"    Tree (Newick) : {tree_nwk}")
    print(f"    Tree image    : {tree_png}")
    print(f"\n    Next step: python pipeline/03_select_representatives.py --config {args.config}")


if __name__ == "__main__":
    main()
