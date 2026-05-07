#!/usr/bin/env python3
"""
run_pipeline.py
===============
BioBuilders 2026 — UPO Structure Retrieval & Classification Pipeline

Usage:
    python run_pipeline.py --config config.yaml [--steps 1,2,3] [--verbose]

Steps:
    1  Structure retrieval   (PDB → AlphaFold DB → flag for Boltz-2)
    2  Long/short UPO classification  (sequence + structural)
    3  Active site QC        (Cys axial, Phe triad, RMSD)

All three steps run by default. Individual steps can be re-run independently
since each reads from the previous step's TSV output.
"""

import argparse
import logging
import sys
import yaml
from pathlib import Path


# ---------------------------------------------------------------------------
# Logging setup
# ---------------------------------------------------------------------------

def setup_logging(verbose: bool, log_dir: Path):
    log_dir.mkdir(parents=True, exist_ok=True)
    level = logging.DEBUG if verbose else logging.INFO
    handlers = [
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(log_dir / "pipeline.log"),
    ]
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=handlers,
    )


# ---------------------------------------------------------------------------
# Step runners
# ---------------------------------------------------------------------------

def run_step1(config: dict):
    """Structure retrieval."""
    from retrieve_structures import run_retrieval, download_references, write_results
    
    output_dir = Path(config["output_dir"])
    struct_dir = Path(config["structure_dir"])

    print("\n" + "="*60)
    print("STEP 1: Structure Retrieval")
    print("="*60)

    # Download reference structures first
    print("\n→ Downloading reference structures (PaDa-I, MroUPO)...")
    ref_paths = download_references(config, struct_dir)

    # Retrieve query structures
    print(f"\n→ Processing {len(config['queries'])} query accessions...")
    results, boltz_todo = run_retrieval(config)

    # Write outputs
    tsv = write_results(results, boltz_todo, output_dir)

    # Summary
    by_source = {}
    for r in results:
        by_source[r.source] = by_source.get(r.source, 0) + 1
    print("\n--- Step 1 Summary ---")
    for src, count in sorted(by_source.items()):
        print(f"  {src}: {count}")
    print(f"  Output: {tsv}")

    return results, ref_paths


def run_step2(config: dict, structure_results=None, ref_paths=None):
    """Long/short classification."""
    from classify_upos import (
        fetch_reference_sequences, classify_by_sequence,
        classify_by_structure, write_classification
    )
    from retrieve_structures import StructureResult
    import csv

    output_dir = Path(config["output_dir"])
    struct_dir = Path(config["structure_dir"])

    print("\n" + "="*60)
    print("STEP 2: Long/Short UPO Classification")
    print("="*60)

    # Load structure results if not passed in
    if structure_results is None:
        tsv_path = output_dir / "structure_sources.tsv"
        if not tsv_path.exists():
            print("ERROR: Run step 1 first (structure_sources.tsv not found)")
            return None
        structure_results = _load_structure_results_from_tsv(tsv_path)

    # Fetch reference sequences
    print("\n→ Fetching reference sequences (AaeUPO, MroUPO) from UniProt...")
    ref_seqs = fetch_reference_sequences(config["references"])

    # Sequence-level classification for all accessions
    print(f"\n→ Sequence-level classification...")
    clf_results = []
    for acc in config["queries"]:
        from classify_upos import fetch_uniprot_sequence
        seq = fetch_uniprot_sequence(acc)
        if seq:
            clf = classify_by_sequence(acc, seq, ref_seqs, config)
        else:
            from classify_upos import ClassificationResult
            clf = ClassificationResult(accession=acc, notes="Sequence fetch failed")
        clf_results.append(clf)

    # Structural classification (TM-align) — only if TMalign is available
    print("\n→ Structural classification (TM-align)...")
    import shutil
    tmalign_available = shutil.which("TMalign") is not None
    if not tmalign_available:
        print("  ⚠ TMalign not found in PATH — skipping structural tier")
        print("    Install: https://zhanggroup.org/TM-align/")
        for clf in clf_results:
            clf.structural_call = "no_tmalign"
    else:
        # Build ref structure path map
        if ref_paths is None:
            ref_paths = {
                "long_upo":  struct_dir / "references" / "5oxu.pdb",
                "short_upo": struct_dir / "references" / "5fuj.pdb",
            }
        # Build structure path map from results
        struct_map = {}
        for r in structure_results:
            if hasattr(r, 'file_path') and r.file_path:
                struct_map[r.accession] = Path(r.file_path)

        for clf in clf_results:
            sp = struct_map.get(clf.accession)
            if sp and sp.exists():
                classify_by_structure(clf, sp, ref_paths, config)

    # Write
    tsv, summary = write_classification(clf_results, output_dir)

    print("\n--- Step 2 Summary ---")
    for clf in clf_results:
        print(f"  {clf.accession}: {clf.final_call} ({clf.confidence})")
    print(f"  Output: {tsv}")
    print(f"  Summary: {summary}")

    return clf_results


def run_step3(config: dict, clf_results=None, ref_paths=None):
    """Active site QC."""
    from active_site_qc import (
        load_structure, run_qc_for_structure, write_qc_results
    )

    output_dir = Path(config["output_dir"])
    struct_dir = Path(config["structure_dir"])

    print("\n" + "="*60)
    print("STEP 3: Active Site QC")
    print("="*60)

    if ref_paths is None:
        ref_paths = {
            "long_upo":  struct_dir / "references" / "5oxu.pdb",
            "short_upo": struct_dir / "references" / "5fuj.pdb",
        }

    # Load structure results
    tsv_path = output_dir / "structure_sources.tsv"
    structure_results = _load_structure_results_from_tsv(tsv_path)

    # Load classification results to pick the right reference per accession
    clf_map = {}
    if clf_results:
        for clf in clf_results:
            clf_map[clf.accession] = clf.final_call

    # Load reference models
    ref_models = {}
    for key, path in ref_paths.items():
        if path.exists():
            model = load_structure(path, key)
            if model:
                ref_models[key] = model

    # QC per structure
    qc_results = []
    phe_triads = {
        "long_upo":  config["references"]["long_upo"].get("phe_triad"),
        "short_upo": config["references"]["short_upo"].get("phe_triad"),
    }

    for r in structure_results:
        if r.source == "BOLTZ_NEEDED":
            print(f"  [{r.accession}] Skipping QC — no structure yet")
            continue

        # Pick reference based on classification
        call = clf_map.get(r.accession, "long")
        ref_key = "long_upo" if call in ("long", "ambiguous") else "short_upo"
        ref_model = ref_models.get(ref_key)

        file_path = Path(r.file_path) if r.file_path else None
        qc = run_qc_for_structure(
            accession=r.accession,
            source=r.source,
            structure_path=file_path,
            ref_model=ref_model,
            ref_phe_triad=phe_triads.get(ref_key),
        )
        qc_results.append(qc)

    tsv = write_qc_results(qc_results, output_dir)

    passed = [q for q in qc_results if q.pass_qc]
    failed = [q for q in qc_results if not q.pass_qc]

    print("\n--- Step 3 Summary ---")
    print(f"  Pass: {len(passed)} | Fail/flag: {len(failed)}")
    for q in failed:
        print(f"  ✗ {q.accession}: {', '.join(q.qc_flags)}")
    print(f"  Output: {tsv}")

    return qc_results


# ---------------------------------------------------------------------------
# Helper: reload TSV results
# ---------------------------------------------------------------------------

def _load_structure_results_from_tsv(tsv_path: Path):
    """Reload StructureResult-like objects from step 1 TSV."""
    from types import SimpleNamespace
    results = []
    with open(tsv_path) as f:
        headers = f.readline().strip().split("\t")
        for line in f:
            vals = line.strip().split("\t")
            r = SimpleNamespace(**dict(zip(headers, vals)))
            results.append(r)
    return results


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="UPO Structure Retrieval & Classification Pipeline"
    )
    parser.add_argument("--config", default="config.yaml",
                        help="Path to config.yaml")
    parser.add_argument("--steps", default="1,2,3",
                        help="Comma-separated steps to run (e.g. 1,2,3)")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    # Load config
    config_path = Path(args.config)
    if not config_path.exists():
        print(f"ERROR: Config file not found: {config_path}")
        sys.exit(1)

    with open(config_path) as f:
        config = yaml.safe_load(f)

    # Setup
    setup_logging(args.verbose, Path(config["logs_dir"]))
    Path(config["output_dir"]).mkdir(parents=True, exist_ok=True)

    steps = [int(s.strip()) for s in args.steps.split(",")]
    print(f"BioBuilders 2026 — UPO Pipeline")
    print(f"Config: {config_path}")
    print(f"Steps:  {steps}")
    print(f"Queries: {len(config['queries'])}")

    # Run steps
    structure_results, ref_paths = None, None
    clf_results = None

    if 1 in steps:
        structure_results, ref_paths = run_step1(config)

    if 2 in steps:
        clf_results = run_step2(config, structure_results, ref_paths)

    if 3 in steps:
        run_step3(config, clf_results, ref_paths)

    print("\n✓ Pipeline complete.")
    print(f"  Results in: {config['output_dir']}/")


if __name__ == "__main__":
    main()
