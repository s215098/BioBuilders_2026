#!/usr/bin/env python3
"""
01_fetch_sequences.py
=====================
Step 1 of the EasyTrack pipeline: retrieve protein sequences.

Three modes (set in config under sequences.mode):
  ncbi    — search NCBI protein database via Entrez and download FASTA
  uniprot — search UniProt REST API and download FASTA
  manual  — skip downloading; validate and copy a user-supplied FASTA

Output
------
  <output_dir>/sequences/raw_sequences.fasta
"""

import argparse
import os
import sys
import time
from pathlib import Path

import yaml
import requests

# ── Biopython Entrez (NCBI) ───────────────────────────────────────────────────
try:
    from Bio import Entrez, SeqIO
    BIOPYTHON_OK = True
except ImportError:
    BIOPYTHON_OK = False


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def load_config(config_path: str) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def count_sequences(fasta_path: str) -> int:
    """Return number of sequences in a FASTA file."""
    count = 0
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                count += 1
    return count


# ─────────────────────────────────────────────────────────────────────────────
# Mode 1: NCBI Entrez
# ─────────────────────────────────────────────────────────────────────────────

def fetch_ncbi(query: str, email: str, database: str, max_results: int,
               output_fasta: str):
    """
    Search NCBI protein database and download sequences as FASTA.
    Requires biopython.
    """
    if not BIOPYTHON_OK:
        print("[ERROR] biopython is not installed. Run: pip install biopython")
        sys.exit(1)

    Entrez.email = email
    print(f"[NCBI] Searching '{database}' for: {query}")
    print(f"[NCBI] Max results: {max_results}")

    # Search
    handle = Entrez.esearch(db=database, term=query, retmax=max_results,
                            usehistory="y")
    record = Entrez.read(handle)
    handle.close()

    total = int(record["Count"])
    ids = record["IdList"]
    webenv = record["WebEnv"]
    query_key = record["QueryKey"]

    print(f"[NCBI] Found {total} total hits; downloading {len(ids)} sequences.")

    if len(ids) == 0:
        print("[ERROR] No sequences found. Try a different query.")
        sys.exit(1)

    # Fetch in batches to avoid timeouts
    batch_size = 200
    all_records = []

    for start in range(0, len(ids), batch_size):
        end = min(start + batch_size, len(ids))
        print(f"[NCBI] Fetching sequences {start + 1}–{end} …")
        fetch_handle = Entrez.efetch(
            db=database,
            rettype="fasta",
            retmode="text",
            retstart=start,
            retmax=batch_size,
            webenv=webenv,
            query_key=query_key,
        )
        batch_text = fetch_handle.read()
        fetch_handle.close()
        all_records.append(batch_text)
        time.sleep(0.4)  # be polite to NCBI servers

    with open(output_fasta, "w") as f:
        f.write("\n".join(all_records))

    print(f"[NCBI] Saved {count_sequences(output_fasta)} sequences → {output_fasta}")


# ─────────────────────────────────────────────────────────────────────────────
# Mode 2: UniProt REST API
# ─────────────────────────────────────────────────────────────────────────────

def fetch_uniprot(query: str, reviewed_only: bool, max_results: int,
                  output_fasta: str):
    """
    Query UniProt REST API and download sequences as FASTA.
    No extra dependencies beyond requests.
    """
    # Build filter
    fields = query
    if reviewed_only:
        fields += " AND reviewed:true"

    url = "https://rest.uniprot.org/uniprotkb/stream"
    params = {
        "query": fields,
        "format": "fasta",
        "size": min(max_results, 500),  # UniProt max per page = 500
    }

    print(f"[UniProt] Searching: {fields}")
    print(f"[UniProt] URL: {url}")

    resp = requests.get(url, params=params, stream=True, timeout=60)
    if resp.status_code != 200:
        print(f"[ERROR] UniProt request failed (HTTP {resp.status_code}).")
        print(f"        Response: {resp.text[:300]}")
        sys.exit(1)

    raw = resp.text.strip()
    if not raw:
        print("[ERROR] UniProt returned empty result. Try a different query.")
        sys.exit(1)

    # Limit to max_results sequences
    lines = raw.split("\n")
    kept_lines = []
    seq_count = 0
    for line in lines:
        if line.startswith(">"):
            seq_count += 1
            if seq_count > max_results:
                break
        kept_lines.append(line)

    with open(output_fasta, "w") as f:
        f.write("\n".join(kept_lines))

    n = count_sequences(output_fasta)
    print(f"[UniProt] Saved {n} sequences → {output_fasta}")


# ─────────────────────────────────────────────────────────────────────────────
# Mode 3: Manual
# ─────────────────────────────────────────────────────────────────────────────

def use_manual_fasta(source_path: str, output_fasta: str):
    """
    Validate a user-supplied FASTA and copy it to the output location.
    """
    source = Path(source_path)
    if not source.exists():
        print(f"[ERROR] Manual FASTA not found: {source_path}")
        sys.exit(1)

    # Validate it looks like a FASTA
    with open(source) as f:
        first_line = f.readline().strip()
    if not first_line.startswith(">"):
        print(f"[ERROR] File does not look like a FASTA file: {source_path}")
        sys.exit(1)

    # Copy
    import shutil
    shutil.copy2(source, output_fasta)

    n = count_sequences(output_fasta)
    print(f"[Manual] Copied {n} sequences from {source_path} → {output_fasta}")


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Step 1: Fetch protein sequences from NCBI, UniProt, or a manual FASTA."
    )
    parser.add_argument("--config", required=True,
                        help="Path to pipeline YAML config file")
    args = parser.parse_args()

    cfg = load_config(args.config)

    # Resolve paths (relative to config file's directory)
    config_dir = Path(args.config).parent.parent  # project root
    out_dir = config_dir / cfg["output_dir"] / "sequences"
    out_dir.mkdir(parents=True, exist_ok=True)
    output_fasta = str(out_dir / "raw_sequences.fasta")

    seq_cfg = cfg.get("sequences", {})
    mode = seq_cfg.get("mode", "manual").lower()

    print("=" * 60)
    print(f"  Step 1 — Sequence Retrieval  [mode: {mode}]")
    print("=" * 60)

    if mode == "ncbi":
        email = seq_cfg.get("ncbi_email", "")
        if not email or email == "your.email@example.com":
            print("[WARNING] No NCBI email set in config. Using a placeholder.")
            print("          Set sequences.ncbi_email in your config file.")
            email = "placeholder@example.com"
        fetch_ncbi(
            query=seq_cfg.get("query", ""),
            email=email,
            database=seq_cfg.get("ncbi_database", "protein"),
            max_results=int(seq_cfg.get("max_results", 200)),
            output_fasta=output_fasta,
        )

    elif mode == "uniprot":
        fetch_uniprot(
            query=seq_cfg.get("query", ""),
            reviewed_only=bool(seq_cfg.get("uniprot_reviewed_only", False)),
            max_results=int(seq_cfg.get("max_results", 200)),
            output_fasta=output_fasta,
        )

    elif mode == "manual":
        manual_path = seq_cfg.get("manual_fasta", "")
        if not manual_path:
            print("[ERROR] mode is 'manual' but no sequences.manual_fasta path is set.")
            print("        Edit your config file and set sequences.manual_fasta.")
            sys.exit(1)
        # Resolve relative to config file location
        if not Path(manual_path).is_absolute():
            manual_path = str(config_dir / manual_path)
        use_manual_fasta(manual_path, output_fasta)

    else:
        print(f"[ERROR] Unknown mode '{mode}'. Choose: ncbi, uniprot, or manual.")
        sys.exit(1)

    print(f"\n[✓] Sequences saved to: {output_fasta}")
    print(f"    Next step: python pipeline/02_build_phylogeny.py --config {args.config}")


if __name__ == "__main__":
    main()
