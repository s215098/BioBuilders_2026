#!/usr/bin/env python3
"""Fetch AlphaFold structure for a given UniProt ID."""

import sys
import urllib.request

def fetch_alphafold(uniprot_id, fmt="pdb"):
    ext = "pdb" if fmt == "pdb" else "cif"
    url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.{ext}"
    outfile = f"{uniprot_id}.{ext}"

    print(f"Fetching {url} ...")
    try:
        urllib.request.urlretrieve(url, outfile)
        print(f"Saved to {outfile}")
    except urllib.error.HTTPError as e:
        print(f"Error: {e.code} - could not fetch structure for '{uniprot_id}'")
        print("Check that the UniProt ID is correct and has an AlphaFold entry.")
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python fetch_alphafold.py <UniProtID> [pdb|cif]")
        print("Example: python fetch_alphafold.py P68871")
        sys.exit(1)

    uniprot_id = sys.argv[1].strip().upper()
    fmt = sys.argv[2].lower() if len(sys.argv) > 2 else "pdb"

    if fmt not in ("pdb", "cif"):
        print("Format must be 'pdb' or 'cif'")
        sys.exit(1)

    fetch_alphafold(uniprot_id, fmt)
