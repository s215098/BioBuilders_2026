#!/usr/bin/env python3
"""
Fetch AlphaFold PDB structures for a UniProt ID.
Workflow:
  1. Query UniProt to find all AlphaFoldDB cross-references for the entry.
  2. For each AlphaFold entry ID, call the AlphaFold API to get the PDB URL,
     structure name, mean pLDDT (globalMetricValue), and number of chains.
  3. Optionally filter to monomers only (--monomer-only).
  4. Among any alternatives for the same fragment, pick the highest mean pLDDT.
  5. Download each PDB file into fetched_structures/, named:
       {description}_AF-{entry_id}_{model_version}.pdb
       e.g. Hemoglobin_subunit_beta_AF-P68871-F1_v4.pdb
  6. Write a manifest TSV summarising all downloads.
"""

import sys
import json
import re
import os
import urllib.request
import urllib.error

UNIPROT_API   = "https://rest.uniprot.org/uniprotkb/{uid}.json"
ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction/{entry_id}"
OUTPUT_DIR    = "fetched_structures"


def fetch_json(url):
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req) as resp:
        return json.loads(resp.read())


def safe_filename(name):
    """Replace spaces with underscores; strip characters invalid in filenames."""
    name = name.replace(" ", "_")
    return re.sub(r'[^\w\-.]', '_', name).strip("_")


def extract_model_version(pdb_url):
    """Pull the version string (e.g. 'v4') out of a pdbUrl like .../AF-X-F1-model_v4.pdb"""
    match = re.search(r'model_(v\d+)\.pdb', pdb_url)
    return match.group(1) if match else "vX"


def get_alphafold_entry_ids(uniprot_id):
    """Return list of AlphaFold entry IDs from UniProt cross-references."""
    url = UNIPROT_API.format(uid=uniprot_id)
    try:
        data = fetch_json(url)
    except urllib.error.HTTPError as e:
        print(f"Error fetching UniProt entry '{uniprot_id}': HTTP {e.code}")
        sys.exit(1)
    refs = data.get("uniProtKBCrossReferences", [])
    af_ids = [r["id"] for r in refs if r.get("database") == "AlphaFoldDB"]
    if not af_ids:
        print(f"No AlphaFoldDB entries found for '{uniprot_id}' in UniProt.")
        sys.exit(0)
    return af_ids


def get_best_alphafold_entry(entry_id, monomer_only=False):
    """
    Fetch all predictions for entry_id from the AlphaFold API.
    Filter by monomer (numChains == 1) if requested.
    Return the candidate with the highest globalMetricValue (mean pLDDT),
    as a dict with keys: pdb_url, name, confidence, num_chains, version.
    Returns None if nothing passes the filters.
    """
    url = ALPHAFOLD_API.format(entry_id=entry_id)
    try:
        results = fetch_json(url)
    except urllib.error.HTTPError as e:
        print(f"  Warning: could not fetch metadata for '{entry_id}': HTTP {e.code}")
        return None

    if not results:
        print(f"  Warning: empty response from AlphaFold API for '{entry_id}'")
        return None

    candidates = []
    for entry in results:
        num_chains = entry.get("numChains", 1)  # default 1 for legacy monomer entries
        if monomer_only and num_chains != 1:
            continue
        pdb_url = entry.get("pdbUrl")
        if not pdb_url:
            continue
        # Support both old (entryId) and new (modelEntityId) API field names
        name = (entry.get("uniprotDescription")
                or entry.get("modelEntityId")
                or entry.get("entryId")
                or entry_id)
        confidence = entry.get("globalMetricValue")  # mean pLDDT, 0–100
        version = extract_model_version(pdb_url)
        candidates.append({
            "pdb_url":    pdb_url,
            "name":       name,
            "confidence": confidence if confidence is not None else 0.0,
            "num_chains": num_chains,
            "version":    version,
        })

    if not candidates:
        reason = "no monomer entries found" if monomer_only else "no usable entries found"
        print(f"  Skipping {entry_id}: {reason}.")
        return None

    # Pick the entry with the highest mean pLDDT
    best = max(candidates, key=lambda c: c["confidence"])
    return best


def download_file(url, dest_path):
    urllib.request.urlretrieve(url, dest_path)


def main():
    import argparse

    # set up arguments:
    parser = argparse.ArgumentParser(
        description="Download AlphaFold PDB structures for a UniProt ID."
    )
    parser.add_argument("uniprot_id",
                        help="UniProt accession, e.g. P68871"
    )
    parser.add_argument(
        "--monomer-only",
        action="store_true",
        help="Skip any multi-chain (complex/homodimer) entries"
    )
    args = parser.parse_args()

    # and get the uniprot id and if searching for only monomers:
    uniprot_id = args.uniprot_id.strip().upper()
    print(f"Looking up AlphaFold structures for UniProt ID: {uniprot_id}")
    if args.monomer_only:
        print("  (monomer-only mode: skipping multi-chain entries)")

    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    af_ids = get_alphafold_entry_ids(uniprot_id)
    print(f"Found {len(af_ids)} AlphaFold entry/entries: {', '.join(af_ids)}\n")

    for entry_id in af_ids:
        print(f"Fetching metadata for {entry_id} ...")
        best = get_best_alphafold_entry(entry_id, monomer_only=args.monomer_only)
        if not best:
            print()
            continue

        # e.g. Hemoglobin_subunit_beta_AF-P68871-F1_v4.pdb
        filename = f"{safe_filename(best['name'])}_AF-{entry_id}_{best['version']}.pdb"
        dest_path = os.path.join(OUTPUT_DIR, filename)
        entry_url = f"https://alphafold.ebi.ac.uk/entry/AF-{entry_id}"

        # print results:
        print(f"  Name       : {best['name']}")
        print(f"  Confidence : {best['confidence']:.1f} mean pLDDT")
        print(f"  Chains     : {best['num_chains']}")
        print(f"  Model ver  : {best['version']}")
        print(f"  Entry URL  : {entry_url}")
        print(f"  Saving     : {dest_path}")

        # And download the file:
        try:
            download_file(best["pdb_url"], dest_path)
            print(f"  Done.\n")
        except urllib.error.HTTPError as e:
            print(f"  Error downloading: HTTP {e.code}\n")

if __name__ == "__main__":
    main()