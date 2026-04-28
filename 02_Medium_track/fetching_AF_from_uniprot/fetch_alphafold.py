#!/usr/bin/env python3
"""
Fetch AlphaFold PDB structures for a UniProt ID.

Workflow:
  1. Query UniProt to find all AlphaFoldDB cross-references for the entry.
  2. For each AlphaFold entry ID, call the AlphaFold API to get the PDB URL
     and the structure name (uniprotDescription or entryId as fallback).
  3. Download each PDB file, named after the structure.
"""

import sys
import json
import urllib.request
import urllib.error
import re


UNIPROT_API   = "https://rest.uniprot.org/uniprotkb/{uid}.json"
ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction/{entry_id}"


def fetch_json(url):
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req) as resp:
        return json.loads(resp.read())


def safe_filename(name):
    """Strip characters that are invalid in filenames."""
    return re.sub(r'[^\w\-. ]', '_', name).strip()


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


def get_alphafold_metadata(entry_id):
    """Return (pdb_url, name) from the AlphaFold API for a given entry ID."""
    url = ALPHAFOLD_API.format(entry_id=entry_id)
    try:
        results = fetch_json(url)
    except urllib.error.HTTPError as e:
        print(f"  Warning: could not fetch AlphaFold metadata for '{entry_id}': HTTP {e.code}")
        return None, None

    if not results:
        print(f"  Warning: empty response from AlphaFold API for '{entry_id}'")
        return None, None

    entry = results[0]
    pdb_url = entry.get("pdbUrl")
    name = entry.get("uniprotDescription") or entry.get("entryId") or entry_id
    return pdb_url, name


def download_file(url, dest_path):
    urllib.request.urlretrieve(url, dest_path)


def main():
    if len(sys.argv) < 2:
        print("Usage: python fetch_alphafold.py <UniProtID>")
        print("Example: python fetch_alphafold.py P68871")
        sys.exit(1)

    uniprot_id = sys.argv[1].strip().upper()

    print(f"Looking up AlphaFold structures for UniProt ID: {uniprot_id}")
    af_ids = get_alphafold_entry_ids(uniprot_id)
    print(f"Found {len(af_ids)} AlphaFold entry/entries: {', '.join(af_ids)}\n")

    for entry_id in af_ids:
        print(f"Fetching metadata for {entry_id} ...")
        pdb_url, name = get_alphafold_metadata(entry_id)

        if not pdb_url:
            print(f"  Skipping {entry_id} (no PDB URL available).\n")
            continue

        filename = safe_filename(name) + ".pdb"
        print(f"  Name   : {name}")
        print(f"  URL    : {pdb_url}")
        print(f"  Saving : {filename}")

        try:
            download_file(pdb_url, filename)
            print(f"  Done.\n")
        except urllib.error.HTTPError as e:
            print(f"  Error downloading: HTTP {e.code}\n")


if __name__ == "__main__":
    main()