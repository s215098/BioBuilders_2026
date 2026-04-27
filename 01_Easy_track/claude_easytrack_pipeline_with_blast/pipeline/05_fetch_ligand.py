#!/usr/bin/env python3
"""
05_fetch_ligand.py
==================
Step 5 of the EasyTrack pipeline: obtain the ligand SDF file.

Two modes:
  auto   — search PubChem by name or CID and download the 3D conformer SDF
  manual — copy a user-supplied SDF file

PubChem returns an SDS with a pre-optimised 3D conformer, which is used
directly by the docking step (ligand prep with RDKit + Meeko happens inside
docking.py, not here).

Input
-----
  config: ligand.name, ligand.pubchem_cid, or ligand.local_sdf

Output
------
  <output_dir>/ligand/<name>.sdf
"""

import argparse
import shutil
import sys
import time
from pathlib import Path

import requests
import yaml


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def load_config(config_path: str) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"


# ─────────────────────────────────────────────────────────────────────────────
# PubChem helpers
# ─────────────────────────────────────────────────────────────────────────────

def search_cid_by_name(name: str) -> int:
    """Resolve a compound name to a PubChem CID."""
    url = f"{PUBCHEM_BASE}/compound/name/{requests.utils.quote(name)}/cids/JSON"
    print(f"[PubChem] Searching for: '{name}'")

    resp = requests.get(url, timeout=20)
    if resp.status_code == 404:
        print(f"[ERROR] Compound '{name}' not found on PubChem.")
        print(f"        Try searching manually at: https://pubchem.ncbi.nlm.nih.gov/")
        sys.exit(1)
    if resp.status_code != 200:
        print(f"[ERROR] PubChem name search failed (HTTP {resp.status_code}).")
        sys.exit(1)

    data = resp.json()
    cids = data["IdentifierList"]["CID"]
    cid = cids[0]
    print(f"[PubChem] Found CID: {cid}  (using first result)")
    return cid


def get_compound_info(cid: int) -> dict:
    """Fetch basic compound info (IUPAC name, molecular formula, MW)."""
    url = f"{PUBCHEM_BASE}/compound/cid/{cid}/property/IUPACName,MolecularFormula,MolecularWeight/JSON"
    resp = requests.get(url, timeout=20)
    if resp.status_code != 200:
        return {}
    data = resp.json()
    props = data.get("PropertyTable", {}).get("Properties", [{}])[0]
    return props


def download_3d_sdf(cid: int, output_sdf: str):
    """
    Download the 3D conformer SDF for a PubChem CID.
    Falls back to 2D + embedding message if no 3D conformer exists.
    """
    url_3d = f"{PUBCHEM_BASE}/compound/cid/{cid}/SDF?record_type=3d"
    print(f"[PubChem] Downloading 3D SDF for CID {cid} …")

    resp = requests.get(url_3d, timeout=30)

    if resp.status_code == 404:
        # No pre-computed 3D conformer — fall back to 2D
        print("[WARNING] No pre-computed 3D conformer on PubChem. Downloading 2D SDF.")
        print("          The docking script will generate a 3D conformation via RDKit.")
        url_2d = f"{PUBCHEM_BASE}/compound/cid/{cid}/SDF"
        resp = requests.get(url_2d, timeout=30)
        if resp.status_code != 200:
            print(f"[ERROR] Could not download SDF (HTTP {resp.status_code}).")
            sys.exit(1)
    elif resp.status_code != 200:
        print(f"[ERROR] PubChem SDF download failed (HTTP {resp.status_code}).")
        sys.exit(1)

    with open(output_sdf, "wb") as f:
        f.write(resp.content)

    size = len(resp.content)
    print(f"[PubChem] SDF saved ({size} bytes) → {output_sdf}")


# ─────────────────────────────────────────────────────────────────────────────
# Manual copy
# ─────────────────────────────────────────────────────────────────────────────

def use_local_sdf(source_path: str, output_sdf: str, config_dir: Path):
    """Validate and copy a user-supplied SDF file."""
    src = Path(source_path) if Path(source_path).is_absolute() else config_dir / source_path

    if not src.exists():
        print(f"[ERROR] Local SDF not found: {source_path}")
        print(f"        Resolved path: {src}")
        sys.exit(1)

    # Basic validation
    with open(src) as f:
        content = f.read(200)
    if "$$$$" not in content and "M  END" not in content:
        print(f"[WARNING] File may not be a valid SDF: {source_path}")

    shutil.copy2(str(src), output_sdf)
    print(f"[Manual] Copied SDF → {output_sdf}")


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Step 5: Fetch or copy ligand SDF file."
    )
    parser.add_argument("--config", required=True,
                        help="Path to pipeline YAML config file")
    args = parser.parse_args()

    cfg = load_config(args.config)
    config_dir = Path(args.config).parent.parent
    out_root = config_dir / cfg["output_dir"]

    lig_cfg = cfg.get("ligand", {})
    ligand_name = lig_cfg.get("name", "ligand")
    pubchem_cid = lig_cfg.get("pubchem_cid", None)
    local_sdf = lig_cfg.get("local_sdf", "")

    # Output directory
    lig_dir = out_root / "ligand"
    lig_dir.mkdir(parents=True, exist_ok=True)

    # Clean name for file naming (remove spaces/special chars)
    safe_name = ligand_name.replace(" ", "_").replace("/", "-")
    output_sdf = str(lig_dir / f"{safe_name}.sdf")

    print("=" * 60)
    print("  Step 5 — Ligand Preparation")
    print("=" * 60)

    if local_sdf:
        # Manual mode
        use_local_sdf(local_sdf, output_sdf, config_dir)

    else:
        # Auto mode via PubChem
        if not pubchem_cid:
            pubchem_cid = search_cid_by_name(ligand_name)
            time.sleep(0.3)  # polite delay

        pubchem_cid = int(pubchem_cid)

        # Print compound info
        info = get_compound_info(pubchem_cid)
        if info:
            print(f"[PubChem] Compound info:")
            print(f"          IUPAC name : {info.get('IUPACName', 'N/A')}")
            print(f"          Formula    : {info.get('MolecularFormula', 'N/A')}")
            print(f"          MW         : {info.get('MolecularWeight', 'N/A')} g/mol")

        download_3d_sdf(pubchem_cid, output_sdf)

    print(f"\n[✓] Ligand SDF ready: {output_sdf}")
    print(f"    Next step: python pipeline/06_docking.py --config {args.config}")


if __name__ == "__main__":
    main()
