#!/usr/bin/env python3
"""
04_prepare_receptor.py
======================
Step 4 of the EasyTrack pipeline: download and prepare the receptor protein
for AutoDock Vina.

What this script does:
  1. Download the PDB file from RCSB (or use a local file)
  2. Clean the structure with pdbfixer:
       - Remove solvent (water)
       - Remove heteroatoms (ligands, ions) that are not part of the protein
       - Add missing heavy atoms and residues
       - Add hydrogens at the specified pH
  3. Keep only the requested chains (or all if unspecified)
  4. Convert the cleaned PDB to PDBQT format using OpenBabel
     (required input format for AutoDock Vina)

Dependencies:
  - pdbfixer + openmm  (conda install -c conda-forge pdbfixer openmm)
  - openbabel          (conda install -c conda-forge openbabel)

Manual alternative:
  If you prefer PyMOL for cleaning, follow the instructions printed at the
  end of this script and skip directly to the PDBQT conversion step.

Input
-----
  PDB ID (downloaded automatically) or local PDB file (via config)

Output
------
  <output_dir>/receptor/<PDB_ID>_clean.pdb
  <output_dir>/receptor/<PDB_ID>_clean_h.pdb
  <output_dir>/receptor/<PDB_ID>_clean_h.pdbqt
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path
from shutil import which

import requests
import yaml


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def load_config(config_path: str) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def check_binary(name: str) -> bool:
    return which(name) is not None


# ─────────────────────────────────────────────────────────────────────────────
# Step A: Download PDB
# ─────────────────────────────────────────────────────────────────────────────

def download_pdb(pdb_id: str, output_path: str):
    """Download a PDB file from RCSB PDB."""
    pdb_id = pdb_id.upper()
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    print(f"[PDB] Downloading {pdb_id} from RCSB …")

    resp = requests.get(url, timeout=30)
    if resp.status_code == 404:
        print(f"[ERROR] PDB ID '{pdb_id}' not found on RCSB.")
        print(f"        Try searching at: https://www.rcsb.org/search")
        sys.exit(1)
    if resp.status_code != 200:
        print(f"[ERROR] Failed to download PDB (HTTP {resp.status_code}).")
        sys.exit(1)

    with open(output_path, "wb") as f:
        f.write(resp.content)

    size_kb = len(resp.content) / 1024
    print(f"[PDB] Saved {pdb_id}.pdb ({size_kb:.1f} KB) → {output_path}")


# ─────────────────────────────────────────────────────────────────────────────
# Step B: Clean with pdbfixer
# ─────────────────────────────────────────────────────────────────────────────

def clean_pdb_with_fixer(raw_pdb: str, clean_h_pdb: str,
                          chains_to_keep: list, ph: float):
    """
    Use pdbfixer to:
      - Keep only specified chains (all if chains_to_keep is empty)
      - Remove water molecules and heteroatoms
      - Find + add missing residues and atoms
      - Add hydrogens at the given pH

    Requires: pip install pdbfixer openmm
    """
    try:
        from pdbfixer import PDBFixer
        from openmm.app import PDBFile
    except ImportError:
        print("[ERROR] pdbfixer and openmm are required for automated receptor prep.")
        print("        Install with:  conda install -c conda-forge pdbfixer openmm")
        print_manual_instructions(raw_pdb)
        sys.exit(1)

    print(f"[pdbfixer] Cleaning structure: {raw_pdb}")

    fixer = PDBFixer(filename=raw_pdb)

    # Remove unwanted chains
    if chains_to_keep:
        all_chains = [c.id for c in fixer.topology.chains()]
        chains_to_remove = [c for c in all_chains if c not in chains_to_keep]
        if chains_to_remove:
            print(f"[pdbfixer] Removing chains: {chains_to_remove}")
            fixer.removeChains(chainIds=chains_to_remove)
        else:
            print(f"[pdbfixer] All chains present: {all_chains} — keeping all.")
    else:
        print("[pdbfixer] No chain filter specified — keeping all chains.")

    # Remove water and heteroatoms (HETATM that are not part of the protein)
    fixer.removeHeterogens(keepWater=False)
    print("[pdbfixer] Removed heterogens and water.")

    # Fill in missing structural elements
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    print("[pdbfixer] Added missing residues and atoms.")

    # Protonate at target pH
    fixer.addMissingHydrogens(ph)
    print(f"[pdbfixer] Added hydrogens at pH {ph}.")

    # Write
    with open(clean_h_pdb, "w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)

    print(f"[pdbfixer] Clean + protonated structure → {clean_h_pdb}")


# ─────────────────────────────────────────────────────────────────────────────
# Step C: Convert to PDBQT with OpenBabel
# ─────────────────────────────────────────────────────────────────────────────

def convert_to_pdbqt(clean_h_pdb: str, pdbqt_out: str):
    """
    Convert the protonated PDB to PDBQT using OpenBabel.

    Flags used:
      -xr  : keep only polar H (remove non-polar H — standard for Vina receptors)
      -h   : ensure H atoms are present

    Requires: conda install -c conda-forge openbabel
    """
    obabel = which("obabel")
    if not obabel:
        print("[ERROR] 'obabel' (OpenBabel) not found on PATH.")
        print("        Install with:  conda install -c conda-forge openbabel")
        print_manual_instructions(clean_h_pdb)
        sys.exit(1)

    print(f"[OpenBabel] Converting to PDBQT …")
    cmd = [obabel, "-ipdb", clean_h_pdb, "-opdbqt", "-O", pdbqt_out, "-xr"]
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"[ERROR] OpenBabel conversion failed:")
        print(result.stderr[:500])
        sys.exit(1)

    if not Path(pdbqt_out).exists() or Path(pdbqt_out).stat().st_size == 0:
        print("[ERROR] OpenBabel produced an empty or missing PDBQT file.")
        print(result.stdout[:300])
        sys.exit(1)

    size_kb = Path(pdbqt_out).stat().st_size / 1024
    print(f"[OpenBabel] PDBQT file ({size_kb:.1f} KB) → {pdbqt_out}")


# ─────────────────────────────────────────────────────────────────────────────
# Manual instructions fallback
# ─────────────────────────────────────────────────────────────────────────────

def print_manual_instructions(pdb_path: str):
    print()
    print("=" * 60)
    print("  MANUAL RECEPTOR PREPARATION (PyMOL + OpenBabel)")
    print("=" * 60)
    print(f"  PDB file: {pdb_path}")
    print()
    print("  1. Open PyMOL and load the PDB file:")
    print(f"       load {pdb_path}")
    print()
    print("  2. Clean the structure (remove non-protein atoms):")
    print("       remove chain B          # keep only chain A")
    print("       remove solvent          # remove water")
    print("       remove hetatm           # remove ligands/ions")
    print()
    print("  3. Add hydrogens:")
    print("       h_add                   # adds all hydrogens")
    print()
    print("  4. Save the cleaned structure:")
    pdb_stem = Path(pdb_path).stem
    pdb_dir = Path(pdb_path).parent
    clean_h = pdb_dir / f"{pdb_stem}_clean_h.pdb"
    pdbqt_out = pdb_dir / f"{pdb_stem}_clean_h.pdbqt"
    print(f"       save {clean_h}")
    print()
    print("  5. Convert to PDBQT using OpenBabel:")
    print(f"       obabel -ipdb {clean_h} -opdbqt -O {pdbqt_out} -xr")
    print()
    print("  6. Then re-run the docking step:")
    print(f"       python pipeline/06_docking.py --config <config>")
    print("=" * 60)


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Step 4: Download and prepare receptor for AutoDock Vina."
    )
    parser.add_argument("--config", required=True,
                        help="Path to pipeline YAML config file")
    args = parser.parse_args()

    cfg = load_config(args.config)
    config_dir = Path(args.config).parent.parent
    out_root = config_dir / cfg["output_dir"]

    rec_cfg = cfg.get("receptor", {})
    pdb_id = rec_cfg.get("pdb_id", "").upper()
    local_pdb = rec_cfg.get("local_pdb", "")
    chains_to_keep = rec_cfg.get("chains_to_keep", [])
    ph = float(rec_cfg.get("ph", 7.4))

    if not pdb_id and not local_pdb:
        print("[ERROR] Neither receptor.pdb_id nor receptor.local_pdb is set in config.")
        sys.exit(1)

    # Output directory
    rec_dir = out_root / "receptor"
    rec_dir.mkdir(parents=True, exist_ok=True)

    # Determine name stem for file naming
    name_stem = pdb_id if pdb_id else Path(local_pdb).stem

    raw_pdb = str(rec_dir / f"{name_stem}_raw.pdb")
    clean_h_pdb = str(rec_dir / f"{name_stem}_clean_h.pdb")
    pdbqt_out = str(rec_dir / f"{name_stem}_clean_h.pdbqt")

    print("=" * 60)
    print("  Step 4 — Receptor Preparation")
    print("=" * 60)

    # A — Get the raw PDB
    if local_pdb:
        local_path = Path(local_pdb) if Path(local_pdb).is_absolute() else config_dir / local_pdb
        if not local_path.exists():
            print(f"[ERROR] Local PDB not found: {local_pdb}")
            sys.exit(1)
        import shutil
        shutil.copy2(str(local_path), raw_pdb)
        print(f"[PDB] Using local file: {local_pdb}")
    else:
        download_pdb(pdb_id, raw_pdb)

    # B — Clean with pdbfixer
    clean_pdb_with_fixer(raw_pdb, clean_h_pdb, chains_to_keep, ph)

    # C — Convert to PDBQT
    convert_to_pdbqt(clean_h_pdb, pdbqt_out)

    print(f"\n[✓] Receptor preparation complete.")
    print(f"    Raw PDB          : {raw_pdb}")
    print(f"    Cleaned + H PDB  : {clean_h_pdb}")
    print(f"    PDBQT (for Vina) : {pdbqt_out}")
    print()
    print_manual_instructions(raw_pdb)
    print(f"\n    Next step: python pipeline/05_fetch_ligand.py --config {args.config}")


if __name__ == "__main__":
    main()
