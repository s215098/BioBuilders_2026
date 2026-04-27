#!/usr/bin/env python3
"""
06_docking.py
=============
Step 6 of the EasyTrack pipeline: run AutoDock Vina molecular docking.

Three docking modes (set in config under docking.mode):

  blind     → search box covers the entire receptor surface.
              No prior knowledge needed. Slower.

  targeted  → search box centred on coordinates you provide.
              Use when you have exact coordinates from a crystal structure.

  residues  → search box centred on named active-site residues.
              Specify residue names (e.g. "HEM", "HIS86") or numbers and the
              script calculates the centre automatically from the receptor file.
              Best when the binding site is known from literature but you do
              not have exact Å coordinates.

Ligand preparation:
  The input is a plain SDF file (from step 5 or your own file).
  The script automatically:
    1. Reads the molecule with RDKit
    2. Adds missing hydrogens
    3. Generates a 3D conformation (ETKDGv3) if the SDF has no 3D coords
    4. Energy-minimises with MMFF94
    5. Converts to PDBQT with Meeko
  No manual conversion needed.

Dependencies:
  pip install vina==1.2.7 meeko==0.7.1 rdkit==2025.9.3

Input
-----
  <output_dir>/receptor/<name>_clean_h.pdbqt
  <output_dir>/ligand/<name>.sdf

Output
------
  <output_dir>/docking/<receptor>_<ligand>_<mode>/
      poses.pdbqt       — all docked poses (Vina PDBQT format)
      poses.pdb         — same, PDB format (load in PyMOL)
      docking_log.txt   — binding affinities per pose
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path
from shutil import which

import numpy as np
import yaml

try:
    from vina import Vina
    VINA_OK = True
except ImportError:
    VINA_OK = False

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors
    RDKIT_OK = True
except ImportError:
    RDKIT_OK = False

try:
    from meeko import MoleculePreparation
    MEEKO_OK = True
except ImportError:
    MEEKO_OK = False


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def load_config(config_path: str) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def check_dependencies():
    missing = []
    if not VINA_OK:
        missing.append("vina  (pip install vina==1.2.7)")
    if not RDKIT_OK:
        missing.append("rdkit  (pip install rdkit==2025.9.3)")
    if not MEEKO_OK:
        missing.append("meeko  (pip install meeko==0.7.1)")
    if missing:
        print("[ERROR] Missing required packages:")
        for m in missing:
            print(f"        - {m}")
        sys.exit(1)


# ─────────────────────────────────────────────────────────────────────────────
# Grid box helpers
# ─────────────────────────────────────────────────────────────────────────────

def get_receptor_center_and_size(pdbqt_file: str, buffer: float = 10.0):
    """Bounding box of the whole receptor (blind docking)."""
    coords = _parse_pdbqt_coords(pdbqt_file)
    if not coords:
        print(f"[ERROR] No ATOM/HETATM records found in: {pdbqt_file}")
        sys.exit(1)
    arr = np.array(coords)
    center = np.mean(arr, axis=0)
    size = np.max(arr, axis=0) - np.min(arr, axis=0) + buffer
    return center.tolist(), size.tolist()


def get_center_from_residues(pdbqt_file: str, residue_specs: list,
                              box_size: list) -> tuple:
    """
    Calculate the docking box centre from named active-site residues.

    residue_specs can be a mix of:
      "HEM"      → match any residue named HEM  (HETATM or ATOM)
      "HIS86"    → match residue named HIS at position 86
      "86"       → match any residue at position 86
      86         → same, integer form

    Returns (center, box_size) where box_size is taken from config.
    """
    specs = []
    for s in residue_specs:
        s = str(s).strip().upper()
        # Pure number → match by residue sequence number only
        if s.isdigit():
            specs.append({"type": "number", "value": int(s)})
        else:
            # Try to split name + number (e.g. "HIS86" → name=HIS, num=86)
            import re
            m = re.match(r"^([A-Z]+)(\d+)$", s)
            if m:
                specs.append({"type": "name_and_number",
                               "name": m.group(1), "number": int(m.group(2))})
            else:
                # Just a name (e.g. "HEM", "HIS", "CYS")
                specs.append({"type": "name", "value": s})

    # Parse PDBQT and collect matching atom coordinates
    matched_coords = []
    matched_residues = set()

    with open(pdbqt_file) as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            try:
                res_name = line[17:21].strip().upper()
                res_num  = int(line[22:26].strip())
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except (ValueError, IndexError):
                continue

            for spec in specs:
                hit = False
                if spec["type"] == "number":
                    hit = (res_num == spec["value"])
                elif spec["type"] == "name":
                    hit = (res_name == spec["value"] or
                           res_name.startswith(spec["value"]))
                elif spec["type"] == "name_and_number":
                    hit = (res_name == spec["name"] and
                           res_num == spec["number"])
                if hit:
                    matched_coords.append([x, y, z])
                    matched_residues.add(f"{res_name}{res_num}")
                    break

    if not matched_coords:
        print("[ERROR] No atoms found matching the specified residues:")
        for s in residue_specs:
            print(f"        - {s}")
        print(f"        Check the residue names/numbers in: {pdbqt_file}")
        print("        Tip: open the PDBQT in PyMOL and check residue names.")
        sys.exit(1)

    arr = np.array(matched_coords)
    center = np.mean(arr, axis=0).tolist()

    print(f"[Residues] Matched residues: {sorted(matched_residues)}")
    print(f"[Residues] Centre calculated from {len(matched_coords)} atoms.")

    return center, box_size


def _parse_pdbqt_coords(pdbqt_file: str) -> list:
    coords = []
    with open(pdbqt_file) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append([x, y, z])
                except ValueError:
                    continue
    return coords


# ─────────────────────────────────────────────────────────────────────────────
# Ligand preparation (RDKit + Meeko)
# ─────────────────────────────────────────────────────────────────────────────

def prepare_ligand_pdbqt(sdf_path: str) -> str:
    """
    Prepare ligand for Vina from an SDF file.

    Steps performed automatically:
      1. Read molecule with RDKit
      2. Add hydrogens
      3. Generate 3D conformation with ETKDGv3 (if SDF has no 3D coords)
      4. Energy-minimise with MMFF94
      5. Convert to PDBQT string with Meeko

    Returns PDBQT string ready to pass to Vina.
    """
    print("[Ligand] Reading SDF …")
    suppl = Chem.SDMolSupplier(sdf_path, removeHs=False)
    mol = next(suppl)

    if mol is None:
        suppl2 = Chem.SDMolSupplier(sdf_path, sanitize=False, removeHs=False)
        mol = next(suppl2)
        if mol is None:
            print(f"[ERROR] Could not read molecule from: {sdf_path}")
            sys.exit(1)
        Chem.SanitizeMol(mol)

    # Report molecule info
    formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
    mw = Descriptors.MolWt(mol)
    n_atoms = mol.GetNumAtoms()
    print(f"[Ligand] Molecule: formula={formula}, MW={mw:.1f}, heavy atoms={n_atoms}")

    # Add hydrogens
    mol = Chem.AddHs(mol)
    print(f"[Ligand] Hydrogens added ({mol.GetNumAtoms()} total atoms)")

    # Check for existing 3D conformation
    has_3d = (mol.GetNumConformers() > 0 and
              any(mol.GetConformer().GetAtomPosition(i).x != 0.0
                  for i in range(min(mol.GetNumAtoms(), 5))))

    if not has_3d:
        print("[Ligand] No 3D coordinates — generating conformation with RDKit ETKDGv3 …")
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        result = AllChem.EmbedMolecule(mol, params)
        if result != 0:
            print("[Ligand] ETKDGv3 failed, trying distance geometry fallback …")
            AllChem.EmbedMolecule(mol, AllChem.ETDG())
    else:
        print("[Ligand] 3D coordinates found in SDF — using existing conformation.")

    # Energy minimisation
    ff_result = AllChem.MMFFOptimizeMolecule(mol)
    if ff_result == 0:
        print("[Ligand] MMFF94 energy minimisation: converged.")
    elif ff_result == 1:
        print("[Ligand] MMFF94 energy minimisation: did not fully converge (acceptable).")
    else:
        print("[Ligand] MMFF94 not applicable — skipped.")

    # Meeko PDBQT conversion
    print("[Ligand] Converting to PDBQT with Meeko …")
    prep = MoleculePreparation()
    prep.prepare(mol)
    pdbqt_str = prep.write_pdbqt_string()
    print("[Ligand] ✓ Ligand ready for Vina.")

    return pdbqt_str


# ─────────────────────────────────────────────────────────────────────────────
# Core docking
# ─────────────────────────────────────────────────────────────────────────────

def run_docking(receptor_pdbqt: str, ligand_sdf: str, out_dir: str,
                mode: str, center: list, box_size: list,
                exhaustiveness: int, n_poses: int,
                active_site_residues: list = None):
    """Run AutoDock Vina. Mode can be 'blind', 'targeted', or 'residues'."""
    os.makedirs(out_dir, exist_ok=True)

    log_path    = os.path.join(out_dir, "docking_log.txt")
    poses_pdbqt = os.path.join(out_dir, "poses.pdbqt")
    poses_pdb   = os.path.join(out_dir, "poses.pdb")

    print("-" * 55)
    print(f"  Mode            : {mode.capitalize()} Docking")
    print(f"  Receptor        : {receptor_pdbqt}")
    print(f"  Ligand (SDF)    : {ligand_sdf}")
    print(f"  Exhaustiveness  : {exhaustiveness}")
    print(f"  Poses requested : {n_poses}")
    print("-" * 55)

    # ── Determine search box ──────────────────────────────────────────────
    if mode == "blind":
        print("[Docking] Calculating bounding box for blind docking …")
        center, box_size = get_receptor_center_and_size(receptor_pdbqt)

    elif mode == "residues":
        if not active_site_residues:
            print("[ERROR] mode='residues' requires docking.active_site_residues "
                  "in config.")
            sys.exit(1)
        print(f"[Docking] Calculating box from residues: {active_site_residues}")
        center, box_size = get_center_from_residues(
            receptor_pdbqt, active_site_residues, box_size
        )

    elif mode == "targeted":
        if not center or len(center) != 3:
            print("[ERROR] mode='targeted' requires docking.center: [x, y, z] "
                  "in config.")
            sys.exit(1)

    else:
        print(f"[ERROR] Unknown docking mode '{mode}'. "
              "Choose: blind, targeted, or residues.")
        sys.exit(1)

    print(f"[Docking] Box centre : {[round(c, 2) for c in center]}")
    print(f"[Docking] Box size   : {[round(b, 2) for b in box_size]}")

    # ── Ligand preparation ────────────────────────────────────────────────
    pdbqt_string = prepare_ligand_pdbqt(ligand_sdf)

    # ── Vina setup + docking ──────────────────────────────────────────────
    print(f"[Docking] Initialising Vina …")
    v = Vina(sf_name="vina")
    v.set_receptor(receptor_pdbqt)
    v.compute_vina_maps(center=center, box_size=box_size)
    v.set_ligand_from_string(pdbqt_string)

    print(f"[Docking] Running search (exhaustiveness={exhaustiveness}) …")
    v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)

    # ── Write results ─────────────────────────────────────────────────────
    v.write_poses(poses_pdbqt, n_poses=n_poses, overwrite=True)
    print(f"[Docking] Poses (PDBQT) → {poses_pdbqt}")

    # Optional PDB conversion for PyMOL
    if which("obabel"):
        try:
            subprocess.run(
                ["obabel", "-ipdbqt", poses_pdbqt, "-opdb", "-O", poses_pdb],
                capture_output=True, check=True
            )
            print(f"[Docking] Poses (PDB)   → {poses_pdb}")
        except Exception:
            pass

    # ── Log file ──────────────────────────────────────────────────────────
    scores = v.energies(n_poses=n_poses)
    with open(log_path, "w") as f:
        f.write(f"Receptor         : {receptor_pdbqt}\n")
        f.write(f"Ligand (SDF)     : {ligand_sdf}\n")
        f.write(f"Mode             : {mode}\n")
        f.write(f"Box centre       : {[round(c, 3) for c in center]}\n")
        f.write(f"Box size         : {[round(b, 3) for b in box_size]}\n")
        if active_site_residues and mode == "residues":
            f.write(f"Active site res. : {active_site_residues}\n")
        f.write(f"Exhaustiveness   : {exhaustiveness}\n\n")
        f.write(f"{'Pose':>4}  {'Affinity (kcal/mol)':>20}  {'RMSD lb':>8}  {'RMSD ub':>8}\n")
        f.write("-" * 48 + "\n")
        for i, score in enumerate(scores, 1):
            f.write(f"  {i:2d}  {score[0]:>20.3f}  {score[1]:>8.3f}  {score[2]:>8.3f}\n")

    print(f"[Docking] Log → {log_path}")
    print(f"\n  ✓ Best binding affinity: {scores[0][0]:.3f} kcal/mol  (pose 1)")

    return poses_pdbqt, scores


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Step 6: Run AutoDock Vina docking (blind / targeted / residues)."
    )
    parser.add_argument("--config", required=True,
                        help="Path to pipeline YAML config file")
    parser.add_argument("--receptor", default=None,
                        help="Override receptor PDBQT path")
    parser.add_argument("--ligand", default=None,
                        help="Override ligand SDF path")
    parser.add_argument("--center", type=float, nargs=3, default=None,
                        help="Override box centre X Y Z")
    parser.add_argument("--box-size", type=float, nargs=3, default=None,
                        help="Override box size X Y Z")
    parser.add_argument("--residues", nargs="+", default=None,
                        help="Override active site residues (e.g. HEM HIS86 36)")
    args = parser.parse_args()

    check_dependencies()

    cfg = load_config(args.config)
    config_dir = Path(args.config).parent.parent
    out_root = config_dir / cfg["output_dir"]

    dock_cfg = cfg.get("docking", {})
    rec_cfg  = cfg.get("receptor", {})
    lig_cfg  = cfg.get("ligand", {})

    # Resolve receptor
    if args.receptor:
        receptor_pdbqt = args.receptor
    else:
        name_stem = rec_cfg.get("pdb_id", "").upper() or \
                    Path(rec_cfg.get("local_pdb", "receptor")).stem
        receptor_pdbqt = str(out_root / "receptor" / f"{name_stem}_clean_h.pdbqt")

    # Resolve ligand
    if args.ligand:
        ligand_sdf = args.ligand
    else:
        lig_name = lig_cfg.get("name", "ligand").replace(" ", "_").replace("/", "-")
        ligand_sdf = str(out_root / "ligand" / f"{lig_name}.sdf")

    # Check files
    for fpath, label, hint in [
        (receptor_pdbqt, "Receptor PDBQT",
         "Run step 4: python pipeline/04_prepare_receptor.py --config <config>"),
        (ligand_sdf, "Ligand SDF",
         "Run step 5: python pipeline/05_fetch_ligand.py --config <config>"),
    ]:
        if not Path(fpath).exists():
            print(f"[ERROR] {label} not found: {fpath}")
            print(f"        {hint}")
            sys.exit(1)

    # Docking parameters
    mode              = dock_cfg.get("mode", "blind").lower()
    center            = args.center   or dock_cfg.get("center", [])
    box_size          = args.box_size or dock_cfg.get("box_size", [20.0, 20.0, 20.0])
    exhaustiveness    = int(dock_cfg.get("exhaustiveness", 32))
    n_poses           = int(dock_cfg.get("n_poses", 10))
    active_site_res   = args.residues or dock_cfg.get("active_site_residues", [])

    # Output dir name includes mode
    rec_stem  = Path(receptor_pdbqt).stem.replace("_clean_h", "")
    lig_stem  = Path(ligand_sdf).stem
    run_name  = f"{rec_stem}_{lig_stem}_{mode}"
    out_dir   = str(out_root / "docking" / run_name)

    print("=" * 60)
    print("  Step 6 — Molecular Docking")
    print("=" * 60)

    run_docking(
        receptor_pdbqt=receptor_pdbqt,
        ligand_sdf=ligand_sdf,
        out_dir=out_dir,
        mode=mode,
        center=center,
        box_size=box_size,
        exhaustiveness=exhaustiveness,
        n_poses=n_poses,
        active_site_residues=active_site_res,
    )

    print(f"\n[✓] Docking complete → {out_dir}")
    print(f"    Next: python pipeline/07_parse_results.py --config {args.config}")


if __name__ == "__main__":
    main()
