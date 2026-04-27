#!/usr/bin/env python3
"""
06b_boltz_vina.py
=================
Optional Step 6b of the EasyTrack pipeline: Boltz-2 AI-guided docking.

This script implements a two-stage approach:
  Stage 1 — Boltz-2 predicts the protein-ligand complex structure from sequence
             and SMILES. This gives a predicted binding pose and affinity score.
  Stage 2 — The predicted ligand position is extracted to define a targeted
             docking box, then AutoDock Vina refines the result with physics.

Why use this instead of (or alongside) step 6?
  - Blind docking (step 6) searches the whole protein surface.
  - Boltz-2 gives a predicted pocket center even when no crystal structure exists.
  - Combining both provides both an AI affinity score and a physics-based
    binding affinity, which together are more informative than either alone.

Core logic is adapted from molecular_docking/src/boltz_vina.py — unchanged
except for config integration and output path management.

Prerequisites:
  pip install boltz          # installs the boltz CLI + model weights (~2 GB)
  boltz requires internet access on first run (MSA server + model download)

Input (from config):
  boltz.sequence  — protein amino acid sequence (single-letter, no spaces)
  boltz.smiles    — SMILES string of the ligand
  ligand SDF      — from step 5 (for the Vina re-docking stage)

Output:
  <output_dir>/boltz/<name>/boltz_input.yaml
  <output_dir>/boltz/<name>/boltz_out/          — raw Boltz outputs
  <output_dir>/boltz/<name>/receptor_clean.pdb  — receptor stripped of AI ligand
  <output_dir>/boltz/<name>/receptor.pdbqt      — receptor ready for Vina
  <output_dir>/boltz/<name>/final_poses.pdbqt   — targeted Vina poses
  <output_dir>/boltz/<name>/boltz_summary.json  — confidence + affinities
"""

import argparse
import glob
import json
import os
import subprocess
import sys
from pathlib import Path
from shutil import which

import yaml


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def load_config(config_path: str) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def log(msg: str):
    print(f"[Boltz] {msg}")


# ─────────────────────────────────────────────────────────────────────────────
# Stage 1a — Run Boltz-2
# (core logic from molecular_docking/src/boltz_vina.py — unchanged)
# ─────────────────────────────────────────────────────────────────────────────

def run_boltz(sequence: str, smiles: str, out_dir: str):
    """
    Write a Boltz YAML input and run boltz predict.
    Returns (structure_path, confidence_score, log10_ic50_or_None).
    """
    log("Generating Boltz YAML input file...")
    yaml_path = f"{out_dir}/boltz_input.yaml"

    yaml_content = f"""sequences:
  - protein:
      id: A
      sequence: {sequence}
  - ligand:
      id: B
      smiles: '{smiles}'
properties:
  - affinity:
      binder: B
"""
    with open(yaml_path, "w") as f:
        f.write(yaml_content.strip())

    log("Running Boltz-2 prediction (MSA server + affinity) — this may take a while...")

    boltz_bin = which("boltz")
    if not boltz_bin:
        print("[ERROR] 'boltz' binary not found on PATH.")
        print("        Install with:  pip install boltz")
        sys.exit(1)

    subprocess.run(
        [boltz_bin, "predict", yaml_path,
         "--out_dir", f"{out_dir}/boltz_out",
         "--use_msa_server",
         "--override"],
        check=True,
        stdout=subprocess.DEVNULL
    )

    boltz_out_dir = f"{out_dir}/boltz_out"

    # Find confidence JSON
    json_files = glob.glob(f"{boltz_out_dir}/**/confidence_*.json", recursive=True)
    if not json_files:
        raise FileNotFoundError(
            f"Boltz finished but no confidence JSON found in {boltz_out_dir}"
        )
    json_path = json_files[0]

    # Find affinity JSON (optional — not all runs produce it)
    affinity_files = glob.glob(f"{boltz_out_dir}/**/affinity_*.json", recursive=True)

    # Find structure file
    cif_files = glob.glob(f"{boltz_out_dir}/**/*.cif", recursive=True)
    pdb_files = glob.glob(f"{boltz_out_dir}/**/*.pdb", recursive=True)

    if cif_files:
        structure_path = cif_files[0]
    elif pdb_files:
        structure_path = pdb_files[0]
    else:
        raise FileNotFoundError(
            f"Boltz finished but no structure file found in {boltz_out_dir}"
        )

    # Read confidence score
    with open(json_path) as f:
        data = json.load(f)
        confidence = data.get("confidence_score", data.get("ptm", 0.0))

    # Read affinity
    log10_ic50 = None
    if affinity_files:
        with open(affinity_files[0]) as f:
            aff_data = json.load(f)
            log10_ic50 = aff_data.get("affinity_pred_value")
            if log10_ic50 is not None:
                log(f"Boltz predicted affinity [log10(IC50)]: {log10_ic50:.2f}")

    log(f"Boltz confidence score: {confidence:.3f}")
    return structure_path, confidence, log10_ic50


# ─────────────────────────────────────────────────────────────────────────────
# Stage 1b — Extract pocket from predicted structure
# (core logic from molecular_docking/src/boltz_vina.py — unchanged)
# ─────────────────────────────────────────────────────────────────────────────

def extract_pocket_and_clean(pdb_file: str, clean_out_path: str):
    """
    Strip the Boltz-placed AI ligand from the predicted structure,
    save receptor-only PDB, and return the pocket center + box size.
    """
    import numpy as np
    from Bio.PDB import PDBParser, PDBIO, Select

    class ReceptorOnlySelect(Select):
        def accept_residue(self, residue):
            if residue.get_resname() in ["UNL", "LIG1", "HOH"]:
                return 0
            return 1

    log("Extracting active site pocket and cleaning AI ligand from structure...")

    if pdb_file.endswith(".cif"):
        from Bio.PDB.MMCIFParser import MMCIFParser
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)

    structure = parser.get_structure("complex", pdb_file)
    coords = []

    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_resname() == "LIG1":
                    for atom in residue:
                        coords.append(atom.get_coord())

    if not coords:
        raise ValueError(
            "Ligand residue 'LIG1' not found in Boltz prediction. "
            "Cannot extract pocket center."
        )

    coords = np.array(coords)
    center = coords.mean(axis=0).tolist()
    raw_span = coords.max(axis=0) - coords.min(axis=0)
    box_size = np.maximum(raw_span + 10.0, [20.0, 20.0, 20.0]).tolist()

    io = PDBIO()
    io.set_structure(structure)
    io.save(clean_out_path, ReceptorOnlySelect())

    log(f"Pocket center: {[round(c, 2) for c in center]}")
    log(f"Box size:      {[round(b, 2) for b in box_size]}")
    log(f"Receptor (no AI ligand) → {clean_out_path}")

    return center, box_size


# ─────────────────────────────────────────────────────────────────────────────
# Stage 1c — Prepare receptor PDBQT from Boltz structure
# (core logic from molecular_docking/src/boltz_vina.py — unchanged)
# ─────────────────────────────────────────────────────────────────────────────

def prep_receptor_meeko(clean_pdb: str, pdbqt_out: str):
    """
    Use Meeko's mk_prepare_receptor.py to prepare the receptor PDBQT.
    This is what boltz_vina.py uses (different from step 4's OpenBabel approach,
    but appropriate here since the Boltz structure is already clean).
    """
    mk_script = which("mk_prepare_receptor.py")
    if not mk_script:
        # Try as a module
        mk_script = None
        for candidate in ["mk_prepare_receptor.py", "mk_prepare_receptor"]:
            if which(candidate):
                mk_script = candidate
                break

    if not mk_script:
        # Fall back to obabel if mk_prepare_receptor is not available
        obabel = which("obabel")
        if obabel:
            log("mk_prepare_receptor.py not found — falling back to OpenBabel.")
            result = subprocess.run(
                [obabel, "-ipdb", clean_pdb, "-opdbqt", "-O", pdbqt_out, "-xr"],
                capture_output=True, text=True
            )
            if result.returncode != 0:
                raise RuntimeError(f"OpenBabel failed: {result.stderr[:200]}")
        else:
            raise RuntimeError(
                "Neither mk_prepare_receptor.py nor obabel found. "
                "Install Meeko fully or conda install -c conda-forge openbabel"
            )
    else:
        log(f"Preparing receptor PDBQT with Meeko ({mk_script})...")
        result = subprocess.run(
            [mk_script, "--read_pdb", clean_pdb, "-p", pdbqt_out],
            capture_output=True, text=True
        )
        if result.returncode != 0:
            raise RuntimeError(f"mk_prepare_receptor failed: {result.stderr[:200]}")

    if not Path(pdbqt_out).exists():
        raise FileNotFoundError(f"Receptor PDBQT not created: {pdbqt_out}")

    log(f"Receptor PDBQT → {pdbqt_out}")


# ─────────────────────────────────────────────────────────────────────────────
# Stage 2 — Targeted Vina docking using Boltz pocket
# (core logic from molecular_docking/src/boltz_vina.py — unchanged)
# ─────────────────────────────────────────────────────────────────────────────

def run_targeted_vina(receptor_pdbqt: str, ligand_sdf: str,
                      center: list, box: list, out_file: str,
                      exhaustiveness: int = 32, n_poses: int = 5):
    """
    Run targeted AutoDock Vina using the pocket center from Boltz.
    Ligand is prepared inline with RDKit + Meeko.
    Returns best binding affinity (kcal/mol).
    """
    from vina import Vina
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from meeko import MoleculePreparation

    log(f"Preparing ligand from {ligand_sdf} (RDKit + Meeko)...")
    suppl = Chem.SDMolSupplier(ligand_sdf)
    mol = next(suppl)
    mol_with_h = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol_with_h, AllChem.ETKDG())
    AllChem.MMFFOptimizeMolecule(mol_with_h)

    preparer = MoleculePreparation()
    preparer.prepare(mol_with_h)
    ligand_pdbqt_string = preparer.write_pdbqt_string()

    log(f"Running targeted Vina at pocket {[round(c, 2) for c in center]}...")
    v = Vina(sf_name='vina')
    v.set_receptor(receptor_pdbqt)
    v.compute_vina_maps(center=center, box_size=box)
    v.set_ligand_from_string(ligand_pdbqt_string)

    v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)
    v.write_poses(out_file, n_poses=n_poses, overwrite=True)

    best_score = v.energies()[0][0]
    log(f"Best Vina affinity: {best_score:.2f} kcal/mol")
    return best_score


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Step 6b: Boltz-2 AI structure prediction + targeted Vina docking."
    )
    parser.add_argument("--config", required=True,
                        help="Path to pipeline YAML config file")
    parser.add_argument("--name", default=None,
                        help="Name tag for this run (default: receptor name from config)")
    # Standalone overrides
    parser.add_argument("--sequence", default=None,
                        help="Protein amino acid sequence (overrides config)")
    parser.add_argument("--smiles", default=None,
                        help="Ligand SMILES string (overrides config)")
    parser.add_argument("--ligand-sdf", default=None,
                        help="Path to ligand SDF (overrides config)")
    args = parser.parse_args()

    cfg = load_config(args.config)
    config_dir = Path(args.config).parent.parent
    out_root = config_dir / cfg["output_dir"]

    boltz_cfg = cfg.get("boltz", {})

    # Resolve sequence
    sequence = args.sequence or boltz_cfg.get("sequence", "")
    if not sequence:
        print("[ERROR] No protein sequence provided.")
        print("        Set boltz.sequence in your config or use --sequence")
        sys.exit(1)
    sequence = sequence.strip().replace(" ", "").replace("\n", "")

    # Resolve SMILES
    smiles = args.smiles or boltz_cfg.get("smiles", "")
    if not smiles:
        print("[ERROR] No ligand SMILES provided.")
        print("        Set boltz.smiles in your config or use --smiles")
        sys.exit(1)

    # Resolve ligand SDF (for the Vina re-docking stage)
    if args.ligand_sdf:
        ligand_sdf = args.ligand_sdf
    else:
        lig_name = cfg.get("ligand", {}).get("name", "ligand").replace(" ", "_")
        ligand_sdf = str(out_root / "ligand" / f"{lig_name}.sdf")

    if not Path(ligand_sdf).exists():
        print(f"[ERROR] Ligand SDF not found: {ligand_sdf}")
        print("        Run step 5 first: python pipeline/05_fetch_ligand.py --config <config>")
        sys.exit(1)

    # Confidence threshold
    threshold = float(boltz_cfg.get("confidence_threshold", 0.60))

    # Vina parameters
    dock_cfg = cfg.get("docking", {})
    exhaustiveness = int(dock_cfg.get("exhaustiveness", 32))
    n_poses = int(dock_cfg.get("n_poses", 5))

    # Name for output directory
    run_name = args.name or cfg.get("receptor", {}).get("pdb_id", "boltz_run")
    boltz_dir = str(out_root / "boltz" / run_name)
    os.makedirs(boltz_dir, exist_ok=True)

    print("=" * 60)
    print("  Step 6b — Boltz-2 + Targeted Vina Docking")
    print("=" * 60)
    print(f"  Sequence length : {len(sequence)} aa")
    print(f"  SMILES          : {smiles[:60]}{'…' if len(smiles) > 60 else ''}")
    print(f"  Ligand SDF      : {ligand_sdf}")
    print(f"  Confidence threshold: {threshold}")
    print(f"  Output dir      : {boltz_dir}")
    print()

    # ── Stage 1: Boltz prediction ─────────────────────────────────────────
    structure_path, confidence, log10_ic50 = run_boltz(sequence, smiles, boltz_dir)

    # Save summary so far
    summary = {
        "run_name": run_name,
        "sequence_length": len(sequence),
        "smiles": smiles,
        "boltz_confidence": confidence,
        "boltz_log10_ic50": log10_ic50,
        "passed_threshold": confidence >= threshold,
        "vina_affinity_kcal_mol": None,
    }

    if confidence < threshold:
        log(f"Confidence {confidence:.3f} is below threshold {threshold}. "
            f"Skipping Vina re-docking.")
        print(f"\n[!] Boltz confidence too low — Vina stage skipped.")
        print(f"    Try lowering boltz.confidence_threshold in your config,")
        print(f"    or check that the sequence and SMILES are correct.")
    else:
        # ── Stage 2: Extract pocket + run targeted Vina ───────────────────
        clean_pdb = str(Path(boltz_dir) / "receptor_clean.pdb")
        receptor_pdbqt = str(Path(boltz_dir) / "receptor.pdbqt")
        final_poses = str(Path(boltz_dir) / "final_poses.pdbqt")

        center, box = extract_pocket_and_clean(structure_path, clean_pdb)
        prep_receptor_meeko(clean_pdb, receptor_pdbqt)

        vina_score = run_targeted_vina(
            receptor_pdbqt, ligand_sdf, center, box, final_poses,
            exhaustiveness=exhaustiveness, n_poses=n_poses
        )
        summary["vina_affinity_kcal_mol"] = vina_score
        summary["pocket_center"] = [round(c, 3) for c in center]
        summary["pocket_box"] = [round(b, 3) for b in box]

    # Save summary JSON
    summary_path = str(Path(boltz_dir) / "boltz_summary.json")
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)

    # ── Final report ──────────────────────────────────────────────────────
    print()
    print("=" * 60)
    print("  Boltz-2 + Vina Results")
    print("=" * 60)
    print(f"  Boltz confidence score  : {confidence:.3f}")
    if log10_ic50 is not None:
        print(f"  Boltz log10(IC50)       : {log10_ic50:.2f}")
    if summary["vina_affinity_kcal_mol"] is not None:
        print(f"  Vina binding affinity   : {summary['vina_affinity_kcal_mol']:.2f} kcal/mol")
        print(f"  Targeted Vina poses     : {final_poses}")
    print(f"  Summary JSON            : {summary_path}")
    print()
    print(f"[✓] Step 6b complete. Output dir: {boltz_dir}")


if __name__ == "__main__":
    main()
