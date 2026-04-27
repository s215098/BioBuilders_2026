#!/usr/bin/env python3
"""
07_parse_results.py
===================
Step 7 of the EasyTrack pipeline: parse AutoDock Vina results and produce
a summary table.

What this script does:
  1. Finds all docking result folders under <output_dir>/docking/
  2. Parses each poses.pdbqt file to extract binding energies per pose
  3. Writes a CSV summary: one row per pose
  4. Prints a nicely formatted terminal table of the best poses

Input
-----
  <output_dir>/docking/*/poses.pdbqt

Output
------
  <output_dir>/docking/summary.csv
  <output_dir>/docking/best_poses.txt   — plain-text table of top-1 poses
"""

import argparse
import csv
import glob
import os
import sys
from pathlib import Path

import yaml


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def load_config(config_path: str) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


# ─────────────────────────────────────────────────────────────────────────────
# PDBQT parser
# ─────────────────────────────────────────────────────────────────────────────

def parse_vina_pdbqt(pdbqt_path: str) -> list:
    """
    Parse a Vina output PDBQT file.
    Returns a list of dicts, one per pose:
      {
        "pose":       int,
        "affinity":   float (kcal/mol),
        "rmsd_lb":    float,
        "rmsd_ub":    float,
      }
    """
    poses = []
    current_pose = None

    with open(pdbqt_path) as f:
        for line in f:
            line = line.rstrip("\n")

            # New model/pose
            if line.startswith("MODEL"):
                parts = line.split()
                pose_id = int(parts[1]) if len(parts) > 1 else len(poses) + 1
                current_pose = {"pose": pose_id, "affinity": None,
                                "rmsd_lb": None, "rmsd_ub": None}

            # Vina energy line: REMARK VINA RESULT  affinity  rmsd_lb  rmsd_ub
            elif line.startswith("REMARK VINA RESULT"):
                parts = line.split()
                if len(parts) >= 6:
                    current_pose["affinity"] = float(parts[3])
                    current_pose["rmsd_lb"] = float(parts[4])
                    current_pose["rmsd_ub"] = float(parts[5])

            elif line.startswith("ENDMDL"):
                if current_pose and current_pose["affinity"] is not None:
                    poses.append(current_pose)
                current_pose = None

    # Fallback: some Vina versions write scores differently
    if not poses:
        with open(pdbqt_path) as f:
            for line in f:
                if "VINA" in line and "kcal/mol" in line:
                    # Try to extract number from the line
                    import re
                    nums = re.findall(r"[-+]?\d+\.\d+", line)
                    if nums:
                        poses.append({
                            "pose": len(poses) + 1,
                            "affinity": float(nums[0]),
                            "rmsd_lb": float(nums[1]) if len(nums) > 1 else 0.0,
                            "rmsd_ub": float(nums[2]) if len(nums) > 2 else 0.0,
                        })

    return poses


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Step 7: Parse docking results into a summary table."
    )
    parser.add_argument("--config", required=True,
                        help="Path to pipeline YAML config file")
    args = parser.parse_args()

    cfg = load_config(args.config)
    config_dir = Path(args.config).parent.parent
    out_root = config_dir / cfg["output_dir"]
    docking_dir = out_root / "docking"

    print("=" * 60)
    print("  Step 7 — Parse Docking Results")
    print("=" * 60)

    if not docking_dir.exists():
        print(f"[ERROR] Docking output directory not found: {docking_dir}")
        print("        Run step 6 first: python pipeline/06_docking.py --config <config>")
        sys.exit(1)

    # Find all poses.pdbqt files
    pdbqt_files = sorted(docking_dir.glob("*/poses.pdbqt"))

    if not pdbqt_files:
        print(f"[ERROR] No poses.pdbqt files found in: {docking_dir}")
        print("        Expected pattern: <docking_dir>/<run_name>/poses.pdbqt")
        sys.exit(1)

    print(f"[Parse] Found {len(pdbqt_files)} docking run(s).")

    # ── Parse all results ──────────────────────────────────────────────────
    all_rows = []
    best_per_run = []

    for pdbqt_file in pdbqt_files:
        run_name = pdbqt_file.parent.name
        poses = parse_vina_pdbqt(str(pdbqt_file))

        if not poses:
            print(f"[WARNING] No poses parsed from: {pdbqt_file}")
            continue

        print(f"\n  Run: {run_name}")
        print(f"  {'Pose':>4}  {'Affinity (kcal/mol)':>22}  {'RMSD lb':>8}  {'RMSD ub':>8}")
        print(f"  {'─' * 4}  {'─' * 22}  {'─' * 8}  {'─' * 8}")

        for p in poses:
            print(f"  {p['pose']:>4}  {p['affinity']:>22.3f}  {p['rmsd_lb']:>8.3f}  {p['rmsd_ub']:>8.3f}")
            all_rows.append({
                "run": run_name,
                "pose": p["pose"],
                "affinity_kcal_mol": p["affinity"],
                "rmsd_lb": p["rmsd_lb"],
                "rmsd_ub": p["rmsd_ub"],
            })

        best = min(poses, key=lambda x: x["affinity"])
        best_per_run.append({
            "run": run_name,
            "best_affinity_kcal_mol": best["affinity"],
            "best_pose": best["pose"],
        })

    # ── Write CSV ──────────────────────────────────────────────────────────
    summary_csv = str(docking_dir / "summary.csv")
    with open(summary_csv, "w", newline="") as f:
        if all_rows:
            writer = csv.DictWriter(f, fieldnames=all_rows[0].keys())
            writer.writeheader()
            writer.writerows(all_rows)
    print(f"\n[✓] Full summary → {summary_csv}")

    # ── Write best poses table ─────────────────────────────────────────────
    best_txt = str(docking_dir / "best_poses.txt")
    with open(best_txt, "w") as f:
        f.write("Best Docking Results (lowest affinity = strongest predicted binding)\n")
        f.write("=" * 70 + "\n")
        f.write(f"{'Run':<45}  {'Best Affinity':>15}  {'Pose':>5}\n")
        f.write("-" * 70 + "\n")
        for row in sorted(best_per_run, key=lambda x: x["best_affinity_kcal_mol"]):
            f.write(f"{row['run']:<45}  {row['best_affinity_kcal_mol']:>13.3f}  {row['best_pose']:>5}\n")

    print(f"[✓] Best poses table → {best_txt}")

    # ── Terminal summary ───────────────────────────────────────────────────
    if best_per_run:
        print()
        print("  Best binding affinities (kcal/mol, lower = stronger binding):")
        print(f"  {'Run':<45}  {'Affinity':>10}")
        print(f"  {'─' * 45}  {'─' * 10}")
        for row in sorted(best_per_run, key=lambda x: x["best_affinity_kcal_mol"]):
            print(f"  {row['run']:<45}  {row['best_affinity_kcal_mol']:>10.3f}")

    print(f"\n  Pipeline complete!")
    print(f"  Visualise results by loading in PyMOL:")
    for pdbqt_file in pdbqt_files:
        poses_pdb = pdbqt_file.parent / "poses.pdb"
        rec_path = out_root / "receptor"
        print(f"    pymol {poses_pdb}  # + add receptor from {rec_path}")


if __name__ == "__main__":
    main()
