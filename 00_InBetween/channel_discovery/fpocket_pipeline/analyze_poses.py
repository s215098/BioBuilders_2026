#!/usr/bin/env python3
"""Analyze Vina docking poses for feasibility and substrate orientation.

For NNBT (or any ligand with a ring and hydroxyl groups), correct placement in
a UPO active site requires:

  1. Feasibility  — ring centroid within --max-ring-dist Å of the heme Fe
                    AND Vina score ≤ --min-score kcal/mol
  2. Orientation  — ring (atom type 'A') closer to Fe than hydroxyl oxygens
                    (type 'OA'); the oxygen-direction tolerance is --flexibility

Orientation check: d_ring < d_oxygen + flexibility
  flexibility=0  → ring must be strictly closer to Fe than oxygens
  flexibility=2  → oxygens may be up to 2 Å closer than ring and still pass
                   (the default; allows a bit of tilt without rejection)

The receptor PDB/PDBQT is needed only to locate the heme Fe.

Examples
--------
    python3 analyze_poses.py docking_out/8AV5_heme_NNBT_targeted_*.pdbqt \\
        --receptor inputs/8AV5_heme.pdb

    python3 analyze_poses.py results.pdbqt --receptor inputs/8AV5_heme.pdb \\
        --min-score -4 --max-ring-dist 10 --flexibility 3
"""
import argparse
import sys
from pathlib import Path

import numpy as np

# import heme-Fe finder from pocket_box (same directory)
sys.path.insert(0, str(Path(__file__).parent))
from pocket_box import find_heme_fe


def parse_vina_pdbqt(path):
    """Parse multi-model Vina PDBQT. Returns list of pose dicts.

    Each dict: {score: float, atoms: [(x, y, z, atom_type), ...]}
    """
    poses = []
    current = None
    with open(path) as fh:
        for line in fh:
            if line.startswith("MODEL"):
                current = {"score": None, "atoms": []}
            elif current is not None and line.startswith("REMARK VINA RESULT"):
                try:
                    current["score"] = float(line.split()[3])
                except (IndexError, ValueError):
                    pass
            elif current is not None and line.startswith("ATOM"):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    atom_type = line[77:].strip()
                    current["atoms"].append((x, y, z, atom_type))
                except (ValueError, IndexError):
                    pass
            elif current is not None and line.startswith("ENDMDL"):
                poses.append(current)
                current = None
    return poses


def pharmacophore_centroids(atoms):
    """Return (ring_centroid, o_centroid) from a pose atom list.

    ring_centroid : mean position of all 'A'  (aromatic carbon) atoms
    o_centroid    : mean position of all 'OA' (H-bond acceptor oxygen) atoms
    Returns None for either group if that atom type is absent.
    """
    ring = np.array([[x, y, z] for x, y, z, t in atoms if t == "A"])
    oxy  = np.array([[x, y, z] for x, y, z, t in atoms if t == "OA"])
    ring_c = ring.mean(0) if len(ring) else None
    o_c    = oxy.mean(0)  if len(oxy)  else None
    return ring_c, o_c


def check_pose(ring_c, o_c, fe, max_ring_dist, flexibility):
    """Evaluate feasibility and orientation.

    Returns
    -------
    feasible  : ring centroid within max_ring_dist of Fe
    oriented  : d_ring < d_oxy + flexibility  (ring faces heme)
    d_ring    : float (Å)
    d_oxy     : float or None
    """
    d_ring = float(np.linalg.norm(ring_c - fe))
    feasible = d_ring <= max_ring_dist

    if o_c is not None:
        d_oxy = float(np.linalg.norm(o_c - fe))
        oriented = d_ring < d_oxy + flexibility
    else:
        d_oxy = None
        oriented = None  # cannot evaluate without oxygen atoms

    return feasible, oriented, d_ring, d_oxy


def analyze_file(pdbqt_path, fe, min_score, max_ring_dist, flexibility, top_n):
    """Analyze all poses in one PDBQT file. Returns list of result dicts."""
    poses = parse_vina_pdbqt(pdbqt_path)
    if not poses:
        return []

    if top_n:
        poses = poses[:top_n]

    results = []
    for i, pose in enumerate(poses, 1):
        score = pose["score"]
        ring_c, o_c = pharmacophore_centroids(pose["atoms"])

        if ring_c is None:
            results.append({
                "pose": i, "score": score,
                "feasible": None, "oriented": None,
                "d_ring": None, "d_oxy": None,
                "pass": False,
                "note": "no aromatic (A) atoms found — check ligand atom types",
            })
            continue

        feasible, oriented, d_ring, d_oxy = check_pose(
            ring_c, o_c, fe, max_ring_dist, flexibility
        )

        score_ok = (score is not None) and (score <= min_score)
        passed = score_ok and feasible and (oriented is not False)

        results.append({
            "pose": i, "score": score,
            "feasible": feasible, "oriented": oriented,
            "d_ring": d_ring, "d_oxy": d_oxy,
            "pass": passed,
            "note": "",
        })

    return results


def print_results(pdbqt_path, results, min_score, max_ring_dist, flexibility):
    n_pass = sum(r["pass"] for r in results)
    print(f"\n{'='*60}")
    print(f"File      : {pdbqt_path}")
    print(f"Thresholds: score ≤ {min_score}  |  d_ring ≤ {max_ring_dist} Å  "
          f"|  flexibility {flexibility} Å")
    print(f"Poses     : {len(results)}   Passing: {n_pass}")
    print(f"{'─'*60}")
    hdr = f"{'Pose':>5}  {'Score':>7}  {'d_ring':>8}  {'d_oxy':>8}  {'Feasible':>9}  {'Oriented':>9}  Pass"
    print(hdr)
    print(f"{'─'*60}")
    for r in results:
        score  = f"{r['score']:7.3f}" if r["score"]  is not None else "    N/A"
        d_ring = f"{r['d_ring']:8.2f}" if r["d_ring"] is not None else "     N/A"
        d_oxy  = f"{r['d_oxy']:8.2f}"  if r["d_oxy"]  is not None else "     N/A"
        feas   = str(r["feasible"])  if r["feasible"]  is not None else "N/A"
        ori    = str(r["oriented"])  if r["oriented"]  is not None else "N/A"
        flag   = "  PASS" if r["pass"] else ""
        note   = f"  [{r['note']}]" if r["note"] else ""
        print(f"{r['pose']:>5}  {score}  {d_ring}  {d_oxy}  {feas:>9}  {ori:>9}{flag}{note}")
    print(f"{'─'*60}")


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("pdbqt", nargs="+", help="Vina output PDBQT file(s)")
    p.add_argument("--receptor", required=True,
                   help="Receptor PDB/PDBQT containing the heme (used to locate Fe)")

    g = p.add_argument_group("thresholds")
    g.add_argument("--min-score", type=float, default=-5.0,
                   help="Maximum (best) Vina score to accept kcal/mol (default: -5.0)")
    g.add_argument("--max-ring-dist", type=float, default=8.0,
                   help="Max ring-centroid-to-Fe distance in Å (default: 8.0)")
    g.add_argument("--flexibility", type=float, default=2.0,
                   help="Tolerance in Å for oxygen direction: d_ring < d_oxy + flexibility "
                        "(default: 2.0 — oxygens must be at least 2 Å farther from Fe than ring, "
                        "or the whole orientation check is only mildly violated)")
    g.add_argument("--top-n", type=int, default=None,
                   help="Only evaluate the top N poses per file (default: all)")

    args = p.parse_args()

    try:
        fe, ref = find_heme_fe(args.receptor)
    except SystemExit as e:
        sys.exit(str(e))
    print(f"Heme reference : {ref} at {np.round(fe, 2).tolist()} (from {args.receptor})")

    any_pass = False
    for path in args.pdbqt:
        results = analyze_file(
            path, fe, args.min_score, args.max_ring_dist,
            args.flexibility, args.top_n,
        )
        if not results:
            print(f"\n[WARN] no poses parsed from {path}")
            continue
        print_results(path, results, args.min_score, args.max_ring_dist, args.flexibility)
        if any(r["pass"] for r in results):
            any_pass = True

    sys.exit(0 if any_pass else 1)


if __name__ == "__main__":
    main()
