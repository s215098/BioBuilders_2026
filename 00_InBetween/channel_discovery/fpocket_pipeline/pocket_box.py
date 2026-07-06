#!/usr/bin/env python3
"""Compute a Vina docking box (center + size) from an fpocket pocket.

Reads the Voronoi vertices (alpha-sphere centers) of a single fpocket pocket
from its '*_vert.pqr' file, computes the bounding box, and prints the
--center / --box_size arguments ready to paste into docking.py.

Optionally clips the box so it never reaches farther than --max-dist Angstroms
from the heme iron (Fe). This trims the part of the box that points away from
the catalytic center without touching the near-heme side. Because a Vina box is
axis-aligned and symmetric about its center, clipping the far edge necessarily
shifts the center toward the heme -- that movement is reported.

Examples
--------
    python3 pocket_box.py 8AV5_heme_out/pockets/pocket1_vert.pqr
    python3 pocket_box.py 8AV5_heme_out/pockets/pocket1_vert.pqr --buffer 8
    python3 pocket_box.py 8AV5_heme_out/pockets/pocket1_vert.pqr \
        --heme-pdb 8AV5_heme.pdb --max-dist 12
"""
import argparse
import sys
import numpy as np


def read_coords(path, resnames=None, atomname=None):
    """Return an (N, 3) array of coordinates from a PQR/PDB(QT) file.

    If `resnames` is given, only those residue names are kept. If `atomname`
    is given, only that atom name is kept.
    """
    coords = []
    with open(path) as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            if resnames is not None and line[17:20].strip() not in resnames:
                continue
            if atomname is not None and line[12:16].strip() != atomname:
                continue
            coords.append(
                [float(line[30:38]), float(line[38:46]), float(line[46:54])]
            )
    return np.array(coords)


def find_heme_fe(pdb_file):
    """Return (coord, label) for the heme iron.

    Search order:
      1. FE atom in a standard heme residue (HEM / HEC / HEB)
      2. Centroid of all atoms in a standard heme residue
      3. Any FE atom in any HETATM record (catches LIG, HME, etc.)
    """
    fe = read_coords(pdb_file, resnames={"HEM", "HEC", "HEB"}, atomname="FE")
    if len(fe):
        return fe[0], "Fe"
    hem = read_coords(pdb_file, resnames={"HEM", "HEC", "HEB"})
    if len(hem):
        return hem.mean(0), "heme centroid"
    # Fallback: any iron in a HETATM record (e.g. residue named LIG)
    fe_any = read_coords(pdb_file, atomname="FE")
    if len(fe_any):
        return fe_any[0], "Fe (HETATM)"
    raise ValueError(f"No heme iron found in {pdb_file}")


def compute_box(pqr_path, buffer=6.0, min_size=None, max_dist=None, heme_ref=None):
    """Compute Vina box programmatically.

    Parameters
    ----------
    pqr_path  : path to pocket *_vert.pqr
    buffer    : Å added to each axis span (default 6.0)
    min_size  : optional minimum box dimension on each axis
    max_dist  : clip far side to this many Å from heme (requires heme_ref)
    heme_ref  : (coord_array, label) as returned by find_heme_fe(), or None

    Returns
    -------
    center : np.ndarray [x, y, z]
    size   : np.ndarray [x, y, z]
    """
    coords = read_coords(pqr_path)
    if not len(coords):
        raise ValueError(f"No ATOM/HETATM records in {pqr_path}")

    mn, mx = coords.min(0), coords.max(0)
    lo = mn - buffer / 2.0
    hi = mx + buffer / 2.0
    orig_center = (lo + hi) / 2.0

    if max_dist is not None and heme_ref is not None:
        fe, _ = heme_ref
        for i in range(3):
            if fe[i] >= orig_center[i]:
                lo[i] = max(lo[i], fe[i] - max_dist)
            else:
                hi[i] = min(hi[i], fe[i] + max_dist)
        if np.any(hi - lo <= 0):
            raise ValueError(
                f"max_dist={max_dist} collapses the box on some axis. Increase it."
            )

    center = (lo + hi) / 2.0
    size = hi - lo
    if min_size is not None:
        size = np.maximum(size, min_size)

    return center, size


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("pqr", help="fpocket pocket vertices file (e.g. pocket1_vert.pqr)")
    p.add_argument("--buffer", type=float, default=6.0,
                   help="Angstroms added to each axis span so the ligand has room "
                        "to explore (default: 6.0)")
    p.add_argument("--min-size", type=float, default=None,
                   help="Optional floor for each box dimension in Angstroms")
    p.add_argument("--max-dist", type=float, default=None,
                   help="Clip the box so no edge sits farther than this many "
                        "Angstroms from the heme Fe (per axis). Requires --heme-pdb "
                        "or --heme-coord. Only shortens the far-from-heme side.")
    p.add_argument("--heme-pdb",
                   help="PDB/PDBQT containing the heme; the Fe (or heme centroid) "
                        "is used as the reference for --max-dist")
    p.add_argument("--heme-coord", type=float, nargs=3, metavar=("X", "Y", "Z"),
                   help="Heme reference coordinate, instead of --heme-pdb")
    args = p.parse_args()

    coords = read_coords(args.pqr)
    if not len(coords):
        sys.exit(f"No ATOM/HETATM records found in {args.pqr}")

    # Work in edge space: lo/hi per axis, buffer split evenly so the geometric
    # center matches the original (size = span + buffer, center = (min+max)/2).
    mn, mx = coords.min(0), coords.max(0)
    lo = mn - args.buffer / 2.0
    hi = mx + args.buffer / 2.0
    orig_center = (lo + hi) / 2.0

    # --- optional heme-distance clipping --------------------------------------
    fe = None
    if args.max_dist is not None:
        if args.heme_coord is not None:
            fe, ref = np.array(args.heme_coord), "given coord"
        elif args.heme_pdb is not None:
            fe, ref = find_heme_fe(args.heme_pdb)
        else:
            sys.exit("--max-dist requires --heme-pdb or --heme-coord")

        for i in range(3):
            if fe[i] >= orig_center[i]:
                # heme on the high side -> the LOW edge is the far one; pull it in
                lo[i] = max(lo[i], fe[i] - args.max_dist)
            else:
                # heme on the low side -> the HIGH edge is the far one; pull it in
                hi[i] = min(hi[i], fe[i] + args.max_dist)
        if np.any(hi - lo <= 0):
            sys.exit("--max-dist is too small: it collapses the box on some axis. "
                     "Increase --max-dist.")

    center = (lo + hi) / 2.0
    size = hi - lo
    if args.min_size is not None:
        # grow symmetrically about the (possibly clipped) center
        size = np.maximum(size, args.min_size)

    # --- report ---------------------------------------------------------------
    print(f"Pocket file     : {args.pqr}")
    print(f"Alpha spheres   : {len(coords)}")
    print(f"Bounding box min: {np.round(mn, 2).tolist()}")
    print(f"Bounding box max: {np.round(mx, 2).tolist()}")
    print(f"Raw span (A)    : {np.round(mx - mn, 2).tolist()}")
    print(f"Buffer per axis : {args.buffer} A")
    if fe is not None:
        shift = center - orig_center
        print(f"Heme ref ({ref}): {np.round(fe, 2).tolist()}")
        print(f"Max dist to heme: {args.max_dist} A  (far side clipped)")
        print(f"Center shift    : {np.round(shift, 2).tolist()}  "
              f"(|{np.round(np.linalg.norm(shift), 2)}| A toward heme)")
    print("-" * 48)
    cx, cy, cz = np.round(center, 2)
    sx, sy, sz = np.round(size, 1)
    print(f"Center  (X Y Z) : {cx} {cy} {cz}")
    print(f"Box size(X Y Z) : {sx} {sy} {sz}")
    print("-" * 48)
    print("Paste into docking.py:")
    print(f"  --center {cx} {cy} {cz} --box_size {sx} {sy} {sz}")


if __name__ == "__main__":
    main()
