#!/usr/bin/env python3
"""Print the protein residues that line an fpocket pocket.

fpocket writes a 'pocketN_atm.pdb' file listing every protein atom that touches
a given pocket. This script collapses those atoms to unique residues and prints
them, ordered by chain then residue number.

Examples
--------
    python3 pocket_residues.py 8AV5_heme_out/pockets/pocket1_atm.pdb
    python3 pocket_residues.py 8AV5_heme_out/pockets/pocket1_atm.pdb --resname PHE
"""
import argparse
import sys


def read_residues(atm_file):
    """Return a list of unique (resname, chain, resnum) tuples, in file order."""
    seen = {}
    with open(atm_file) as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            resname = line[17:20].strip()
            chain = line[21].strip() or "_"
            resnum = line[22:26].strip()
            key = (chain, resnum, resname)
            seen.setdefault(key, (resname, chain, resnum))
    return list(seen.values())


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("atm", help="fpocket pocket atoms file (e.g. pocket1_atm.pdb)")
    p.add_argument("--resname", action="append", default=None,
                   help="Only show this residue type (repeatable, e.g. "
                        "--resname PHE --resname TYR)")
    p.add_argument("--no-hetatm", action="store_true",
                   help="Exclude non-standard residues such as HEM")
    args = p.parse_args()

    try:
        residues = read_residues(args.atm)
    except FileNotFoundError:
        sys.exit(f"File not found: {args.atm}")
    if not residues:
        sys.exit(f"No ATOM/HETATM records found in {args.atm}")

    HETATM_NONAA = {"HEM", "HEC", "HEB", "HOH", "WAT"}
    if args.resname:
        wanted = {r.upper() for r in args.resname}
        residues = [r for r in residues if r[0].upper() in wanted]
    if args.no_hetatm:
        residues = [r for r in residues if r[0].upper() not in HETATM_NONAA]

    # sort by chain, then numeric residue number
    residues.sort(key=lambda r: (r[1], int(r[2])))

    print(f"Pocket file : {args.atm}")
    print(f"Residues    : {len(residues)}")
    print("-" * 28)
    counts = {}
    for resname, chain, resnum in residues:
        print(f"  {resname} {chain}{resnum}")
        counts[resname] = counts.get(resname, 0) + 1
    print("-" * 28)
    summary = "  ".join(f"{k}:{v}" for k, v in sorted(counts.items()))
    print(f"By type: {summary}")


if __name__ == "__main__":
    main()
