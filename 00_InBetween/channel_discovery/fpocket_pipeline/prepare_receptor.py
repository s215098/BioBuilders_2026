#!/usr/bin/env python3
"""Prepare receptor PDB or mmCIF file(s) for AutoDock Vina.

Converts ATOM and HETATM records to PDBQT format: maps elements to AutoDock
atom types and writes zero partial charges (sufficient for Vina's built-in
scoring function, which does not use receptor partial charges).

Handles two input formats automatically:
  - Standard PDB column format (coordinates at columns 31-54)
  - mmCIF space-separated format as produced by RFdiffusion/ProteinMPNN
    (coordinates at space-separated fields 18-20)

Designed for heme-containing enzymes (UPOs, P450s): the heme cofactor and its
iron are preserved and correctly typed (Fe), while crystallographic waters are
removed by default.

Limitations
-----------
- Hydrogens already present are passed through; missing hydrogens are NOT added.
  For best results, run through a protonation tool (reduce, openbabel --addh)
  before using this script.
- Partial charges are set to 0.000. Vina's scoring function does not require
  receptor charges, so this is fine for screening purposes.

Examples
--------
    # Single PDB or CIF file
    python3 prepare_receptor.py -i inputs/8AV5_heme.pdb
    python3 prepare_receptor.py -i inputs/design.cif

    # Batch: all .pdb and .cif files in a directory
    python3 prepare_receptor.py -i inputs/

    # Keep crystallographic waters
    python3 prepare_receptor.py -i inputs/ --keep-water
"""
import argparse
import sys
from pathlib import Path

_ELEMENT_TO_ADTYPE = {
    "C":  "C",
    "N":  "N",
    "O":  "OA",
    "S":  "SA",
    "H":  "HD",
    "P":  "P",
    "F":  "F",
    "CL": "Cl",
    "BR": "Br",
    "I":  "I",
    "FE": "Fe",
    "MG": "Mg",
    "CA": "Ca",
    "MN": "Mn",
    "ZN": "Zn",
    "CU": "Cu",
    "NA": "Na",
    "K":  "K",
}

_WATER_RESNAMES = {"HOH", "WAT", "H2O", "DOD"}


def _guess_element(atom_name: str) -> str:
    """Derive element from atom name when PDB element column is absent."""
    name = atom_name.strip().lstrip("0123456789 ")
    alpha = "".join(c for c in name if c.isalpha()).upper()
    if len(alpha) >= 2 and alpha[:2] in _ELEMENT_TO_ADTYPE:
        return alpha[:2]
    return alpha[:1] if alpha else "C"


def _mmcif_to_pdb_base(fields: list, serial: int):
    """Reconstruct a 66-char PDB-format base line from mmCIF space-separated fields.

    Expected mmCIF _atom_site field order (RFdiffusion/ProteinMPNN output):
      [0]  record (ATOM/HETATM)
      [1]  element
      [2]  atom_name
      [3]  alt_loc  ('.' = none)
      [4]  resname
      [5]  chain
      [7]  resnum
      [18] x
      [19] y
      [20] z

    Returns (base_str, element_str) or (None, None) on parse failure.
    """
    try:
        record    = fields[0]
        element   = fields[1].upper()
        atom_name = fields[2]
        resname   = fields[4][:3]
        chain     = fields[5][:1]
        resnum    = int(fields[7])
        x         = float(fields[18])
        y         = float(fields[19])
        z         = float(fields[20])
    except (IndexError, ValueError):
        return None, None

    # PDB atom name: 1-char elements get a leading space (e.g., " CA ")
    aname = f" {atom_name:<3s}" if len(atom_name) < 4 else f"{atom_name:<4s}"

    base = (
        f"{record:<6s}{serial:5d} {aname} {resname:<3s} {chain}"
        f"{resnum:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.10"
    )
    return base[:66].ljust(66), element


def convert(input_file: Path, output_pdbqt: Path, remove_water: bool = True) -> bool:
    """Convert one PDB or mmCIF file to PDBQT.

    Returns True on success, False on failure.
    """
    lines_out = []
    skipped_water = 0
    skipped_altloc = 0
    serial = 0

    try:
        with open(input_file) as fh:
            for line in fh:
                if not line.startswith(("ATOM", "HETATM")):
                    continue

                fields = line.split()
                is_mmcif = len(fields) > 15

                if is_mmcif:
                    # mmCIF space-separated: altloc is fields[3], '.' means none
                    altloc = fields[3] if len(fields) > 3 else "."
                    if altloc not in (".", " ", "A"):
                        skipped_altloc += 1
                        continue
                    resname = fields[4].upper() if len(fields) > 4 else ""
                    if remove_water and resname in _WATER_RESNAMES:
                        skipped_water += 1
                        continue
                    serial += 1
                    base, element = _mmcif_to_pdb_base(fields, serial)
                    if base is None:
                        continue
                else:
                    # Standard PDB column format
                    altloc = line[16] if len(line) > 16 else " "
                    if altloc not in (" ", "A"):
                        skipped_altloc += 1
                        continue
                    resname = line[17:20].strip().upper()
                    if remove_water and resname in _WATER_RESNAMES:
                        skipped_water += 1
                        continue
                    element = line[76:78].strip().upper() if len(line) > 76 else ""
                    if not element:
                        element = _guess_element(line[12:16])
                    serial += 1
                    base = line[:66].rstrip().ljust(66)

                ad_type = _ELEMENT_TO_ADTYPE.get(element, element.capitalize() or "C")
                lines_out.append(f"{base}    {0.000:6.3f} {ad_type:<2s}\n")

    except OSError as e:
        print(f"Error reading {input_file}: {e}", file=sys.stderr)
        return False

    if not lines_out:
        print(f"Warning: no ATOM/HETATM records written for {input_file.name}",
              file=sys.stderr)
        return False

    try:
        output_pdbqt.write_text("".join(lines_out))
    except OSError as e:
        print(f"Error writing {output_pdbqt}: {e}", file=sys.stderr)
        return False

    notes = []
    if skipped_water:
        notes.append(f"{skipped_water} water(s) removed")
    if skipped_altloc:
        notes.append(f"{skipped_altloc} alt-conformation atom(s) removed")
    note_str = f"  ({', '.join(notes)})" if notes else ""
    print(f"  {input_file.name} -> {output_pdbqt.name}{note_str}")
    return True


def prepare_directory(input_dir: Path, output_dir: Path,
                      remove_water: bool = True, overwrite: bool = False):
    """Batch-convert all .pdb and .cif files in input_dir."""
    files = sorted(input_dir.glob("*.pdb")) + sorted(input_dir.glob("*.cif"))
    if not files:
        print(f"No .pdb or .cif files found in {input_dir}", file=sys.stderr)
        return

    output_dir.mkdir(parents=True, exist_ok=True)
    ok = failed = skipped = 0

    for f in files:
        out = output_dir / f.with_suffix(".pdbqt").name
        if out.exists() and not overwrite:
            print(f"  [SKIP] {f.name}  (PDBQT already exists; use --overwrite to redo)")
            skipped += 1
            continue
        if convert(f, out, remove_water=remove_water):
            ok += 1
        else:
            failed += 1

    print(f"\n{ok} converted, {skipped} skipped, {failed} failed.")


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("-i", "--input", type=Path, required=True,
                   help="PDB/CIF file or directory of PDB/CIF files")
    p.add_argument("-o", "--output", type=Path, default=None,
                   help="Output PDBQT file or directory. "
                        "Defaults to same location as input with .pdbqt extension.")
    p.add_argument("--keep-water", action="store_true",
                   help="Do not remove crystallographic water molecules")
    p.add_argument("--overwrite", action="store_true",
                   help="Re-convert even if the PDBQT already exists")
    args = p.parse_args()

    remove_water = not args.keep_water

    if args.input.is_dir():
        out_dir = args.output if args.output else args.input
        prepare_directory(args.input, out_dir,
                          remove_water=remove_water, overwrite=args.overwrite)

    elif args.input.is_file():
        out = args.output if args.output else args.input.with_suffix(".pdbqt")
        if out.exists() and not args.overwrite:
            print(f"[SKIP] {out} already exists (use --overwrite to redo)")
            sys.exit(0)
        if not convert(args.input, out, remove_water=remove_water):
            sys.exit(1)

    else:
        sys.exit(f"Input not found: {args.input}")


if __name__ == "__main__":
    main()
