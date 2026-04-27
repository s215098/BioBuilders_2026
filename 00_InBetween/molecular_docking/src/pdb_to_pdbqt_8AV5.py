#!/usr/bin/env python3
"""
Standalone PDB to PDBQT Converter
Safely handles standard amino acids (ATOM) and cofactors like Heme (HETATM).
"""

import argparse
import sys
from pathlib import Path

# Mapping of standard elements to AutoDock Vina specific atom types
_ELEMENT_TO_ADTYPE = {
    "C": "C", "N": "N", "O": "OA", "S": "SA", "H": "HD",
    "P": "P", "F": "F", "CL": "Cl", "BR": "Br", "I": "I",
    "FE": "Fe", "MG": "Mg", "CA": "Ca", "MN": "Mn", "ZN": "Zn",
    "CU": "Cu", "NA": "Na", "K": "K",
}

def convert_pdb_to_pdbqt(input_pdb: Path, output_pdbqt: Path):
    """
    Minimal PDB to PDBQT conversion.
    Retains ATOM and HETATM records (including Heme), assigns a default charge (0.000),
    and maps elements to AutoDock atom types.
    """
    lines = []
    
    try:
        with open(input_pdb, 'r') as f:
            for line in f:
                # Capture both standard residues and heteroatoms (like Heme)
                if line.startswith(("ATOM", "HETATM")):
                    
                    # Extract the element symbol (columns 77-78 in standard PDB format)
                    element = line[76:78].strip().upper() if len(line) > 76 else ""
                    
                    # Get the AutoDock type. If it's empty, default to 'C'.
                    # It specifically maps "FE" to "Fe" so Vina doesn't crash.
                    ad_type = _ELEMENT_TO_ADTYPE.get(element, element.capitalize() if element else "C")
                    
                    # PDBQT strictly requires the standard PDB format up to column 66
                    base = line[:66].rstrip().ljust(66)
                    
                    # Append the dummy charge (0.000) and the AutoDock atom type
                    lines.append(f"{base}    {0.000:6.3f} {ad_type:<2s}\n")
                    
        # Write everything to the output file
        output_pdbqt.write_text("".join(lines))
        print(f"✓ Successfully converted: {input_pdb.name} -> {output_pdbqt.name}")
        
    except Exception as e:
        print(f"Error during conversion: {e}", file=sys.stderr)
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(
        description="Convert a PDB file to PDBQT format, safely handling Heme (HETATM)."
    )
    
    parser.add_argument(
        "-i", "--input", 
        type=Path, 
        required=True, 
        help="Path to the input PDB file"
    )
    parser.add_argument(
        "-o", "--output", 
        type=Path, 
        required=True, 
        help="Path to save the output PDBQT file"
    )
    
    args = parser.parse_args()
    
    if not args.input.exists():
        print(f"Error: Input file '{args.input}' does not exist.", file=sys.stderr)
        sys.exit(1)
        
    convert_pdb_to_pdbqt(args.input, args.output)

if __name__ == "__main__":
    main()
