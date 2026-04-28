import argparse
from ast import arg
import os
import numpy as np
from pathlib import Path
from vina import Vina
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation

# def get_receptor_center_and_size(pdb_file):
#     """Calculates the center and bounding box size of the receptor for blind docking."""
#     coordinates = []
#     with open(pdb_file, 'r') as f:
#         for line in f:
#             if line.startswith("ATOM") or line.startswith("HETATM"):
#                 if line[16] == 'A' or line[16] == ' ' :
#                     x = float(line[30:38])
#                     y = float(line[38:46])
#                     z = float(line[46:54])
#                     coordinates.append([x, y, z])
    
#     coords = np.array(coordinates)
#     center = np.mean(coords, axis=0)
    
#     # (Max - Min) per axis gives the full span of the protein
#     size = np.max(coords, axis=0) - np.min(coords, axis=0)
    
#     # Add a 10 Angstrom buffer
#     size = size + 10.0
    
#     return center.tolist(), size.tolist()

def get_receptor_center_and_size(pdbqt_file):
    """Calculates the center and bounding box size of the receptor for blind docking."""
    coordinates = []
    with open(pdbqt_file, 'r') as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coordinates.append([x, y, z])
                except ValueError:
                    continue  # skip malformed lines
    
    coords = np.array(coordinates)
    center = np.mean(coords, axis=0)
    size = np.max(coords, axis=0) - np.min(coords, axis=0) + 10.0
    
    return center.tolist(), size.tolist()

def main():
    parser = argparse.ArgumentParser(description="Unified script for Blind and Targeted molecular docking using Vina, RDKit, and Meeko.")
    
    # Required Arguments
    parser.add_argument('-r', '--receptor', required=True, help="Path to the receptor PDBQT file")
    parser.add_argument('-o', '--out_dir', required=True, help="Directory to save the resulting PDBQT poses")
    
    # Adding optional argument for docking with SMILES string of the ligand
    parser.add_argument('-l', '--ligand', required=False, help="Path to the ligand SDF file")
    parser.add_argument('-ls', '--ligand_smiles', required=False, help="SMILES string of the ligand (optional)")

    # Optional Arguments for Targeted Docking
    parser.add_argument('--center', type=float, nargs=3, help="Center coordinates (X Y Z). If omitted, blind docking is triggered.")
    parser.add_argument('--box_size', type=float, nargs=3, default=[15.0, 15.0, 15.0], help="Box dimensions (X Y Z). Default: 15 15 15. Used only if --center is provided.")
    
    # Optional Vina Parameters
    parser.add_argument('--exhaustiveness', type=int, default=64, help="Exhaustiveness of the search (default: 64)")
    parser.add_argument('--n_poses', type=int, default=20, help="Number of poses to output (default: 20)")

    args = parser.parse_args()

    # 1. Setup paths and output directory
    os.makedirs(args.out_dir, exist_ok=True)
    rec_path = Path(args.receptor)
    
    if args.ligand_smiles:
        lig_name = Chem.MolToInchiKey(Chem.MolFromSmiles(args.ligand_smiles))[:8]

    else: 
        lig_path = Path(args.ligand)
        lig_name = lig_path.stem
    
    # Extract base names for dynamic file naming (e.g., '6T0Y' and 'BADGE')
    rec_name = rec_path.stem.replace('_clean_h_meeko', '') 

    # 2. Determine Mode & Grid Box Parameters
    if args.center:
        mode = "targeted"
        center_coords = args.center
        box_dims = args.box_size
        # Format filename with coordinates to 1 decimal place to avoid messy file names
        coords_str = f"{center_coords[0]:.1f}_{center_coords[1]:.1f}_{center_coords[2]:.1f}"
        out_filename = f"{rec_name}_{lig_name}_{mode}_{coords_str}.pdbqt"
    else:
        mode = "blind"
        print(f"Calculating grid box for {mode} docking...")
        center_coords, box_dims = get_receptor_center_and_size(args.receptor)
        out_filename = f"{rec_name}_{lig_name}_{mode}.pdbqt"

    out_filepath = os.path.join(args.out_dir, out_filename)

    print("-" * 40)
    print(f"Mode:           {mode.capitalize()} Docking")
    print(f"Receptor:       {args.receptor}")
    print(f"Ligand:         {args.ligand}")
    print(f"Center Coords:  {[round(c, 3) for c in center_coords]}")
    print(f"Box Dimensions: {[round(b, 3) for b in box_dims]}")
    print(f"Output File:    {out_filepath}")
    print("-" * 40)

    # 3. Initialize Vina & Receptor
    v = Vina(sf_name='vina')
    v.set_receptor(str(rec_path))

    # 4. Prepare Ligand (RDKit -> Meeko)
    print("Preparing ligand with RDKit and Meeko...")

    if args.ligand_smiles:
        mol = Chem.MolFromSmiles(args.ligand_smiles)
        if mol is None:
            raise ValueError(f"Failed to parse SMILES: {args.ligand_smiles}")

    else:
        suppl = Chem.SDMolSupplier(str(lig_path))
        mol = next(suppl)
        if mol is None:
            raise ValueError(f"Failed to read a valid molecule from {lig_path}")
        
    mol_with_h = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol_with_h, AllChem.ETKDG())
    AllChem.MMFFOptimizeMolecule(mol_with_h)
    
    preparer = MoleculePreparation()
    preparer.prepare(mol_with_h)
    pdbqt_string = preparer.write_pdbqt_string()

    # 5. Compute Maps & Set Ligand
    print("Computing Vina maps...")
    v.compute_vina_maps(center=center_coords, box_size=box_dims)
    v.set_ligand_from_string(pdbqt_string)

    # 6. Run Docking
    print(f"Docking (Exhaustiveness: {args.exhaustiveness}, Poses: {args.n_poses})...")
    v.dock(exhaustiveness=args.exhaustiveness, n_poses=args.n_poses)

    # 7. Write Results
    v.write_poses(out_filepath, n_poses=args.n_poses, overwrite=True)
    print(f"Docking complete. Results saved to {out_filepath}")

if __name__ == "__main__":
    main()
