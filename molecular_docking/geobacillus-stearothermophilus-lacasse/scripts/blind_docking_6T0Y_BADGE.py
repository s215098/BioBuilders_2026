
from vina import Vina
import numpy as np

def get_receptor_center_and_size(pdb_file):
    coordinates = []
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                # Extract x, y, z from standard PDB/PDBQT columns
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coordinates.append([x, y, z])
    
    coords = np.array(coordinates)
    center = np.mean(coords, axis=0) # Geometric center
    
    # Calculate dimensions for blind docking
    # (Max - Min) per axis gives the full span of the protein
    size = np.max(coords, axis=0) - np.min(coords, axis=0)
    
    # Add a buffer (e.g., 5-10 Angstroms) so the ligand can 
    # fully bind to surface residues
    size = size + 10.0
    
    return center.tolist(), size.tolist()

# Usage in your script:
receptor_path = '../molecules/6T0Y_clean_h_meeko.pdbqt'
center_coords, box_dims = get_receptor_center_and_size(receptor_path)

print(f"Calculated Center: {center_coords}")
print(f"Calculated Box Size: {box_dims}")


v = Vina(sf_name='vina')
v.set_receptor(receptor_path)
v.set_ligand_from_file('../molecules/BADGE_meeko.pdbqt')

# Use the calculated values
v.compute_vina_maps(center=center_coords, box_size=box_dims)

# Blind docking typically requires higher exhaustiveness 
# because the search space is much larger than a pocket
v.dock(exhaustiveness=64, n_poses=20)
v.write_poses('../molecules/6T0Y_BADGE_blind_docking_results.pdbqt', n_poses=10, overwrite=True)