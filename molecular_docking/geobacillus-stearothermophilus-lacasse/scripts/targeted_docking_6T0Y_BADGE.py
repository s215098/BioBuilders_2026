
from vina import Vina

# Usage in your script:
receptor_path = '../molecules/6T0Y_clean_h_meeko.pdbqt'

# if you want to focus on the active site open Pymol click on a residue near the area of interest 
# and then in Pymol's command line type "center (sele)" and then "get_position". You will get the coordinates
# of the selected residue which can be used in the "center_coords" array below
center_coords = [20, 9, 15]
box_dims = [20, 20, 20]

v = Vina(sf_name='vina')
v.set_receptor(receptor_path)
v.set_ligand_from_file('../molecules/BADGE_meeko.pdbqt')

# Use the calculated values
v.compute_vina_maps(center=center_coords, box_size=box_dims)

v.dock(exhaustiveness=64, n_poses=20)
v.write_poses('../molecules/6T0Y_BADGE_targeted_docking_results.pdbqt', n_poses=10, overwrite=True)