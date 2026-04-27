
from vina import Vina
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation

v = Vina(sf_name='vina')

# Usage in your script:
receptor_path = '../molecules/6T0Y_clean_h_meeko.pdbqt'
v.set_receptor(receptor_path)

# if you want to focus on the active site open Pymol click on a residue near the area of interest 
# and then in Pymol's command line type "center (sele)" and then "get_position". You will get the coordinates
# of the selected residue which can be used in the "center_coords" array below

# best pose
#center_coords = [20, 9, 15]

# best pose
# center_coords = [  23.796,   2.233,  10.547]
# second best pose
# center_coords = [  15.610,  -3.834,  22.129]
# third best pose
center_coords = [  18.802,  32.168,  14.933]

box_dims = [15, 15, 15]

suppl = Chem.SDMolSupplier('../molecules/BADGE.sdf')
mol = next(suppl)
mol_with_h = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol_with_h, AllChem.ETKDG())
AllChem.MMFFOptimizeMolecule(mol_with_h)
preparer = MoleculePreparation()
preparer.prepare(mol_with_h)
pdbqt_string = preparer.write_pdbqt_string()
# v.set_ligand_from_file('../molecules/BADGE_meeko.pdbqt')

# Use the calculated values
v.compute_vina_maps(center=center_coords, box_size=box_dims)

v.set_ligand_from_string(pdbqt_string)

v.dock(exhaustiveness=64, n_poses=20)
v.write_poses('../molecules/6T0Y_BADGE_targeted_docking_results_3.pdbqt', n_poses=20, overwrite=True)