
# Creating Virtual Environments

## Protocol to create the "requirements.txt" and "environment.yml" files

1. "python3 -m pip freeze > requirements.txt"

2. conda env export --no-builds > environment.yml

## Protocol to create and activate the enrivonments in your local machine

1. "conda create -n biobuilders python=3.12.11"

2. "conda activate biobuilders"

3. "python3 -m pip install -r requirements.txt"


# Protocol for Molecular Docking of Geobacillus Laccase (6T0Y)

## Prepare Receptor for Docking (working under the molecules/ directory)

1. Open 6T0Y with Pymol

2. In the commands prompt type: "remove solvent", to remove the water molecules

3. Then type "remove resn ZN" to remove the ZINC cofactors.

3. Finally type "save 6T0Y_clean.pdb", to save the custom molecule

4. If additionally you want to protonate you structure type "h_add".

5. And then save with "save 6T0Y_clean_h.pdb" and close Pymol

6. In the command line now type: "mk_prepare_receptor.py -i 6T0Y_clean_h.pdb -o 6T0Y_clean_h_meeko -p --default_altloc A -a", to convert the pdb format to pdbqt.

## Prepare Ligand for Docking (working under the molecules/ directory)

1. It is supposed you already have downloaded the 3D structure of BADGE from pubchem (BADGE.sdf).

2. Type in the command line: "mk_prepare_ligand.py -i BADGE.sdf -o BADGE_meeko.pdbqt", to convert the sdf to pdbqt format.

## Perform Molecular Docking (working under the scripts/ directory)

1. Run the docking script: "python3 blind_docking_6T0Y_BADGE.py"

2. Change directory to "molecules" where output of the script is written.

3. Open 6TOY with Pymol and then in the same Pymol session open the "6T0Y_BADGE_blind_docking_results.pdbqt" file

4. Iterate through different docking poses using the arrows in the bottom right panel