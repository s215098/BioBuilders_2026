# Molecular Docking

To see how candidate enzymes are predicted to bind with the ligands.

# Principals
Delta G etc...

# Methods

## Getting the ligand - (Example for BADGE)

**Ligand download and preparation**

Download the correct compound from PubChem, one of the most popular databases of small molecules (https://pubchem.ncbi.nlm.nih.gov/ ). Search for the compound's name in the search bar and download the 3D conformer as .sdf file, and name the downloaded file as [ligand].sdf. Visualize the sdf file in PyMol.

![BADGE.sdf file visualized in PyMol](images/BADGE_sdf.png)

**Creation of the ligand file**

The 3D structure of the ligand must be converted in .mol2 file. To do so load, the .sdf into PyMOL and use the export Molecule function in PyMOL and save it as ligand.mol2.

## Prepare enzyme

**Protein Receptor download and preparation**
Download the PDB file of the BADGE protein. First, visualize the 3D structure in PyMOL and analyze in which region of the protein the ligand compound is interacting (commonly referred as “pocket”) and study its binding mode.


The 3D structure needs firstly to be prepared for the docking protocol. To this point, you should firstly remove from the PDB file all non-essential protein atoms, solvent or duplicated chains. I.e. doing:
- remove chain B #remove chain B
- remove resn S #remove residue S.
- remove solvent #remove solvent
- remove hetatm #remove heteroatoms

When done rename the file to [enzyme]_cleaned.pdb

## Perform Docking
**Run docking calculation on SwissDock**

SwissDock supports two state-of-the-art docking programs, i.e. Attracting Cavities and AutoDock Vina, not blind docking, and requires knowledge of possible binding regions of the ligand on the 3D structure of the target protein. This information should be manually provided to SwissDock on the web server interface under “**Define search space”**. In our case, we know the pocket of interest in the protein from the original crystal structure and literature, being the active site. You will perform the docking calculation using *Autodock Vina* with the prepared input files at https://www.swissdock.ch/ . In the first step under “**Submit ligand**”, you should upload the prepared ligand structure (ligand.mol2) and click on “**Prepare Ligand”**. Next, under “**Submit a target”** you should upload the prepared receptor structure (3udh_clean.pdb). Then, you should define the search box for the docking procedure based on the knowledge of the active site from the original crystal structure. Finally, you should select the exhaustivity, under “**Select parameters”**. As a project name for the docking calculation, you must follow [enzyme]_[ligand]_[search box center]_[search box size]_[exhaustivity], making it easy to retrieve the correct results. When finished, you can download your results as ZIP files.


**Analysis of the docking results**

To examine the results, first download the output files from the web server. Once you have downloaded the ZIP, open with PyMOL all the binding poses file **vina_dock.pdbqt** and the receptor file **system.pdbqt**, along with the reference structure **[enzyme].pdb**. Examine in detail each of the docking binding poses and explore the surrounding environment of the ligand and identify the protein residues involved in the binding. Compare the binding poses with respect to the X-ray structure. While investigating the different poses of the ligands, also consider the score associated with each binding pose (from the webserver or within the vina_dock.pdbqt).