
# Investigating UPOs as potential enzymes for breakdown of epoxy resins

## Task list
- [ ] Here's a task that we need to do
- [x] Here's a task that we have completed

## About UPOs - chat quick overview

### Key properties
Use H₂O₂ as the only cosubstrate
Catalyze:
Hydroxylation of C–H bonds
Epoxidation of alkenes
Ether and aromatic oxidation
Work on bulky, hydrophobic, xenobiotic compounds
Often extracellular → easier recovery and handling


### Why are UPOs interesting for epoxy resin breakdown?
Epoxy resins (even simplified single-chain models) are hard to degrade because they typically contain:
Aromatic rings (e.g. bisphenol-like structures)
Ether bonds (C–O–C)
Highly stable C–C backbones
UPOs are relevant because they can attack exactly these features.
Potential reactions UPOs could catalyze on epoxy-like compounds
1. Aromatic hydroxylation
UPOs can hydroxylate phenyl rings:
Adds –OH groups
Makes the molecule more polar
Creates “handles” for:
Further enzymatic attack
Spontaneous chemical instability
Downstream metabolism
This is huge for epoxy monomers derived from bisphenols.
2. Ether bond destabilization (indirect cleavage)
UPOs typically don’t “snap” ether bonds directly, but:
Oxidation near the ether linkage
Formation of phenols or alcohols
Increased susceptibility to:
Hydrolysis
Secondary enzymes (laccases, esterases, hydrolases)
For an iGEM project, this fits well into a multi-enzyme degradation pathway concept.
3. Epoxide ring modification (model systems)
If your single-chain epoxy model still contains epoxide-like structures:
UPOs can oxidize adjacent carbons
This can:
Destabilize the epoxide
Promote ring opening
Reduce cross-linking potential
That’s perfect for proof-of-concept degradation studies.
Why UPOs are attractive for iGEM specifically
From a synthetic biology / iGEM perspective, UPOs check a lot of boxes:
Practical advantages
Only require H₂O₂ (easy to supply enzymatically)
Can be expressed in yeast or fungi (and sometimes engineered bacteria)
Act on non-natural substrates (plastics, resins, pollutants)
Conceptual strength
Fit well into:
“Plastic-degrading enzyme toolbox”
Sustainable materials research
Enzyme cascade designs
Strong link to bioremediation & circular bioeconomy

## Methods - what did we do

### Uniprot search
Initially, an advanced search was carried out on Uniprot, to find potential Unspecific Peroxygenases (UPOs) and investigate their functions. UPOs are fungal enzymes that include a heme group in their active sites. 

Advanced search: “(protein_name:"unspecific peroxygenase")”
Got 64 hits. All unreviewed from the TrEMBL.

Downloaded fasta file with all results:
"uniprotkb_unspecific_peroxygenase_2026_02_07.fasta"

https://pubs.acs.org/doi/10.1021/acschembio.4c00504

#### Looking at structures for one of them:
PDB from first hit: 8RNJ
Downloaded from pdb: "Download Files" -> "PDBx/mmCIF Format". Called 8RNJ.cif
Opened in PyMol and cleaned with these commands:

PyMOL>remove solvent
 Remove: eliminated 444 atoms in model "8RNJ".

PyMOL>remove chain B
 Remove: eliminated 1951 atoms in model "8RNJ".

Selected all organic residues except "HEM" and then:
PyMOL>set_name sele, organic_except_heme

PyMOL>cmd.hide("everything","organic_except_heme")

Saved the cleaned structure via:
 Save: wrote "~/BioBuilders_2026/UPOs/8RNJ_cleaned.cif".

![alt text](8RNJ_cleaned.png)


### PDB search
searched advanced search: structure title has exact phrase Unspecific Peroxygenase
resulting in 50 structures
downloaded to pdb files in batch located in this folder: PDBsearch_structures_07_02_2026
pdb files were unzipped via gunzip UPOs/PDBsearch_structures_07_02_2026/* in terminal

Maybe do the same with these - or qualitative investigation of which are better


### Doing Multiple Sequence Alignment

Using the instructions from Kristians course:


## Other possible methods / things to do.

**Maybe try Molecular Docking of cleaned 8RNJ with the substrate**
to see if it can tell us something about the binding of the ligand.

From my notes:
“used to predict the three-dimensional (3D) structure of the complex between a protein receptor and small molecule (i.e., a drug) and the affinity between them”

webserver: https://www.swissdock.ch/ 

**Preparing the pdb file:**

Folding sequence with AF2 to include the missing regions
Truncate the PDB file by removing the long loop tail in the end 
have to remove all heteroatoms - DONT THINK ITS THE CASE FOR US SINCE HEME GROUP IS PART OF ACTIVE SITE???
‘remove hetatm'

Download compound from PubChem and convert to .mol2 file via export Molecule function (DONE)

**Steps for docking:**

Autodock Vina på swissdock.ch

“Submit ligand”

“Submit target”

define the search box based on knowledge about active site.

“Select parameters” = maximum 10. - is something about how much it searches - it's not exhaustive search.
“Increase the sampling exhaustivity to increase the amount of computational effort.”

Download results as zip files.

**Steps for analyzing results:**

- open with PyMOL all the binding poses file **vina_dock.pdbqt** and the receptor file **system.pdbqt** along with reference structure. ****
- Examine each of the docking binding poses and explore the surrounding environment of the ligand and identify the protein residues involved in the binding.
- Compare the binding poses with respect to the reference structure.
- Also consider the score associated with each binding pose (from the webserver or within the vina_dock.pdbqt).