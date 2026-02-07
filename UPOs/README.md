
# Investigating UPOs as potential enzymes for breakdown of epoxy resins

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


## Uniprot search
Initially, an advanced search was carried out on Uniprot, to find potential Unspecific Peroxygenases (UPOs) and investigate their functions. UPOs are fungal enzymes that include a heme group in their active sites. 

Advanced search: “(protein_name:"unspecific peroxygenase")”
Got 64 hits. All unreviewed from the TrEMBL.

Downloaded fasta file with all results:
"uniprotkb_unspecific_peroxygenase_2026_02_07.fasta"

### Looking at structures for one of them:
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