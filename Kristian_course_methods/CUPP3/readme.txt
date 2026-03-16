README FOR RUNNING CUPP CLUSTERING AND VISUALIZATION
This document explains how to run the CUPP workflow using run.py.
It is written entirely as plain text so it can be placed directly in a README file.


INPUT FASTA FILE
You must provide a FASTA file with the protein sequences you want to cluster.
Two options:
A: Edit run.py. Change curr_fam to the name you want, for example curr_fam = "PF14346". The script will then expect a file named PF14346.faa.
B: Put your FASTA file in the same folder as the scripts and call it OWN1.faa. Then run.py will automatically pick it up.


HOMOLOGY REDUCTION TO 90 PERCENT (MANDATORY)
Before running CUPP, you must ensure the input sequences are homology-reduced to 90 percent identity.
This is done using MMseqs2. The recommended tool is at:
https://toolkit.tuebingen.mpg.de/tools/mmseqs2
When you start run.py, it asks:
"Have you homology reduced the sequences to 90 percent using mmseqs2?"
If you answer n, the script stops.
This requirement ensures consistent clustering quality.


CLEANING FASTA HEADERS
After confirming homology reduction, run.py cleans your FASTA headers.
It removes prefixes like sp| or tr|, strips long descriptions, and normalizes the header format so that accessions match UniProt and iTOL requirements.
The cleaned FASTA is saved as NAME_cleaned.faa.
All later steps use the cleaned file.


CUPP CLUSTERING
run.py then runs CUPPclustering_DIRECT.py with:
clustering enabled
domain mode off
CD-HIT disabled
minimum CUPP group size = 5
This produces the CUPP pool JSON:
CUPP/CUPPpools/<name>_CUPPpool.json</name>


BUILDING THE CUPP LIBRARY
The script then runs CUPPclustering_DIRECT.py -recom to build a CUPP library from the pool.
The result appears at:
CUPP/CUPPlibrary/8x2_90_CUPP_lib_<name>_CUPPlibrary.json</name>


CUPP PREDICTION
Next, CUPPprediction_DIRECT.py is run.
It uses the cleaned FASTA and the library to assign CUPP groups to each input protein.


GENERATING THE VISUALIZATION AND ITOL DATASETS
Finally the script runs CUPPvisualization.py.
This uses:
the CUPP pool
the cleaned FASTA file
the CUPP tree (<name>_fa8x2_90.tree)
It produces a folder:
CUPP/itol/<name>/additional/</name></name>


This folder contains all iTOL datasets generated from CUPP groups and from UniProt metadata.


CUPP GROUP VISUALIZATION FILES
The following datasets are created for CUPP groups:
color strip
text labels
single representative text label
symbol layers
arrow (connection) dataset
The arrow dataset uses iTOL’s DATASET_CONNECTION format and draws arrows from a single anchor sequence to all other sequences in the same CUPP group.


UNIPROT METADATA VISUALIZATION
The script automatically downloads UniProt metadata and generates the following datasets:
Protein descriptions
GO terms
Pfam domains
Taxonomic ranks including species, genus, family, order, class, phylum, kingdom, superkingdom
Each of these is provided as:
color strip
text labels
single text label
symbol dataset
arrow dataset connecting one anchor to all group members


MOUSEOVER POPUP BOXES (POPUP_INFO)
The script creates:
popup_info.txt
This file uses iTOL’s POPUP_INFO template.
Each entry contains a small HTML card with:
Accession
Protein name
Organism
Sequence length
Protein existence code
EC number if present
Pfam domains (shortened if long)
GO terms (shortened if long)
Direct link to UniProt
Direct link to AlphaFold if the protein is available in AlphaFold
If metadata text is too long, it is automatically truncated.


HOW TO RUN THE WORKFLOW
Open a terminal and type:
python3 run.py
The script will:
check homology reduction
clean the FASTA
run all CUPP clustering and prediction steps
generate iTOL datasets


The final visualization files are located under:
CUPP/itol/<name>/
Import the tree into iTOL
Drag all dataset files onto the tree in iTOL (https://itol.embl.de/upload.cgi)

COMMON PROBLEMS
If iTOL rejects a dataset:
check that the SEPARATOR matches the dataset
TEXT datasets must use COMMA
POPUP_INFO uses TAB
Connections (arrows) use COMMA

Make sure the IDs in the datasets exactly match the leaf names of your tree.
The visualization script already filters IDs to only those present in the tree.
Long metadata can cause issues. The script automatically shortens long entries.