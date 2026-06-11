#!/opt/anaconda3/bin/python

"""
nomalization.py ...
"""

import pandas as pd
import argparse
import numpy as np


# --- Argument parsing ---
parser = argparse.ArgumentParser(description="Normalize ASA values from STRIDE output")
parser.add_argument('--stride', required=True, help="Path to normalization .txt file")
parser.add_argument('--column',
                    choices=['Theoretical', 'Empirical', 'Miller et al', 'Rose et al'],
                    default='Empirical',
                    help="Normalization values to use (default: Empirical)")
parser.add_argument('--out', default='normalized_ASA.txt', help='Output file name (default: normalized_ASA.txt)')
parser.add_argument('--pdb', help='Path to input PDB file')
args = parser.parse_args()



# --- Column mapping ---
col_map = {
    'Theoretical':  'Theoretical',
    'Empirical':    'Empirical',
    'Miller et al': 'Miller et al. (1987)',
    'Rose et al':   'Rose et al. (1985)'
}

norm_col = col_map[args.column]



# --- Parse STRIDE AREA lines ---
records = []

# open the file with the stride results: 
with open(args.stride) as f:
    for line in f:

        # The relevant lines (the results we want) are starting with "ASG":
        if line.startswith("ASG"):

            parts = line.split()
            residue   = parts[1]          # e.g. GLY
            chain     = parts[2]          # e.g. A
            res_num   = int(parts[3])     # e.g. 3
            area      = float(parts[-2])   # solvent accessible area
            records.append({
                'Residue':  residue,
                'Chain':    chain,
                'ResNum':   res_num,
                'Area':     area
            })

stride_df = pd.DataFrame(records)


# # --- Load normalization table ---
# # Define the file with the normalization values, from paper. 
norm_df = pd.read_csv("normalization.txt", sep="\t")
print("Norm_df: ", norm_df)



# --- Merge on residue name ---

three_to_full = {
    'ALA': 'Alanine',
    'ARG': 'Arginine',
    'ASN': 'Asparagine',
    'ASP': 'Aspartate',
    'CYS': 'Cysteine',
    'GLU': 'Glutamate',
    'GLN': 'Glutamine',
    'GLY': 'Glycine',
    'HIS': 'Histidine',
    'ILE': 'Isoleucine',
    'LEU': 'Leucine',
    'LYS': 'Lysine',
    'MET': 'Methionine',
    'PHE': 'Phenylalanine',
    'PRO': 'Proline',
    'SER': 'Serine',
    'THR': 'Threonine',
    'TRP': 'Tryptophan',
    'TYR': 'Tyrosine',
    'VAL': 'Valine'
}


# Convert 3-letter codes to full names to match normalization table
stride_df['Residue'] = stride_df['Residue'].map(three_to_full)

# .merge() joins two dataframes like a SQL JOIN
# on='Residue' means match rows where the Residue name is the same in both dataframes (e.g. 'GLY' matches 'GLY')
# how='left' means keep ALL rows from stride_df (left), even if no match is found in norm_df
# If a residue has no match in norm_df, the normalization column will just be NaN
merged = stride_df.merge(norm_df[['Residue', norm_col]], on='Residue', how='left')




# # --- Normalize ---

# # Divides every value in the 'Area' column by the corresponding normalization value in the same row
# # pandas does this element-wise automatically (row 1 / row 1, row 2 / row 2, etc.)
# # The result is stored as a new column called 'Normalized_ASA'
merged['Normalized_ASA'] = merged['Area'] / merged[norm_col]




# --- Saving results to csv ---
# Save to a tab-separated file
outfile_name = str(args.out) + ".csv"
merged[['ResNum', 'Residue', 'Chain', 'Area', norm_col, 'Normalized_ASA']].to_csv(outfile_name, sep='\t', index=False)




# --- Coloring PyMol if called
if args.pdb:
        # --- Write PyMOL coloring script ---
    pml_out =  outfile_name.replace('.csv', '.pml')

    with open(pml_out, 'w') as pml:
        # Load the structure
        pml.write(f"load {args.pdb}\n")
        pml.write("bg_color white\n")
        pml.write("hide everything\n")
        pml.write("show cartoon\n\n")
        
        # Set all B-factors to 0 first
        pml.write("alter all, b=0\n")
        
        # Write a b-factor value for each residue
        for _, row in merged.iterrows():
            val = round(row['Normalized_ASA'], 3)
            pml.write(f"alter (chain {row['Chain']} and resi {row['ResNum']}), b={val}\n")
        
        # Color by B-factor (0=buried, 1=exposed)
        pml.write("spectrum b, blue_white_red, minimum=0, maximum=1\n") # Cool to warm (blue=buried, red=exposed) - default

        pml.write("rebuild\n")

    print(f"PyMOL script written to {pml_out}")