#!/usr/bin/env python

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
import argparse
import sys

# --- Argument parsing ---
parser = argparse.ArgumentParser(description="Plot normalized solvation")
parser.add_argument('--input', required=True, help="Path to normalized result file (.csv)")
parser.add_argument('--stride', required=True, help="Path to STRIDE output file")
parser.add_argument('--output', default='normalized.png', help='Output file name (default: normalized.png)')
args = parser.parse_args()

# --- Load normalized ASA file ---
df = pd.read_csv(args.input, sep='\t')

# Verify expected columns are present
required_cols = {'ResNum', 'Normalized_ASA'}
if not required_cols.issubset(df.columns):
    print(f"Error: input file must contain columns: {required_cols}")
    print(f"Found columns: {set(df.columns)}")
    sys.exit(1)

# --- Parse secondary structure from STRIDE file ---
ss_map = {}
with open(args.stride) as f:
    for line in f:
        if line.startswith("ASG"):
            parts = line.split()
            res_num = int(parts[3])
            ss      = parts[5]   # H, E, T, C, G, B etc.
            ss_map[res_num] = ss

if not ss_map:
    print("Warning: no ASG lines found in STRIDE file — check that the file is correct")
    sys.exit(1)

# Map secondary structure onto dataframe rows
df['SS'] = df['ResNum'].map(ss_map)

# --- Color each secondary structure type ---
colors = {
    'H': '#e05c5c',   # alpha helix - red
    'G': '#e07b5c',   # 310 helix - orange
    'E': '#5c7be0',   # strand - blue
    'B': '#8fb3f5',   # bridge - light blue
    'T': '#f0c040',   # turn - yellow
    'C': '#aaaaaa',   # coil - gray
}

# .tolist() converts the pandas Series to a plain list,
# which is required for older matplotlib versions
bar_colors = df['SS'].map(colors).fillna('#aaaaaa').tolist()

# --- Plot ---
fig, ax = plt.subplots(figsize=(18, 4))
ax.bar(df['ResNum'], df['Normalized_ASA'], color=bar_colors, width=1.0)

ax.set_xlabel('Residue number')
ax.set_ylabel('Normalized ASA')
ax.set_title('Solvent accessibility per residue')
ax.set_xlim(df['ResNum'].min() - 1, df['ResNum'].max() + 1)
ax.set_ylim(0, 1.2)

# Dashed reference line at 1.0 (fully exposed)
ax.axhline(1.0, color='black', linewidth=0.5, linestyle='--', alpha=0.4)

# --- Legend ---
legend_labels = {
    'H': 'Alpha helix',
    'G': '310 helix',
    'E': 'Strand',
    'B': 'Bridge',
    'T': 'Turn',
    'C': 'Coil'
}

# Only include secondary structure types that actually appear in the data
present_ss = df['SS'].dropna().unique()
patches = [
    mpatches.Patch(color=colors[k], label=v)
    for k, v in legend_labels.items()
    if k in present_ss
]
ax.legend(handles=patches, loc='upper right', fontsize=8)

plt.tight_layout()
plt.savefig(args.output, dpi=150)
print(f"Plot saved to {args.output}")
plt.show()