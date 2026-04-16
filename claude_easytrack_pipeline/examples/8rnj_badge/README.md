# Example Run: 8RNJ (Unspecific Peroxygenase) + BADGE

This folder documents the example run used to validate the pipeline.

## What is 8RNJ?

8RNJ is the crystal structure of an **Unspecific Peroxygenase (UPO)** — a fungal heme-thiolate enzyme that can oxidise a wide range of substrates using H₂O₂ as the oxidant. UPOs are of interest to our project because of their ability to degrade aromatic compounds, including potential bisphenol-derived pollutants.

- PDB: https://www.rcsb.org/structure/8RNJ
- Organism: *Marasmius wettsteinii*
- EC: 1.11.2.1

## What is BADGE?

Bisphenol A diglycidyl ether (BADGE) is an epoxy resin monomer and environmental contaminant. It is structurally related to BPA and has been shown to be cytotoxic and potentially endocrine-disrupting. Testing whether UPOs can bind (and therefore potentially oxidise) BADGE is directly relevant to our project goal.

- PubChem CID: 2286
- Molecular formula: C₂₁H₂₄O₄
- MW: 340.4 g/mol

## Running this example

```bash
# From the pipeline root directory:
conda activate easytrack
bash run_pipeline.sh config/example_8rnj_badge.yml
```

## Expected outputs

After the run completes, results are in `results/8rnj_badge/`:

- `sequences/raw_sequences.fasta` — ~200 UPO/peroxygenase sequences from NCBI
- `phylogeny/tree.png` — phylogenetic tree image
- `phylogeny/selections.txt` — one representative per cluster
- `receptor/8RNJ_clean_h.pdbqt` — receptor ready for Vina
- `ligand/BADGE.sdf` — BADGE 3D structure from PubChem
- `docking/8RNJ_BADGE_blind/poses.pdbqt` — 10 binding poses
- `docking/summary.csv` — binding affinities for all poses

## Interpreting results

Binding affinity is reported in **kcal/mol**. More negative = stronger predicted binding.

| Affinity | Interpretation |
|---|---|
| < -10 kcal/mol | Very strong (likely true binder) |
| -8 to -10 kcal/mol | Strong |
| -6 to -8 kcal/mol | Moderate |
| > -6 kcal/mol | Weak |

Note: these are **predicted** values. Experimental validation (e.g. kinetics, binding assays) is needed to confirm activity.

## Notes on the 8RNJ structure

The raw PDB file contains two chains (A and B). Our config keeps only chain A (`chains_to_keep: ["A"]`). The active site contains a heme group with a cysteine ligand (typical of UPOs). For targeted docking, the heme iron coordinates can be used as the box center.
