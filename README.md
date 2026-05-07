# DTU BioBuilders 2026 — Epoxy Resin Degradation Pipeline

Computational pipeline for engineering unspecific peroxygenases (UPOs) to degrade epoxy resin monomers, developed for iGEM 2026.

---

## Repository structure

```
Pipeline/
│
├── 01_Easy_track/
│   ├── claude_easytrack_pipeline/             — Claude-assisted sequence screening
│   └── claude_easytrack_pipeline_with_blast/  — same + BLAST homology search
│
├── 02_Medium_track/
│   ├── automated_rational_mutagenesis_pipeline/  — Boltz-2 + Claude mutation design
│   ├── enzyme_discovery/                          — novel UPO discovery from sequence databases
│   ├── fetching_AF_from_uniprot/                  — AlphaFold structure retrieval by UniProt ID
│   └── Kristian_course_methods/                   — structure-based methods from course
│
├── 03_Advanced_track/
│   └── DeNovoDesign/                          — de novo enzyme design
│
├── 02_Medium_track/molecular_dynamics/
│   └── MCPB_gamess/                           — Full MD pipeline for PaDa-I UPO (5OXU)
│       └── README.md                          ← start here for MD
│
├── PerEnzyme/
│   ├── UPOs/                                  — per-enzyme configs and results for UPOs
│   └── laccases/                              — per-enzyme configs and results for laccases
│
├── OtherStuff/
│   ├── BioBrick/                              — BioBrick submission materials
│   └── streamlit_apps/                        — interactive visualisation apps
│
├── pipeline.py                                — top-level AI-rational design loop (Boltz-2 + Claude)
├── dock_nnbt.py                               — AutoDock Vina docking for NNBT substrate
└── setup.sh                                   — dependency installation guide
```

---

## Tracks

**Easy track** — sequence-based screening and BLAST homology search to find candidate enzymes.

**Medium track** — structure-based rational mutagenesis using Boltz-2 co-folding and Claude for reasoning. Iteratively proposes mutations, evaluates binding affinity predictions, and selects variants.

**Advanced track / MD** — quantum-mechanically parameterized molecular dynamics of PaDa-I UPO. Full MCPB.py + GAMESS + OpenMM pipeline with explicit heme-Fe force field. See [`MD/README.md`](MD/README.md).

---

## Getting started

Each subdirectory has its own `README.md` with setup instructions. If you add new work, please add a short README to your folder explaining what's in there and how to run it.

For the MD pipeline specifically: [`02_Medium_track/molecular_dynamics/MCPB_gamess/README.md`](02_Medium_track/molecular_dynamics/MCPB_gamess/README.md).
