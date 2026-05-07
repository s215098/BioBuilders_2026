# UPO Enzyme Discovery Pipeline
## DTU BioBuilders 2026

Automated pipeline for retrieving and classifying UnSpecific Peroxygenase (UPO)
structures from CUPP clade representatives, preparing them for docking.

---

## What does it do?

Given a list of **UniProt accessions**, the pipeline:

1. **Retrieves 3D structures** — tries experimental PDB entries first (best
   resolution), then AlphaFold DB, then flags anything left over for local
   Boltz-2 prediction.
2. **Classifies each enzyme** as long-UPO or short-UPO — using sequence
   length, pairwise alignment against AaeUPO (long) and MroUPO (short), and
   optionally TM-score structural comparison (requires TMalign).
3. **Quality-checks the active site** — verifies the Cys axial heme ligand
   (the only universal UPO feature checked as pass/fail), and reports
   informational annotations (aromatic residues near Fe, channel occupancy,
   CA-RMSD) that are recorded but do not affect the pass/fail result.

After the pipeline, `graft_heme.py` can place heme into predicted structures
(Boltz-2 / AlphaFold models that lack it) by alignment to a reference.

---

## Input

Edit **`config.yaml`** before running. There are two things to configure:

### 1 — Query accessions

```yaml
queries:
  - A0ACD6BA37   # clade 1
  - A0A8H5CVF9   # clade 10
  # ... one UniProt accession per UPO clade representative
```

Replace or extend this list with your own accessions. The pipeline fetches
sequences and structures from the internet at runtime — no local files needed
upfront.

### 2 — Reference structures (usually no change needed)

```yaml
references:
  long_upo:
    pdb_id: "5OXU"   # PaDa-I (evolved AaeUPO), 1.47 Å
    chain: "A"
  short_upo:
    pdb_id: "5FUJ"   # MroUPO, chain A
    chain: "A"
```

These are downloaded automatically and used for classification and active-site
QC. Only change them if you want to use different reference structures.

---

## Installation

```bash
pip install biopython requests numpy pyyaml tqdm

# Optional but strongly recommended — enables structural classification (Step 2 Tier 3):
# TMalign — https://zhanggroup.org/TM-align/
# Place the tmalign binary in your PATH or /usr/local/bin/
```

---

## Usage

Run from the `files/` directory (where `config.yaml` lives):

```bash
# Full pipeline (all 3 steps)
python run_pipeline.py --config config.yaml

# Run individual steps (useful when re-running after Boltz-2 predictions)
python run_pipeline.py --config config.yaml --steps 1
python run_pipeline.py --config config.yaml --steps 2,3

# Verbose output (shows debug-level details)
python run_pipeline.py --config config.yaml --verbose
```

---

## Outputs

All outputs are written to `results/`:

| File | Description |
|---|---|
| `structure_sources.tsv` | Source, PDB ID, resolution, file path per accession |
| `boltz_todo.txt` | Accessions where no structure was found — need Boltz-2 |
| `classification.tsv` | Length, alignment scores, TM-scores, final call per accession |
| `classification_summary.txt` | Human-readable long/short/ambiguous summary |
| `active_site_qc.tsv` | Cys axial (pass/fail), aromatics near Fe, RMSD, channel occupancy per structure |
| `logs/pipeline.log` | Full run log |
| `structures/pdb/` | Downloaded experimental PDB files |
| `structures/alphafold/` | Downloaded AlphaFold DB models |
| `structures/references/` | Reference structures (5OXU, 5FUJ) |
| `structures/boltz/` | Boltz-2 predictions (placed here manually — see below) |
| `structures/grafted/` | Heme-grafted structures ready for docking |

---

## Boltz-2 predictions (when needed)

After Step 1, check `results/boltz_todo.txt`. If it is non-empty, those
accessions need local Boltz-2 structure prediction before steps 2–3 can use
them.

```bash
# Install Boltz-2
pip install boltz

# FASTA files for all pending accessions are written to:
#   results/structures/boltz_input/<accession>.fasta

# Run prediction (batch)
boltz predict results/structures/boltz_input/ \
      --out_dir results/structures/boltz/ \
      --use_msa_server

# Optional: improve channel geometry by supplying pocket constraints from 5OXU
boltz predict results/structures/boltz_input/ \
      --pocket_conditioning results/structures/references/5oxu.pdb \
      --out_dir results/structures/boltz/

# After prediction, re-run classification and QC
python run_pipeline.py --config config.yaml --steps 2,3
```

Boltz-2 output PDB files must be placed (or symlinked) at:
`results/structures/boltz/<ACCESSION>.pdb`

---

## Heme grafting (after structure prediction)

Predicted structures (AlphaFold, Boltz-2) do not contain heme. Run
`graft_heme.py` to place the heme cofactor by structural alignment to a
reference before docking:

```bash
python graft_heme.py
```

Grafted structures are saved to `results/structures/grafted/`.

---

## Classification confidence

The classification uses up to three tiers and combines them by majority vote:

| Tier | Method | Requires |
|---|---|---|
| 1 | Sequence length (long ~330–420 aa, short ~220–310 aa) | UniProt accession |
| 2 | Pairwise alignment vs AaeUPO (long) and MroUPO (short) | UniProt accession |
| 3 | TM-score vs 5OXU / 5FUJ | TMalign + 3D structure |

- **High confidence** — all available tiers agree
- **Medium confidence** — sequence length + alignment agree, or one clear tier
- **Low confidence / ambiguous** — tiers disagree or only one tier available

Without TMalign, most calls will be low confidence. Installing TMalign is
recommended before treating any classification as final.

---

## Active site QC design

The QC step distinguishes between checks that are universal and checks that are
specific to the reference enzyme:

| Check | Type | Rationale |
|---|---|---|
| Cys axial ligand (SG within 3 Å of Fe) | **Pass/fail** | The Cys-thiolate Fe bond is the biochemical hallmark of all UPOs/P450s — absence means the structure is broken or the enzyme is misclassified |
| Heme present (experimental structures only) | **Pass/fail** | A crystallographic structure without heme is missing a critical cofactor — predicted structures (AlphaFold, Boltz-2) are expected to lack heme and pass QC regardless (`INFO_NO_HEME_PREDICTED_OK`) |
| Aromatics near Fe (Phe/Trp/Tyr within 8 Å) | Informational | AaeUPO has a characteristic Phe triad (Phe69/121/199) lining the substrate channel, but whether other clades conserve the same residues at equivalent positions is unknown — the count is recorded for manual inspection |
| Channel residue count (within 8 Å of Fe) | Informational | A rough fold-sanity proxy; threshold not calibrated for divergent clades |
| CA-RMSD to reference | Informational | Divergent clades will naturally have high RMSD to PaDa-I or MroUPO — this does not indicate a broken structure |

Informational flags appear in the `qc_flags` column prefixed with `INFO_` and
do not affect the `pass_qc` result.

---

## Reference structures

| Enzyme | PDB | Type | Used for |
|---|---|---|---|
| PaDa-I (evolved AaeUPO) | 5OXU (chain A) | Long-UPO | Classification, active site QC |
| MroUPO | 5FUJ (chain A) | Short-UPO | Classification, active site QC |

---

## Next steps after this pipeline

1. **Heme grafting** — run `graft_heme.py` for all predicted structures
2. **Docking** — AutoDock Vina or GNINA, box centred on grafted Fe atom
3. **Pose filtering** — productive pose = substrate C within ≤ 4.5 Å of Fe
4. **CUPP tree overlay** — map long/short calls and docking scores back to the
   phylogenetic tree
