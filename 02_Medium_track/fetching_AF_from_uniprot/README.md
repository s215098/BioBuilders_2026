# fetch_alphafold.py

Downloads AlphaFold PDB structures for a given UniProt ID. Structures are saved to a `fetched_structures/` subfolder.

## How it works

1. Queries the **UniProt REST API** for all AlphaFold DB cross-references for the entry.
2. For each AlphaFold entry ID, calls the **AlphaFold EBI API** to get the PDB download URL, structure name, mean pLDDT confidence score, and number of chains.
3. If multiple predictions exist for the same entry, picks the one with the **highest mean pLDDT**.
4. Downloads each PDB file into `fetched_structures/`, named as:
   ```
   {Description}_AF-{EntryID}_{ModelVersion}.pdb
   ```
   e.g. `Hemoglobin_subunit_beta_AF-P68871-F1_v4.pdb`

## Requirements

- Python 3.6+ (standard library only — no `pip install` needed)

## Usage

```bash
python fetch_alphafold.py <UniProtID> [--monomer-only]
```

| Argument         | Description                                                      |
|------------------|------------------------------------------------------------------|
| `UniProtID`      | UniProt accession number, e.g. `P68871`                          |
| `--monomer-only` | Skip multi-chain entries (homodimers, complexes). Recommended for most pipelines. |

## Examples

```bash
# Hemoglobin subunit beta — single fragment
python fetch_alphafold.py P68871

# Titin — very long protein, downloaded as multiple fragments (F1, F2, ...)
python fetch_alphafold.py Q8WZ42

# Monomers only — skip any homodimer or complex entries
python fetch_alphafold.py P68871 --monomer-only
```

## Example output

```
Looking up AlphaFold structures for UniProt ID: P68871
Found 1 AlphaFold entry/entries: P68871

Fetching metadata for P68871 ...
  Name       : Hemoglobin subunit beta
  Confidence : 96.4 mean pLDDT
  Chains     : 1
  Model ver  : v4
  Entry URL  : https://alphafold.ebi.ac.uk/entry/AF-P68871-F1
  Saving     : fetched_structures/Hemoglobin_subunit_beta_AF-P68871-F1_v4.pdb
  Done.

Manifest saved: fetched_structures/P68871_manifest.tsv
```

## Output files

```
fetched_structures/
└── Hemoglobin_subunit_beta_AF-P68871-F1_v4.pdb
```

## Finding a structure on the AlphaFold website

The version number in the filename (`v4`, `v6`, etc.) is **not** part of the searchable identifier on the website. Use the `entry_url` column in the manifest to go directly to the entry page, e.g.:

```
https://alphafold.ebi.ac.uk/entry/AF-P68871-F1
```

The website always shows the latest model version with the interactive 3D viewer and pLDDT/PAE plots.


## Notes
- **Fragments**: large proteins are split into overlapping fragments (`F1`, `F2`, ...). All fragments are downloaded; they are complementary, not alternative predictions.
- **API transition**: the AlphaFold EBI API is renaming some fields (e.g. `entryId` → `modelEntityId`) with a sunset date of June 2026. This script handles both field names automatically.