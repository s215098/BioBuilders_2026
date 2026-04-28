# fetch_alphafold.py

A minimal command-line script to download an AlphaFold protein structure using a UniProt ID.

## Requirements

- Python 3 (standard library only — no `pip install` needed)

## Usage

```bash
python fetch_alphafold.py <UniProtID> [pdb|cif]
```

The format argument is optional and defaults to `pdb`.

## Examples

```bash
# Download hemoglobin subunit beta as PDB (default)
python fetch_alphafold.py P68871

# Download as mmCIF format
python fetch_alphafold.py P68871 cif
```

The file is saved in the current directory as `<UniProtID>.pdb` or `<UniProtID>.cif`.

## Notes

- Structures are fetched from the [AlphaFold EBI database](https://alphafold.ebi.ac.uk/) (model v4).
- Not every UniProt ID has an AlphaFold entry — the script will tell you if one is not found.
- For multi-domain proteins, AlphaFold may have multiple fragments (F1, F2, …). This script only fetches fragment 1 (F1).
