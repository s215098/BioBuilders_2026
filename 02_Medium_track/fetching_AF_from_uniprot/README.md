# fetch_alphafold.py

Downloads AlphaFold PDB structures for a given UniProt ID, using the names listed in UniProt's structure section.

## How it works

1. Queries the **UniProt REST API** for the entry and finds all cross-references where `source = "AlphaFold DB"`.
2. For each AlphaFold entry ID found, calls the **AlphaFold EBI API** to get the direct PDB download URL and the structure name (`uniprotDescription`).
3. Downloads each PDB file, named after the structure (e.g. `Hemoglobin_subunit_beta.pdb`).

## Requirements

- Python 3 (standard library only — no `pip install` needed)

## Usage

```bash
python fetch_alphafold.py <UniProtID>
```

## Examples

```bash
# Hemoglobin subunit beta — single fragment
python fetch_alphafold.py P68871

# Titin — multiple fragments (F1, F2, ...)
python fetch_alphafold.py Q8WZ42
```

### Example output

```
Looking up AlphaFold structures for UniProt ID: P68871
Found 1 AlphaFold entry/entries: P68871

Fetching metadata for P68871 ...
  Name   : Hemoglobin subunit beta
  URL    : https://alphafold.ebi.ac.uk/files/AF-P68871-F1-model_v4.pdb
  Saving : Hemoglobin subunit beta.pdb
  Done.
```

## Notes

- Files are saved in the **current working directory**.
- For large proteins with multiple fragments (e.g. Titin), all fragments are downloaded.
- If a UniProt ID has no AlphaFold entry, the script exits with a clear message.