# Docking CLI Usage

## Blind Docking

Run this command to automatically calculate the bounding box for the entire receptor:

```bash
python docking_cli.py \
  --receptor path/to/receptor.pdbqt \
  --ligand path/to/ligand.sdf \
  --out_dir path/to/output_dir/
```

Targeted Docking

Run this command to restrict the search space to a specific active site using known coordinates:

```bash
python docking_cli.py \
  --receptor path/to/receptor.pdbqt \
  --ligand path/to/ligand.sdf \
  --out_dir path/to/output_dir/ \
  --center 18.8 32.2 14.9 \
  --box_size 15.0 15.0 15.0
```