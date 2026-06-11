# plot_normalized.py

Plots normalized solvent accessibility (ASA) values per residue as a bar chart, colored by secondary structure type. Takes the output from `normalization.py` and the original STRIDE file as input.

---

## Requirements

- Python 3.9+
- matplotlib
- pandas

Install dependencies:
```bash
pip install matplotlib pandas
```

---

## Usage

```bash
python plot_normalized.py --input <normalized.csv> --stride <stride_output.txt> [options]
```

### Arguments

| Argument | Required | Default | Description |
|---|---|---|---|
| `--input` | Yes | — | Path to normalized ASA `.csv` file (output from `normalization.py`) |
| `--stride` | Yes | — | Path to original STRIDE output file (used to extract secondary structure) |
| `--output` | No | `normalized.png` | Output plot filename |

---

## Example

```bash
python plot_normalized.py --input 6EKZ_normalized.csv --stride 6EKZ_out.txt --output 6EKZ_plot.png
```

---

## Output

A `.png` bar chart with one bar per residue showing:

- Bar height — normalized ASA value (0 = buried, 1 = fully exposed)
- Bar color — secondary structure type (see legend below)
- Dashed line at 1.0 — reference for fully exposed

### Color legend

| Color | Secondary structure |
|---|---|
| Red | Alpha helix (H) |
| Orange | 310 helix (G) |
| Blue | Strand (E) |
| Light blue | Bridge (B) |
| Yellow | Turn (T) |
| Gray | Coil (C) |

Only secondary structure types present in the data are shown in the legend.

---

## Notes

- The `--input` file must contain `ResNum` and `Normalized_ASA` columns — these are produced automatically by `normalization.py`
- Values slightly above 1.0 are normal and expected — reference normalization values are averages
- The script exits with an informative error if the input file is missing required columns or if no ASG lines are found in the STRIDE file
- If running in the base conda environment on Mac, make sure matplotlib is installed there or activate the correct environment first