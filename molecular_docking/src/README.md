# Docking CLI Usage
To be run from the src folder. Or you would need to provide where the script is.
## Blind Docking

Run this command to automatically calculate the bounding box for the entire receptor:

```bash
python docking.py \
  --receptor path/to/receptor.pdbqt \
  --ligand path/to/ligand.sdf \
  --out_dir path/to/output_dir/
```

Targeted Docking

Run this command to restrict the search space to a specific active site using known coordinates:

```bash
python docking.py \
  --receptor path/to/receptor.pdbqt \
  --ligand path/to/ligand.sdf \
  --out_dir path/to/output_dir/ \
  --center 18.8 32.2 14.9 \
  --box_size 15.0 15.0 15.0
```


Env:

```bash
# Create the environment
conda create -n bio_docking python=3.12 -y
conda activate bio_docking

# Install build dependencies via conda-forge
conda install -c conda-forge swig boost-cpp -y

pip install "numpy<2.0" 

# Install the rest of your specific versions
# Note: RDKit 2025+ and Meeko 0.7+ are happy with Python 3.12
export BOOST_INCLUDE=$CONDA_PREFIX/include
export BOOST_LIB=$CONDA_PREFIX/lib

pip install \
  vina==1.2.7 \
  meeko==0.7.1 \
  rdkit==2025.9.5 \
  py3Dmol==2.5.4 \
  stmol==0.0.9 \
  ipython-genutils==0.2.0 \
  scipy==1.16.1 \
  gemmi==0.7.4 \
  prody==2.6.1
```
