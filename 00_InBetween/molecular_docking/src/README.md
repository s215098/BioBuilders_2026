# Docking CLI Usage (For Autodock Vina)
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

# Boltz2 and Vina ENV setup

Boltz needs to be run on a GPU, so you need to use the DTU HPC.
https://www.hpc.dtu.dk/?page_id=2129

First you need to setup miniforge on the DTU HPC if using that.

https://www.hpc.dtu.dk/?page_id=3816


```bash
conda create -n boltz_vina python=3.12 -y
```

```bash
conda activate boltz_vina
```

```bash
pip install vina meeko rdkit biopython numpy pandas boltz
```

```bash
pip uninstall -y torch torchvision torchaudio
```

```bash
pip3 install torch --index-url https://download.pytorch.org/whl/cu128
```






# Boltz2 only

Run this command to automatically calculate the bounding box for the entire receptor:

```bash
python boltz_affinity.py \
  -n "name_for_run" \
  -seq "sequence_of_the_enzyme" \
  -smi "smile_string_of_ligand_or_substrate" \
```

If you run this example command for BADGE and 

```bash
python boltz_affinity.py \
  -n "Boltz2_affinity" \
  -seq "MNPQTLPVFPDLDIFSPEYACNREKYAARALRDYPLHFYKPLNMWIVSKHKDVRSALFTPQVFSSVAFGLLPPPDDIAPRVPDLYTDVHLPSMDPPEHTKLRVPVQQALLPGRLVGKDEVVRRIANELIDTFIDKGECDLLHDFSYKLALYLIVDMLGLPKERAEDYHRWSNCFFQLFTPKVPERADARFFVPMPEEVLRQIWEDLAEANDYLREVVENLDRNPGNNMLSNLLQLREPDGSRTITISANVRNALEFGAAGHDTTATLIAHLTYFVLTTPDLKDTLTEDPSLIPAAISETLRRRGSVDGLFRRTLSDVELCGQKIESGSIVYLDLTAANLDPDVFPEPETFRLNRDNIKEMVSFGYGRHVCAGQYLSRIEAKAAYEELMRRIPNMRLADGFKLEYMPSVATTVLKGLPLVWDKN" \
  -smi "CC(C)(C1=CC=C(C=C1)OCC2CO2)C3=CC=C(C=C3)OCC4CO4" \
```
You get:

```bash
[EPOXY-SCREEN] Boltz Predicted Affinity [log10(IC50)]: 0.58  
[EPOXY-SCREEN] Boltz Confidence: 0.895
```

And the structure predictions are stored in the output folder.


# Boltz_Vina


## Docking Pipeline Overview

The script executes four autonomous steps:

### 1. Boltz-2
- Fold the enzyme's 3D structure from its raw amino acid sequence with the ligand  
- Give an Affinity and confidence score

---

### 2. Clean & Target (Biopython)
Processes the Boltz2 structure by:
- Removing the ligand and water molecules  
- Creating an empty binding pocket  
- Dynamically calculating the X, Y, Z coordinates for the docking grid box  

---

### 3. Prep (Meeko / RDKit)
Prepares both enzyme and substrate for simulation:
- Use Meeko for conversion, from pdb and sdf to pdbqt (required for Vina)

---

### 4. Dock & Score (AutoDock Vina)
Performs the docking simulation:
- Runs a physics-based thermodynamic docking process  
- Positions the substrate inside the binding pocket (binding pocket found with Boltz2
- Outputs the binding affinity as Gibbs Free Energy (ΔG) in kcal/mol (lower the better)


If both scores are good, high confidence and high affinity score from boltz2 and low (negative) ΔG from Vina it is a good candidate for further analysis (simulating actual degradation).

```bash
python boltz_vina.py \
  -n "name_for_run" \
  -seq "sequence_of_the_enzyme" \
  -smi "smile_string_of_ligand_or_substrate" \
  -sdf "path/to/sdf/of/ligand"
```

Example:

```bash
python boltz_vina.py \
  -n "BADGE_Test_01" \
  -seq "MNPQTLPVFPDLDIFSPEYACNREKYAARALRDYPLHFYKPLNMWIVSKHKDVRSALFTPQVFSSVAFGLLPPPDDIAPRVPDLYTDVHLPSMDPPEHTKLRVPVQQALLPGRLVGKDEVVRRIANELIDTFIDKGECDLLHDFSYKLALYLIVDMLGLPKERAEDYHRWSNCFFQLFTPKVPERADARFFVPMPEEVLRQIWEDLAEANDYLREVVENLDRNPGNNMLSNLLQLREPDGSRTITISANVRNALEFGAAGHDTTATLIAHLTYFVLTTPDLKDTLTEDPSLIPAAISETLRRRGSVDGLFRRTLSDVELCGQKIESGSIVYLDLTAANLDPDVFPEPETFRLNRDNIKEMVSFGYGRHVCAGQYLSRIEAKAAYEELMRRIPNMRLADGFKLEYMPSVATTVLKGLPLVWDKN" \
  -smi "CC(C)(C1=CC=C(C=C1)OCC2CO2)C3=CC=C(C=C3)OCC4CO4" \
  -sdf "/work3/s225191/igem/BADGE.sdf"
```

You get:

```bash
[EPOXY-SCREEN] Boltz Predicted Affinity [log10(IC50)]: 0.58  
[EPOXY-SCREEN] Boltz Confidence: 0.895  

[EPOXY-SCREEN] Extracting active site and stripping AI ligand...  
[EPOXY-SCREEN] Preparing Receptor PDBQT (using official Meeko flags)...  
[EPOXY-SCREEN] Receptor PDBQT successfully saved to outputs/BADGE_Test_01/receptor.pdbqt  

[EPOXY-SCREEN] Minimizing custom substrate with RDKit & Meeko...  

[EPOXY-SCREEN] Executing Native Vina at coordinates [-1.16, 3.69, 5.63]...  

Computing Vina grid ... done.  
Performing docking (random seed: 1235712316) ...  

0%   10   20   30   40   50   60   70   80   90   100%  
|----|----|----|----|----|----|----|----|----|----|  
***************************************************  

mode |   affinity | dist from best mode  
     | (kcal/mol) | rmsd l.b.| rmsd u.b.  
-----+------------+----------+----------  
   1       -8.461          0          0  
   2       -8.445    0.08459      9.047  
   3       -8.402      1.337      1.957  
   4       -8.397      1.314      8.571  
   5       -8.102      3.967       7.25  

==================================================  
 PIPELINE SUCCESS  
 Target: BADGE_Test_01  
 Vina Binding Affinity: -8.46 kcal/mol  
 Best pose saved to: outputs/BADGE_Test_01/final_poses.pdbqt  
==================================================
```

You can then visualize the final output.



