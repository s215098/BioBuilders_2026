#!/bin/bash

#BSUB -J mutatex_6EKZ
#BSUB -n 4
#BSUB -R "rusage[mem=16GB]"
#BSUB -W 24:00
#BSUB -o mutatex_%J.log
#BSUB -e mutatex_%J.err

# Activate the shared mutatex environment
conda activate /shared/envs/mutatex-env

# Make sure you're in the specific protein folder:
cd /Users/kristinetoftjohansen/Desktop/Biobuilders_2026/02_Medium_track/Mutagenesis/MutateX/6EKZ/

# Run MutateX with the specified parameters
mutatex 6EKZ_clean.pdb \
-x ~/foldx/foldx \
-f suite5 \
-p 4 \
-m mutation_list.txt \
-q position_list.txt \
-R repair_runfile_template.txt \
-M mutate_runfile_template.txt \
-b /Users/kristinetoftjohansen/foldx/rotabase.txt