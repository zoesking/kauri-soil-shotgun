#!/bin/bash -e

#SBATCH --account       aut03725
#SBATCH --job-name      host_filt_bbmap_index
#SBATCH --partition     milan
#SBATCH --time          00:20:00
#SBATCH --mem           32GB
#SBATCH --cpus-per-task 12
#SBATCH --output        %x_%j.out

# Working directory
#/nesi/project/aut03725/reference

# Load BBMap module
module purge
module load BBMap/39.01-GCC-11.3.0

# Build indexed reference file via BBMap
bbmap.sh ref=hg19_main_mask_ribo_animal_allplant_allfungus.fa.gz  
