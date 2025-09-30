#!/bin/bash -e

#SBATCH --account       aut03725
#SBATCH --job-name      map-index
#SBATCH --time          04:00:00
#SBATCH --mem           25G
#SBATCH --cpus-per-task 4
#SBATCH --out           slurmlog/%j.out

# Load modules
module purge
module load Bowtie2/2.4.5-GCC-11.3.0

# Create results directory
dirname=$(basename "$PWD")
prefix=cds-map-index
results="${dirname}-${prefix}"

if [ ! -d "${results}" ]; then
mkdir ${results}
fi

# Build index
bowtie2-build non-redundant-cds.fna ${results}/cds-map-index
