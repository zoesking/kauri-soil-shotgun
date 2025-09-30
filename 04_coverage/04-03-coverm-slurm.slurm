#!/bin/bash -e

#SBATCH --account       aut03725
#SBATCH --job-name      coverm-coverage
#SBATCH --time          05:00:00
#SBATCH --mem           300GB
#SBATCH --cpus-per-task 4
#SBATCH --output        slurmlog/%j.out

# Load modules
module purge
module load CoverM/0.7.0-GCC-12.3.0

coverm contig -t ${SLURM_CPUS_PER_TASK} -b *.bam -m tpm -o coverm-cds-tpm-cov-tab.txt
