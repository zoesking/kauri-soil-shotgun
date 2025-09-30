#!/bin/bash -e

#SBATCH --account       aut03725
#SBATCH --job-name      pigz-compress
#SBATCH --time          5-00:00:00
#SBATCH --mem           2G             
#SBATCH --cpus-per-task 12
#SBATCH --out           slurmlog/%j.out

# Load modules
module purge
module load pigz/2.7

# Run compression

find -name '*.fastq' | xargs pigz -v -p ${SLURM_CPUS_PER_TASK}

## End of Script ##