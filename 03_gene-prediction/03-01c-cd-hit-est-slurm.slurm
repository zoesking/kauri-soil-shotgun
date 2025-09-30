#!/bin/bash -e

#SBATCH --account       aut03725
#SBATCH --job-name      cd-hit-est
#SBATCH --time          5-00:00:00
#SBATCH --mem           60GB
#SBATCH --cpus-per-task 64
#SBATCH --output        slurmlog/%j.out
#SBATCH --mail-type	TIME_LIMIT_80
#SBATCH --mail-user	zoe.s.king@autuni.ac.nz

# Load modules
module purge
module load CD-HIT/4.8.1-GCC-11.3.0

# Run CD-Hit-EST
cd-hit-est -i all-cds-no-meta-rename.fna -o non-redundant-cds.fna -c 0.95 -T $SLURM_CPUS_PER_TASK -M $SLURM_MEM_PER_NODE -G 0 -aS 0.9 -g 1 -r 1 -d 0
