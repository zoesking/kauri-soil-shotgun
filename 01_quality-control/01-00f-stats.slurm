#!/bin/bash -e

#SBATCH --account       aut03725
#SBATCH --job-name      stats
#SBATCH --time          04:00:00
#SBATCH --mem           1G
#SBATCH --cpus-per-task 10
#SBATCH --output        slurmlog/%j.out

# Load modules
module purge
module load SeqKit/2.4.0

find -name *.fastq.gz | xargs seqkit stats -T -a -b --threads $SLURM_CPUS_PER_TASK >> stats.txt