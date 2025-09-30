#!/bin/bash

#SBATCH --account       aut03725
#SBATCH --job-name      diamond-trial
#SBATCH --time          00:50:00
#SBATCH --mem           10GB
#SBATCH --cpus-per-task 6
#SBATCH --output        slurmlog/%j.out

# Load modules 
module purge
module load DIAMOND/2.1.10-GCC-12.3.0

# Output directory
dirname=$(basename "$PWD")
prefix=diamond
results="${dirname}-${prefix}"

if [ ! -d "${results}" ]; then
mkdir ${results}
fi

# Run diamond
diamond blastx \
 -d /nesi/nobackup/aut03725/databases/NCyc-db/NCyc_100.dmnd \
 --evalue 0.00001 \
 --max-target-seqs 1 \
 -q non-redundant-cds.fna \
 --out ${results}/non-redundant-cds-ncyc.diamond \
 --outfmt 6 \
 -p ${SLURM_CPUS_PER_TASK} \
 --tmpdir /nesi/nobackup/aut03725/qiime-temp/temp
