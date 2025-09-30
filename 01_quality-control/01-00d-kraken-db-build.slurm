#!/bin/bash -e 

#SBATCH --account 	aut03725
#SBATCH --job-name	kraken2-build
#SBATCH --cpus-per-task 12
#SBATCH --mem 		1700G
#SBATCH --time		7-00:00:00
#SBATCH --out		slurmlog/%j.out


module purge
module load Kraken2/2.1.2-GCC-11.3.0
module load BLAST/2.13.0-GCC-11.3.0

DB=/nesi/nobackup/aut03725/ncbi-databases/arthropoda-db-2023-10-31

kraken2-build --threads ${SLURM_CPUS_PER_TASK} --build --db $DB
