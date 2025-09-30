#!/bin/bash -e 

#SBATCH --account       aut03725
#SBATCH --job-name      kraken-custom-db
#SBATCH --cpus-per-task 12
#SBATCH --mem 	        12G
#SBATCH --time          1-00:00:00
#SBATCH --out           slurmlog/%j.out

module purge
module load Kraken2/2.1.2-GCC-11.3.0
module load BLAST/2.13.0-GCC-11.3.0

#define the path to the databases in the same way you have done for 
DB=/nesi/nobackup/aut03725/ncbi-databases/oomycota-db-2023-10-31

#run kraken2 to create custom database
find ./data/ -name '*.fna' -print0 | xargs -0 -I{} -n1 kraken2-build --add-to-library {} --db $DB --threads ${SLURM_CPUS_PER_TASK}
