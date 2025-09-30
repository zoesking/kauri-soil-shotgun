#!/bin/bash -e 

#SBATCH --account 	aut03725
#SBATCH --job-name 	fastqc-batch
#SBATCH --cpus-per-task 36
#SBATCH --mem 		36G
#SBATCH --time 		24:00:00
#SBATCH --out 		slurmlog/%j.out 

module purge
module load FastQC/0.11.9


find -name '*.gz' |  xargs fastqc --threads ${SLURM_CPUS_PER_TASK}
