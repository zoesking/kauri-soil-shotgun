#!/bin/bash -e

#SBATCH --account       aut03725
#SBATCH --job-name      seqtk-subset
#SBATCH --time          00:30:00
#SBATCH --mem           50G             
#SBATCH --cpus-per-task 4
#SBATCH --out           slurmlog/%j.out

# Load modules
module purge
module load seqtk/1.4-GCC-11.3.0

# Run command
# Remember to use the same random seed to keep pairing
seqtk sample -s589 NRT-PI-A26_S57_R1_001.fastq.gz 60000000 > NRT-PI-A26_S57_R1_001_subset.fastq
seqtk sample -s589 NRT-PI-A26_S57_R2_001.fastq.gz 60000000 > NRT-PI-A26_S57_R2_001_subset.fastq