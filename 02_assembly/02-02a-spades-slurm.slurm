#!/bin/bash -e

#SBATCH --account       aut03725      
#SBATCH --job-name      spades-assembly      
#SBATCH --time          24:00:00          
#SBATCH --mem           70G           
#SBATCH --cpus-per-task 16
#SBATCH --output        slurmlog/%j.out
#SBATCH --partition     hugemem,milan

# Load modules
module purge
module load SPAdes/3.15.4-gimkl-2022a-Python-3.10.5

# Extract the name of current working direcory (Not the full path)
current_dir=${PWD##*/}

# Remove the suffix -trim-phix-trim-human-filt-kraken2-bact
ID=$(basename ${current_dir} -trim-phix-trim-human-filt-kraken2-bact)

# Run SPAdes
for forward in *R1*fastq.gz
do
    reverse=${forward/R1/R2}
    echo "R1 sample is ${forward}"
    echo "R2 sample is ${reverse}"

    spades.py --meta -k 33,55,77,99,121 -t $SLURM_CPUS_PER_TASK \
          -1 ${forward} -2 ${reverse} \
          -o ${ID}-spades-bact/

done
