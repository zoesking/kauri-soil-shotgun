#!/bin/bash -e

#SBATCH --account       aut03725
#SBATCH --job-name      megahit_assembly
#SBATCH --time          24:00:00
#SBATCH --mem           40GB
#SBATCH --cpus-per-task 16
#SBATCH --output        slurmlog/%j.out

# Load modules
module purge
module load MEGAHIT/1.2.9-gimkl-2022a-Python-3.10.5

# Extract the name of current working directory (Not the full path)
current_dir=${PWD##*/}

# Remove the suffix -trim-phix-trim-human-filt-kraken2-bact
ID=$(basename ${current_dir} -trim-phix-trim-human-filt-kraken2-bact)

# Run MegaHit

for forward in *R1*fastq.gz
do
    reverse=${forward/R1/R2}
    echo "R1 sample is ${forward}"
    echo "R2 sample is ${reverse}"

        megahit -1 ${forward} -2 ${reverse} \
                -o ${ID}-megahit-assembly -t ${SLURM_CPUS_PER_TASK} --min-contig-len 1000 \
        	--tmp-dir /nesi/nobackup/aut03725/qiime-temp/temp
done
