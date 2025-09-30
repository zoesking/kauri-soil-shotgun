#!/bin/bash -e

#SBATCH --account       aut03725
#SBATCH --job-name      kraken-extract-reads
#SBATCH --time          08:00:00
#SBATCH --mem           12G 
#SBATCH --cpus-per-task 2
#SBATCH --out           slurmlog/%j.out

# Extract the name of current working direcory (Not the full path)
current_dir=${PWD##*/}

# Remove the suffix
ID=$(basename ${current_dir} -plant)

# Taxon to remove
taxa="annelida"

# Load modules
module purge
module load Python/3.11.3-gimkl-2022a

for FILE in $(ls *_R1.kraken-plant-extract.fastq.gz | sed 's/_R1.kraken-plant-extract.fastq.gz//'); do
        /nesi/nobackup/aut03725/software/extract_kraken_reads.py \
        -k /nesi/nobackup/aut03725/mg-outputs/filter-output/kraken2-output/${taxa}/${ID}-${taxa}/${FILE}.kraken${taxa}_output.txt \
        -r /nesi/nobackup/aut03725/mg-outputs/filter-output/kraken2-output/${taxa}/${ID}-${taxa}/${FILE}.kraken${taxa}_report.txt \
        -s1 ${FILE}_R1.kraken-exclu-plant-extract.fastq -s2 ${FILE}_R2.kraken-exclu-plant-extract.fastq \
        -t 6340 --include-children --exclude \
        --fastq-output -o ${FILE}_R1.kraken-ex-pl-an.fastq -o2 ${FILE}_R2.kraken-ex-pl-an.fastq
done

## End of Script ##