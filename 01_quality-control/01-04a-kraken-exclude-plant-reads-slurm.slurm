#!/bin/bash -e

#SBATCH --account       aut03725
#SBATCH --job-name      kraken-extract-reads
#SBATCH --time          04:00:00
#SBATCH --mem           12G
#SBATCH --cpus-per-task 2
#SBATCH --out           slurmlog/%j.out

# Set up variables
dirname=$(basename "$PWD")
prefix=kraken2-plant
results="${dirname}-${prefix}"

# Load modules
module purge
module load Python/3.11.3-gimkl-2022a

for FILE in $(ls *_R1_001_humanFilt.fastq.gz | sed 's/_R1_001_humanFilt.fastq.gz//'); do
	/nesi/nobackup/aut03725/software/extract_kraken_reads.py \
	-k ${results}/${FILE}.krakenplant_output.txt \
	-r ${results}/${FILE}.krakenplant_report.txt \
	-s1 ${FILE}_R1_001_humanFilt.fastq.gz -s2 ${FILE}_R2_001_humanFilt.fastq.gz \
	-t 33090 --include-children --exclude \
	--fastq-output -o ${results}/${FILE}_R1_kraken-ex-pl.fastq -o2 ${results}/${FILE}_R2_kraken-ex-pl.fastq
done

## End of Script ##