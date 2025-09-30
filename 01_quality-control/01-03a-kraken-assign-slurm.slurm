#!/bin/bash -e

#SBATCH --account       aut03725
#SBATCH --job-name      kraken-archaea
#SBATCH --time          00:15:00
#SBATCH --mem           45G             
#SBATCH --cpus-per-task 8
#SBATCH --out           slurmlog/%j.out

# Load module
module purge
module load Kraken2/2.1.2-GCC-11.3.0
module load BLAST/2.13.0-GCC-11.3.0

# Set up variables
dirname=$(basename "$PWD")
prefix=kraken2-oomycota
results="${dirname}-${prefix}"

# Make results directory
if [ ! -d "${results}" ]; then
mkdir ${results}
fi

# Path to Kraken2 database
#export KRAKEN2_DEFAULT_DB=/nesi/project/aut03725/databases/mykraken-202309/archaea
export KRAKEN2_DEFAULT_DB=/nesi/nobackup/aut03725/ncbi-databases/oomycota-db-2023-10-31

for FILE in $(ls *_R1_001.humanFilt.fastq.gz | sed 's/_R1_001.humanFilt.fastq.gz//'); do
    kraken2 --db ${KRAKEN2_DEFAULT_DB} \
   --threads ${SLURM_CPUS_PER_TASK} \
   --use-names  \
   --report-zero-counts \
   --report ${results}/${FILE}.krakenoomycota_report.txt \
   --output ${results}/${FILE}.krakenoomycota_output.txt \
   --paired \
   ${FILE}_R1_001.humanFilt.fastq.gz ${FILE}_R2_001.humanFilt.fastq.gz
done

echo "kraken2 done"

## End of Script ##