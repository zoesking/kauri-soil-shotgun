#!/bin/bash -e

#SBATCH --account       aut03725
#SBATCH --job-name      prodigal
#SBATCH --time          2:30:00
#SBATCH --mem           1G
#SBATCH --cpus-per-task 1
#SBATCH --output        slurmlog/%j.out

# Load modules
module purge
module load prodigal/2.6.3-GCCcore-7.4.0

# Output directory
dirname=$(basename "$PWD")
prefix=prodigal
results="${dirname}-${prefix}"

if [ ! -d "${results}" ]; then
mkdir ${results}
fi

# Extract sample code
for forward in ../*_R1.kraken-excl.fastq.gz;
do
sample_id=$(basename ${forward} _R1.kraken-excl.fastq.gz)
echo ${sample_id}; 
done

# Run prodigal
prodigal -i final.contigs.fa -p meta \
         -d ${results}/${sample_id}.genes.fna \
         -a ${results}/${sample_id}.genes.faa \
         -o ${results}/${sample_id}.genes.gbk
