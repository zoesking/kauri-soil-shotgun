#!/bin/bash -e

#SBATCH --account       aut03725
#SBATCH --job-name      cds-map-reads
#SBATCH --time          04:00:00
#SBATCH --mem           18G
#SBATCH --cpus-per-task 10
#SBATCH --output        bowtie-slurmlog/%j.out

# Load modules
module purge
module load Bowtie2/2.4.5-GCC-11.3.0 SAMtools/1.15.1-GCC-11.3.0

# Extract the name of current working directory (Not the full path)
current_dir=${PWD##*/}

# Remove the suffix
ID=$(basename ${current_dir} -trim-phix-trim-human-filt-kraken2-all-excl)

# Create results directory
prefix=cds-map-reads
results="${ID}-${prefix}"

if [ ! -d "${results}" ]; then
mkdir ${results}
fi

# Extract sample code
for forward in *_R1.kraken-excl.fastq.gz;
do
sample_id=$(basename ${forward} _R1.kraken-excl.fastq.gz)
echo ${sample_id}; 
done

# Mapping reads
bowtie2 --minins 200 --maxins 800 --threads $SLURM_CPUS_PER_TASK --sensitive \
        -x /nesi/nobackup/aut03725/mg-outputs/filter-output/kraken2-output/exclude-all/non-redundant-catalog/non-redundant-catalog-cds-map-index/cds-map-index \
        -1 ${sample_id}_R1.kraken-excl.fastq.gz -2 ${sample_id}_R2.kraken-excl.fastq.gz \
        -S ${results}/${sample_id}-cds-map-reads.sam

# Sort and compress results
samtools sort -@ $SLURM_CPUS_PER_TASK -o ${results}/${sample_id}-cds-map-reads.bam ${results}/${sample_id}-cds-map-reads.sam


