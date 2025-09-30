#!/bin/bash -e 

#SBATCH --account       aut03725
#SBATCH --job-name      bbduk-trimming
#SBATCH --cpus-per-task 8
#SBATCH --mem           6G 
#SBATCH --time          00:15:00
#SBATCH --out           slurmlog/%j.out
#SBATCH --hint          nomultithread

module purge
module load BBMap/39.01-GCC-11.3.0

dirname=$(basename "$PWD")
prefix=trim
results="${dirname}-${prefix}"
mkdir "$results"

bbmap_adapters=/nesi/project/aut03725/adapters/bbmap-adapters.fa

if [ ! -d "${results}" ]; then
mkdir ${results}
fi


for file in *_R1_001.fastq.gz
do
    base=$(basename ${file} _R1_001.fastq.gz)
    
    srun bbduk.sh qtrim=rl trimq=25 \
    in1=${base}_R1_001.fastq.gz in2=${base}_R2_001.fastq.gz \
    out1=${results}/${base}_R1_001_trim.fastq.gz out2=${results}/${base}_R2_001_trim.fastq.gz \
    ref=${bbmap_adapters} \
    threads=${SLURM_CPUS_PER_TASK}
done
