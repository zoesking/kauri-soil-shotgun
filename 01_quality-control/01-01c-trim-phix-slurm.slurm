#!/bin/bash -e 

#SBATCH --account       aut03725
#SBATCH --job-name      bbduk-phix-trimming
#SBATCH --cpus-per-task 8
#SBATCH --mem           6G 
#SBATCH --time          00:15:00
#SBATCH --out           slurmlog/%j.out
#SBATCH --hint          nomultithread


module purge
module load BBMap/39.01-GCC-11.3.0

dirname=$(basename "$PWD")
prefix=phix-trim
results="${dirname}-${prefix}"
mkdir "$results"

phix_adapter=/nesi/project/aut03725/adapters/phix_adapters.fa

if [ ! -d "${results}" ]; then
mkdir ${results}
fi

for file in *_R1_001_trim.fastq.gz
do
        base=$(basename ${file} _R1_001_trim.fastq.gz)
        
        srun bbduk.sh in1=${base}_R1_001_trim.fastq.gz \
                       in2=${base}_R2_001_trim.fastq.gz \
                       out1=${results}/${base}_R1_001_phixtrim.fastq.gz \
                       out2=${results}/${base}_R2_001_phixtrim.fastq.gz \
                       ref=${phix_adapter} \
                       k=25 hdist=1 stats=${base}_paired_stats.txt \
                       threads=${SLURM_CPUS_PER_TASK}
done