#!/bin/bash -e 

#SBATCH --account       aut03725
#SBATCH --job-name      human-filt-bbmap
#SBATCH --cpus-per-task 12
#SBATCH --mem           30G 
#SBATCH --time          00:40:00
#SBATCH --out           slurmlog/%j.out
#SBATCH --hint          nomultithread


module purge
module load BBMap/39.01-GCC-11.3.0

dirname=$(basename "$PWD")
prefix=human-filt
results="${dirname}-${prefix}"
mkdir "$results"


if [ ! -d "${results}" ]; then
mkdir ${results}
fi

export _JAVA_OPTIONS=-Djava.io.tmpdir=/nesi/nobackup/aut03725/qiime-temp

for file in *_R1_001.phixtrim.fastq.gz
do
    base=$(basename ${file} _R1_001.phixtrim.fastq.gz)
    
    srun bbmap.sh -Xmx27g -t=${SLURM_CPUS_PER_TASK} overwrite=t\
            minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 qtrim=rl trimq=10 untrim \
            in1=${base}_R1_001.phixtrim.fastq.gz \
            in2=${base}_R2_001.phixtrim.fastq.gz \
            path=/nesi/project/aut03725/reference/ \
            out1=${results}/${base}_R1_001.humanFilt.fastq.gz \
            out2=${results}/${base}_R2_001.humanFilt.fastq.gz
done
