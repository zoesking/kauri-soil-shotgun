#!/bin/bash -e

#SBATCH --account       aut03725
#SBATCH --job-name      eggNOG-mapper
#SBATCH --time          3-00:00:00
#SBATCH --mem           65GB
#SBATCH --cpus-per-task 14
#SBATCH --output        slurmlog/%j.out

# Create results directory
dirname=$(basename "$PWD")
prefix=eggNOG-fast
results="${dirname}-${prefix}"

if [ ! -d "${results}" ]; then
mkdir ${results}
fi

# Load modules
module purge
module load eggnog-mapper/2.1.12-gimkl-2022a

# Run eggNOG mapper
emapper.py -i non-redundant-cds.fna --itype CDS --translate \
	 -o non-redundant-cds --output_dir ${results} \
	--cpu $SLURM_CPUS_PER_TASK \
	--data_dir $DATA_PATH \
	-m diamond \
	--dbmem \
	--temp_dir /nesi/nobackup/aut03725/qiime-temp/temp