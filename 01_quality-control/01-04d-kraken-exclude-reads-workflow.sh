#!/bin/bash -e 

#write directory names matching "LIC" prefix to a file 
find . -type d -name "LIC*human-filt-kraken2-plant" > dir-list.txt

#copy the slurm script to invididual directories
cat dir-list.txt | xargs -n 1 -I \{} echo cp kraken-exclude-reads-slurm.slurm \{}/ > copy-slurm.sh

#add executable permissions to copy script and execute using relative path 
chmod +x copy-slurm.sh && ./copy-slurm.sh 


#execute `sbatch` command for each input directory'
for dirnames in $(cat dir-list.txt)
do
 cd $dirnames || continue 
 sbatch kraken-exclude-reads-slurm.slurm
 cd "$OLDPWD"
done

## End of Script ##