#!/bin/bash -e 

#write directory names matching "LIC" prefix to a file 
find . -type d -name "LIC*trim" > dir-list-phix.txt

#copy the slurm script to invididual directories
cat dir-list-phix.txt | xargs -n 1 -I \{} echo cp trim-phix-slurm.slurm \{}/ > copy-slurm.sh

#add executable permissions to copy script and execute using relative path 
chmod +x copy-slurm.sh && ./copy-slurm.sh 


#execute `sbatch` command for each input directory'
for dirnames in $(cat dir-list-phix.txt)
do
 cd $dirnames || continue 
 sbatch trim-phix-slurm.slurm
 cd "$OLDPWD"
done

