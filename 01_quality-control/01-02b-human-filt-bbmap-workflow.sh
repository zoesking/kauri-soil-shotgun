#!/bin/bash -e 

#write directory names matching "LIC" prefix to a file 
find . -type d -name "LIC*phix-trim" > dir-list-human.txt

#copy the slurm script to invididual directories
cat dir-list-human.txt | xargs -n 1 -I \{} echo cp human-filt-bbmap-slurm.slurm \{}/ > copy-slurm.sh

#add executable permissions to copy script and execute using relative path 
chmod +x copy-slurm.sh && ./copy-slurm.sh 


#execute `sbatch` command for each input directory'
for dirnames in $(cat dir-list-human.txt)
do
 cd $dirnames || continue 
 sbatch human-filt-bbmap-slurm.slurm
 cd "$OLDPWD"
done
