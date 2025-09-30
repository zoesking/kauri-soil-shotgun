## ---------------------------
##
## Script name: combine-stats-tables.R
##
## Purpose of script: Combine and view stats on QC stage 
## of metagenome pipeline
##
## Author: Zoe King
##
## Date Created: 2025-02-28
##
## Email: zoe.s.king@autuni.ac.nz OR zasking@gmail.com
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------


#****************
# Packages----
#****************
library(tidyverse)


#***********
# Data----
#***********
raw <- read_tsv("metagenome-96-raw-stats.txt")
trim <- read_tsv("trim-stats.txt")
filt <- read_tsv("filt-excl-all-stats.txt")

#*****************************************
# Remove duplicated sample from raw df----
#*****************************************
grep('IP-D49*', raw$file)
raw$file[47:50]

raw2 <- raw[-c(47:48),]
grep('IP-D49*', raw2$file)
raw2$file[47:48]

raw_files <- str_replace(raw2$file, "_subset.fastq.gz", ".fastq.gz")
raw2$file <- raw_files

#****************
# Tidy data----
#****************
raw_tidy <- raw2 %>% 
  rename_with( .fn = function(.x){paste0("raw_reads_", .x)}) %>% 
  select(raw_reads_file, raw_reads_num_seqs, raw_reads_sum_len) %>% 
  mutate(sample = str_replace(raw2$file, "_001.fastq.gz", "")) %>% # Need to create a column that is consistent across the three dataframes. This removes the specific QC step sample suffix.
  select(-raw_reads_file)

trim_tidy <- trim %>% 
  rename_with( .fn = function(.x){paste0("trim_reads_", .x)}) %>% 
  select(trim_reads_file, trim_reads_num_seqs, trim_reads_sum_len) %>% 
  mutate(sample = str_replace(trim$file, "_001.phixtrim.fastq.gz", "")) %>% 
  select(-trim_reads_file)
  
filt_tidy <- filt %>% 
  rename_with( .fn = function(.x){paste0("filt_reads_", .x)}) %>% 
  select(filt_reads_file, filt_reads_num_seqs, filt_reads_sum_len) %>% 
  mutate(sample = str_replace(filt$file, ".kraken-excl.fastq.gz", "")) %>% 
  select(-filt_reads_file)

#****************
# Join tables----
#****************
all <- right_join(raw_tidy, trim_tidy, by = "sample") %>% # Right join only matches values from "trim_tidy" (removes all unwanted samples from raw_tidy)
  left_join(., filt_tidy, by = "sample") %>% 
  select(sample, raw_reads_num_seqs, raw_reads_sum_len, trim_reads_num_seqs, trim_reads_sum_len, filt_reads_num_seqs, filt_reads_sum_len)

# Calculate percentage of reads that remain after each major QC step
all_perc <- all %>% 
  mutate(trim_perc_reads_remain = round(rowSums(all[,4]/all[,2])*100,1)) %>% 
  mutate(filt_perc_reads_remain = round(rowSums(all[,6]/all[,2])*100,1)) %>% 
  select(sample, raw_reads_num_seqs, raw_reads_sum_len, trim_reads_num_seqs, trim_reads_sum_len, trim_perc_reads_remain, filt_reads_num_seqs, filt_reads_sum_len, filt_perc_reads_remain)

# Quick look at table
summarise(all_perc)
View(all_perc)

write.table(all_perc, file = "all-qc-stats.tsv", sep = "\t", col.names = T, row.names = F)

## End of Script ##