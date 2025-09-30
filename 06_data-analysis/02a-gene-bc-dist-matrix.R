## ---------------------------
##
## Script name: gene-bc-dist-matrix.R
##
## Purpose of script: Generate Bray-Curtis distance matrix
##
## Author: Zoe King
##
## Date Created: 2025-09-30
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
library(vegan)

#***********
# Data----
#***********
tpm_tab <- vroom("eggnog-coverage-tpm.tsv") %>% 
  rename_with(~ str_remove(.x, "_.*")) %>% 
  column_to_rownames(var = "Contig")

#*********************
# Distance Matrix----
#*********************
# Bray-Curtis distance matrix
bcdist_css <- vegdist(t(tpm_tab), method = "bray") # Transpose so samples are rows

saveRDS(bcdist_css, file = "eggnog-bc-distance-matrix.rds")

## End of Script ##