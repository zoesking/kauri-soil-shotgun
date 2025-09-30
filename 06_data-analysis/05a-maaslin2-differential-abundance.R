## ---------------------------
##
## Script name: maaslin2-differential-abundance.R
##
## Purpose of script:
##
## Author: Zoe King
##
## Date Created: 2025-06-29
##
## Email: zoe.s.king@autuni.ac.nz OR zasking@gmail.com
##
## ---------------------------
##
## Notes:
## When running the eggNOG dataset this should be run on z6
## via a tmux session - it took 2 weeks to run.
##
## ---------------------------

#****************
# Packages----
#****************
library(tidyverse)
library(Maaslin2)

#***********
# Data----
#***********
# Using the tpm coverage information as I think this is similar to the output from HUMAnN that has been normalized with the humann_renorm_table script. https://forum.biobakery.org/t/rna-and-dna-cpm-as-covariates/1617/3
dat <- read_tsv("eggnog-coverage-tpm.tsv") %>% 
  rename_with(~ str_remove(.x, "_.*")) %>% 
  column_to_rownames(var = "Contig") %>% 
  t() # Features as columns, samples as rows

# Metadata
meta <- read_csv("../Master-data-file-shotgun.csv") %>%  # contains sample, tree_health, canopy_score
  filter(type == "Tree") %>% 
  filter(tree_code != "UC-D2") %>% 
  select(tree_code, canopy_score, tree_health, location, site) %>% 
  rename("sample_code" = "tree_code") %>% 
  mutate(tree_health = recode(
    tree_health,
    "healthy" = "Healthy",
    "defoliated" = "Defoliated",
    "dead" = "Dead",
  )) %>% 
  mutate(health_status = case_when(
    tree_health == "Healthy" ~ "Healthy",
    tree_health == "Defoliated" ~ "Unhealthy",
    tree_health == "Dead" ~ "Unhealthy"
  )) %>% 
  column_to_rownames(var = "sample_code")

#**************************
# Running Maaslin2----
#**************************
fit_data = Maaslin2(input_data     = dat, 
                    input_metadata = meta, 
                    min_prevalence = 0.1, # Prevalence filtering to only include features in >10% of samples (test only features with at least 10% non-zero values)
                    normalization  = "NONE",
                    transform      = "NONE",
                    output         = "maaslin2_eggNOG_output", 
                    fixed_effects  = "health_status")

## End of Script ##