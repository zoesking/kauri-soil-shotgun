## ---------------------------
##
## Script name: gene-richness-df.R
##
## Purpose of script: Generate richness dataframes
##
## Author: Zoe King
##
## Date Created: 2025-09-20
##
## Email: zoe.s.king@autuni.ac.nz OR zasking@gmail.com
##
## ---------------------------
##
## Notes:
## The eggNOG TPM table is big!! 7.17 GB!!
## Do not try to load on your local machine 
## unless you have enough RAM or it will be unhappy :(
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
### TPM gene count tables====
# The eggNOG table is 7.17 GB! Make sure you have enough RAM before loading it in
# egg <- vroom::vroom("gene-catalog-data/tpm-tables/eggnog-coverage-tpm.tsv") %>% 
#   filter(!COG_category %in% c("-", "X")) %>% 
#   column_to_rownames(var = "COG_category") %>% 
#   t() # samples are rows, taxa are columns
ncyc <- read_tsv("gene-catalog-data/tpm-tables/ncyc-coverage-tpm.tsv") %>% 
  rename_with(~ str_remove(.x, "_.*")) %>% 
  column_to_rownames(var = "Contig") %>% 
  t() # samples are rows, taxa are columns
cazy <- vroom::vroom("gene-catalog-data/tpm-tables/cazy-coverage-tpm.tsv") %>% 
  rename_with(~ str_remove(.x, "_.*")) %>% 
  column_to_rownames(var = "Contig") %>% 
  t() # samples are rows, taxa are columns

### Metadata====
meta <- read_csv("Master-data-file-shotgun.csv") %>%  # contains sample, tree_health, canopy_score
  filter(type == "Tree") %>% 
  filter(tree_code != "UC-D2") %>% 
  dplyr::select(tree_code, canopy_score, tree_health, location, site) %>% 
  rename("sample_id" = "tree_code") %>% 
  mutate(tree_health = recode(
    tree_health,
    "healthy" = "Healthy",
    "defoliated" = "Defoliated",
    "dead" = "Dead",
  ),health_status = case_when(
    tree_health == "Healthy" ~ "Healthy",
    tree_health == "Defoliated" ~ "Unhealthy",
    tree_health == "Dead" ~ "Unhealthy"
  )) %>% 
  mutate(tree_health = fct_relevel(tree_health, "Healthy", "Defoliated", "Dead"))


#*********************
# eggNOG richness----
#*********************
# Generate richness conts
egg_richness <- specnumber(egg)
egg_richness_df <- data.frame(sample_id = rownames(egg), richness = egg_richness)

# Plot box plots of richness counts
# Join with metadata
egg_combined_df <- egg_richness_df %>%
  left_join(meta, by = "sample_id")

saveRDS(egg_combined_df, file = "gene-catalog-data/R-scripts/egg-richness-df.rds")

#*********************
# NCyc richness----
#*********************
# Generate richness conts
ncyc_richness <- specnumber(ncyc)
ncyc_richness_df <- data.frame(sample_id = rownames(ncyc), richness = ncyc_richness)

# Plot box plots of richness counts
# Join with metadata
ncyc_combined_df <- ncyc_richness_df %>%
  left_join(meta, by = "sample_id")

saveRDS(ncyc_combined_df, file = "gene-catalog-data/R-scripts/ncyc-richness-df.rds")

#*********************
# CAZy richness----
#*********************
# Generate richness conts
cazy_richness <- specnumber(cazy)
cazy_richness_df <- data.frame(sample_id = rownames(cazy), richness = cazy_richness)

# Plot box plots of richness counts
# Join with metadata
cazy_combined_df <- cazy_richness_df %>%
  left_join(meta, by = "sample_id")

saveRDS(cazy_combined_df, file = "gene-catalog-data/R-scripts/cazy-richness-df.rds")

## End of Script ##