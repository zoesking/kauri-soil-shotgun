## ---------------------------
##
## Script name: gene-sets.R
##
## Purpose of script: Generate list of gene sets across health
## status groups (healthy, defoliated, and dead) for the three
## annotation sets (eggNOG, NCyc, CAZy), in preparation for
## plotting Venn diagrams.
##
## Author: Zoe King
##
## Date Created: 2025-09-19
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

#***********
# Data----
#***********
### TPM gene tables====
# This file is 7.17 GB, make sure you have enough RAM to load! 
# egg <- vroom::vroom("gene-catalog-data/tpm-tables/eggnog-coverage-tpm.tsv") %>% 
#   rename_with(~ str_remove(.x, "_.*")) %>% 
#   column_to_rownames(var = "Contig")
ncyc <- read_tsv("gene-catalog-data/tpm-tables/ncyc-coverage-tpm.tsv") %>% 
  rename_with(~ str_remove(.x, "_.*")) %>% 
  column_to_rownames(var = "Contig")
cazy <- vroom::vroom("gene-catalog-data/tpm-tables/cazy-coverage-tpm.tsv") %>% 
  rename_with(~ str_remove(.x, "_.*")) %>% 
  column_to_rownames(var = "Contig")
### Metadata====
meta <- read_csv("Master-data-file-shotgun.csv") %>%  # contains sample, tree_health, canopy_score
  filter(type == "Tree") %>% 
  filter(tree_code != "UC-D2") %>% 
  dplyr::select(tree_code, canopy_score, tree_health, location, site) %>% 
  mutate(tree_health = recode(
    tree_health,
    "healthy" = "Healthy",
    "defoliated" = "Defoliated",
    "dead" = "Dead",
  ), health_status = case_when(
    tree_health == "Healthy" ~ "Healthy",
    tree_health == "Defoliated" ~ "Unhealthy",
    tree_health == "Dead" ~ "Unhealthy"
  ))
# Get unique health groups
health_groups <- unique(meta$tree_health)

#*********************
# eggNOG gene sets----
#*********************
# Initialize a list to hold gene sets per group
gene_sets <- list()

# Loop through health groups and collect genes with TPM > 0
for (group in health_groups) {
  samples_in_group <- meta %>%
    filter(tree_health == group) %>%
    pull(tree_code)
  
  # Subset the TPM matrix
  group_tpm <- egg[, colnames(egg) %in% samples_in_group, drop = FALSE]
  
  # Find genes present in at least one sample in the group
  detected_genes <- rownames(group_tpm)[rowSums(group_tpm > 0) > 0]
  
  # Store in list
  gene_sets[[group]] <- detected_genes
}

saveRDS(gene_sets, file = "gene-catalog-data/gene-stats/eggnog-gene-sets-venn.rds")

#*********************
# NCyc gene stats----
#*********************
# Initialize a list to hold gene sets per group
gene_sets <- list()

# Loop through health groups and collect genes with TPM > 0
for (group in health_groups) {
  samples_in_group <- meta %>%
    filter(tree_health == group) %>%
    pull(tree_code)
  
  # Subset the TPM matrix
  group_tpm <- ncyc[, colnames(ncyc) %in% samples_in_group, drop = FALSE]
  
  # Find genes present in at least one sample in the group
  detected_genes <- rownames(group_tpm)[rowSums(group_tpm > 0) > 0]
  
  # Store in list
  gene_sets[[group]] <- detected_genes
}

saveRDS(gene_sets, file = "gene-catalog-data/gene-stats/ncyc-gene-sets-venn.rds")

#*********************
# CAZy gene sets----
#*********************
# Initialize a list to hold gene sets per group
gene_sets <- list()

# Loop through health groups and collect genes with TPM > 0
for (group in health_groups) {
  samples_in_group <- meta %>%
    filter(tree_health == group) %>%
    pull(tree_code)
  
  # Subset the TPM matrix
  group_tpm <- cazy[, colnames(cazy) %in% samples_in_group, drop = FALSE]
  
  # Find genes present in at least one sample in the group
  detected_genes <- rownames(group_tpm)[rowSums(group_tpm > 0) > 0]
  
  # Store in list
  gene_sets[[group]] <- detected_genes
}

saveRDS(gene_sets, file = "gene-catalog-data/gene-stats/cazy-gene-sets-venn.rds")

## End of Script ##