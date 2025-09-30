## ---------------------------
##
## Script name: gene-var-part-prep.R
##
## Purpose of script: Run variation partitioning analysis
## and prep data for plotting.
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
### Distance matrix====
egg <- readRDS("rds-objects/eggnog-bc-distance-matrix.rds")
ncyc <- readRDS("rds-objects/ncyc-bc-distance-matrix.rds")
cazy <- readRDS("rds-objects/cazy-bc-distance-matrix.rds")
# Check that rows are in the same order
identical(labels(egg), labels(ncyc)) #TRUE
identical(labels(ncyc), labels(cazy)) #TRUE
### Metadata====
meta <- read_csv("../../Master-data-file-shotgun.csv") %>%  # contains sample, tree_health, canopy_score
  filter(type == "Tree") %>% 
  filter(tree_code != "UC-D2") %>% 
  rename("sample_id" = "tree_code") %>% 
  mutate(tree_health = recode(
    tree_health,
    "healthy" = "Healthy",
    "defoliated" = "Defoliated",
    "dead" = "Dead",
  ), health_status = case_when(
    tree_health == "Healthy" ~ "Healthy",
    tree_health == "Defoliated" ~ "Unhealthy",
    tree_health == "Dead" ~ "Unhealthy"
  )) %>% 
  column_to_rownames(var = "sample_id")

# Subset metadata and standardise if necessary
# Chemistry data (all)
chem <- meta[,c("pH", "electrical_conductivity", "total_carbon", "total_hydrogen", "total_nitrogen", "carbon_to_nitrogen", "bulk_density", "moisture_factor", "water_holding_capacity")] # Select variables to subset metadata dataframe
chem_std <- decostand(chem, method = "standardize") # Standardise data
chem_std <- chem_std[order(match(rownames(chem_std), rownames(as.matrix(ncyc)))), ] # Re-order rows
# Check rownames match between dataframes
if (!all(rownames(chem_std) == rownames(as.matrix(ncyc)))) {
  stop("Row names are not in the same order.")
}

### Dummy variables====
# Site and plot
plot_dummy <- meta %>%
  rownames_to_column(var = "SampleID") %>% 
  dplyr::select(location, SampleID) %>% 
  mutate(plot_num = recode(
    location,
    "Infected-Cascades" = "C_Plot_1",
    "Uinfected-Cascades" = "C_Plot_2",
    "Infected-Huia" = "H_Plot_1",
    "Uninfected-Huia" = "H_Plot_2",
    "Infected-Piha" = "P_Plot_1",
    "Uninfected-Piha" = "P_Plot_2"
  )) %>% 
  dplyr::select(-location) %>% 
  mutate(dummy = 1) %>% 
  pivot_wider(names_from = plot_num, values_from = dummy, values_fill = 0) %>% 
  dplyr::select(-last_col()) %>% 
  data.frame(row.names = "SampleID")

# Check rownames match
plot_dummy <- plot_dummy[order(match(rownames(plot_dummy), rownames(as.matrix(ncyc)))), ]
if (!all(rownames(plot_dummy) == rownames(as.matrix(ncyc)))) {
  stop("Row names are not in the same order.")
}

# Tree health
health <- meta %>%
  mutate(health_dummy = if_else(health_status == "Unhealthy", 1, 0)) %>% 
  dplyr::select(health_dummy)
# Check rownames match
health <- health[order(match(rownames(health), rownames(as.matrix(ncyc)))), ]
if (!all(rownames(health) == rownames(as.matrix(ncyc)))) {
  stop("Row names are not in the same order.")
}

#************************************
# Variation Partitioning analysis----
#************************************
varvenn_env_3_egg <- varpart(egg, chem_std, plot_dummy, health)
saveRDS(varvenn_env_3_egg, "rds-objects/eggnog-bray-varpart.rds")

varvenn_env_3_ncyc <- varpart(ncyc, chem_std, plot_dummy, health)
saveRDS(varvenn_env_3_ncyc, "rds-objects/ncyc-bray-varpart.rds")

varvenn_env_3_cazy <- varpart(cazy, chem_std, plot_dummy, health)
saveRDS(varvenn_env_3_cazy, "rds-objects/cazy-bray-varpart.rds")

# Quick view
# plot(varvenn_env_3_ncyc)

## End of Script ##