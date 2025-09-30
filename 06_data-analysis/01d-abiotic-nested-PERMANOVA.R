## ---------------------------
##
## Script name: abiotic_nested_PERMANOVA.R
##
## Purpose of script: Nested PERMANOVA of physicochemical 
## variables across plots and sites
##
## Author: Zoe King
##
## Date Created: 2025-07-26
##
## Email: zoe.s.king@autuni.ac.nz OR zasking@gmail.com
##
## ---------------------------
##
## Notes:
## In metadata file plot = "Asymptomatic" and "Symptomatic",
## location = "Uninfected-Cascades", "Infected-Cascades", etc.
##
## ---------------------------

#****************
# Packages----
#****************
library(tidyverse)
library(vegan)
library(pairwiseAdonis)

#***********
# Data----
#***********
meta <- read_csv("../../Master-data-file-shotgun.csv") %>% 
  filter(type == "Tree") %>% 
  filter(tree_code != "UC-D2") %>% 
  mutate(health_status = case_when(
    tree_health == "healthy" ~ "Healthy",
    tree_health == "defoliated" ~ "Unhealthy",
    tree_health == "dead" ~ "Unhealthy"
  ))

# Select physicochemical data
physico <- meta %>% 
  select(tree_code, pH, electrical_conductivity, total_carbon, total_hydrogen, total_nitrogen, carbon_to_nitrogen, bulk_density, moisture_factor, water_holding_capacity) %>% 
  column_to_rownames(var = "tree_code")

#*********************
# Nested PERMANOVA----
#*********************
# Scale the data
physico_scale <- scale(physico) 
# Check rows are in the same order
table(rownames(physico) == meta$tree_code) # TRUE

# Nested permutation strategy
# Context of samples:
  # 3 sites
  # 2 plots nested within each site
  # 16 trees sampled per plot (1 sample per tree)
# Permutations must be structured to respect the nested design while allowing 
# samples to be shuffled meaningfully 
# Permutations should be blocked at the site level to restrict swapping samples across sites
# Allow free permutations within each site (shuffle samples freely but only among trees/plots within the same site)
# Modle nested structure in the formula as site/location to test site and plot effects (location is the name of the column in the metadata that contains plot information - this is confusing but ¯\_(ツ)_/¯)

set.seed(413)
# Euclidean distance matrix
physico_eu_dist <- vegdist(physico_scale, method = "euclidean")

# Define permutation scheme
cntrl <- how(within = Within(type = "free"), blocks = meta$site, nperm = 999)

# Run adonis2
adonis2(physico_eu_dist ~ site/location, data = meta,
        permutations = cntrl, by = "terms")

# PERMANOVA tree health
adonis2(physico_eu_dist ~ health_status, data = meta)

# Pairwise PERMANOVA
# Between sites
pairwise.adonis(physico_eu_dist, meta$site, p.adjust.m = "holm")

# Between plots within each site
# Loop over each Site
unique_sites <- unique(meta$site)

for (site in unique_sites) {
  cat("Running for Site:", site, "\n")
  
  # Subset metadata and distance matrix
  idx <- meta$site == site
  sub_meta <- meta[idx, ]
  sub_dist <- as.matrix(physico_eu_dist)[idx, idx]
  
  # Convert to 'dist' class
  sub_dist <- as.dist(sub_dist)
  
  # Run pairwise PERMANOVA on Plot within this Site
  pw <- pairwise.adonis(sub_dist, sub_meta$location, p.adjust.m = "holm")
  
  print(pw)
}

## End of Script ##