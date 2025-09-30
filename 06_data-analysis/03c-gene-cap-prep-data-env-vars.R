## ---------------------------
##
## Script name: gene-cap-prep-data-env-vars.R
##
## Purpose of script: Prepare data for constrained ordination
## with env fit analysis
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
### Distance matrices====
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
  #select(tree_code, canopy_score, tree_health, location, site) %>% 
  column_to_rownames(var = "tree_code") %>% 
  mutate(tree_health = recode(
    tree_health,
    "healthy" = "Healthy",
    "defoliated" = "Defoliated",
    "dead" = "Dead",
  ), health_status = case_when(
    tree_health == "Healthy" ~ "Healthy",
    tree_health == "Defoliated" ~ "Unhealthy",
    tree_health == "Dead" ~ "Unhealthy"
  ), plot_num = recode(
    location,
    "Infected-Cascades" = "Plot 1",
    "Uinfected-Cascades" = "Plot 2",
    "Infected-Huia" = "Plot 1",
    "Uninfected-Huia" = "Plot 2",
    "Infected-Piha" = "Plot 1",
    "Uninfected-Piha" = "Plot 2"
  ))

### Extract physicochemical environmental data====
chem_meta <- meta %>%
  rownames_to_column(var = "SampleID") %>% 
  dplyr::select(
    SampleID,
    pH,
    electrical_conductivity,
    total_carbon,
    total_hydrogen,
    total_nitrogen,
    carbon_to_nitrogen,
    bulk_density,
    moisture_factor,
    water_holding_capacity
  ) %>%
  data.frame(row.names = "SampleID") # Need to remove sample names as column to standardise data
chem_std <- decostand(chem_meta, method = "standardize") # Standardise data

# Check rows are in the same order
chem_std <- chem_std[order(match(rownames(chem_std), rownames(as.matrix(ncyc)))), ] # Re-order rows
# Sanity check that rows are in the same order
if (!all(rownames(as.matrix(ncyc)) == rownames(chem_std))) {
  stop("Row names are not in the same order.")
}

### Set seed for analyses====
set.seed(413)

#****************
# eggNOG data----
#****************
### Constrained ordination====
cap_egg <- capscale(egg ~ site/location, data = meta)

saveRDS(cap_egg, "rds-objects/eggnog-bray-cap-ordination.rds")

# cap_egg <- readRDS("constrained-ordination/egg-cap-ordination.rds")
# summary(cap_egg)
# RsquareAdj(cap_egg)

# Get sample scores (site scores)
site_scores_egg <- scores(cap_egg, display = "sites") %>% as.data.frame()
site_scores_egg$tree_code <- rownames(site_scores_egg)

saveRDS(site_scores_egg, "rds-objects/eggnog-bray-cap-scores.rds")

### Envfit analysis====
# Prep data
# Sanity check that rows are in the same order
chem_std <- chem_std[order(match(rownames(chem_std), rownames(site_scores_egg))), ] # Re-order rows
table(rownames(chem_std) == rownames(site_scores_egg)) # Sanity check that rows are in the same order # TRUE

fit <- envfit(cap_egg, chem_std, permutations = 999)

### Prepare env data for plotting====
# Join metadata to scores
plot_data_egg <- meta %>% 
  rownames_to_column(var = "tree_code") %>% 
  left_join(site_scores_egg, by = "tree_code") %>% 
  column_to_rownames(var = "tree_code")

saveRDS(plot_data_egg, file = "rds-objects/eggnog-bray-ord-meta-df.rds")

# Apply factor to allow arrow lengths to make full use of plot area (refer to linked stackoverflow thread for more explanation)
arrow_factor <- ordiArrowMul(fit)
fit_scores <- as.data.frame(vegan::scores(fit, display = "vectors"))*arrow_factor
fit_scores <- cbind(fit_scores, envvar = rownames(fit_scores), Pvalues = fit$vectors$pvals, r_square = fit$vectors$r)
fit_scores <- subset(fit_scores, Pvalues < 0.05)

# Rename physicochemical variables
fit_scores <- fit_scores %>% 
  mutate(envname = recode(
    envvar,
    "electrical_conductivity" = "Electrical conductivity", 
    "total_carbon" = "Total carbon",
    "total_hydrogen" = "Total hydrogen",
    "total_nitrogen" = "Total nitrogen",
    "carbon_to_nitrogen" = "C:N ratio",
    "bulk_density" = "Bulk density",
    "moisture_factor" = "Soil moisture",
    "water_holding_capacity" = "Water holding capacity",
    "elevation_mean" = "Elevation",
    "longitude" = "Longitude",
    "latitude" = "Latitude",
    "x_coordinates" = "Easting",
    "y_coordinates" = "Northing"
  ))

saveRDS(fit_scores, file = "rds-objects/eggnog-bray-cap-envfit-scores.rds")

#****************
# NCyc data----
#****************
### Constrained ordination====
cap_ncyc <- capscale(ncyc ~ site/location, data = meta)

saveRDS(cap_ncyc, "rds-objects/ncyc-bray-cap-ordination.rds")

# cap_ncyc <- readRDS("constrained-ordination/ncyc-cap-ordination.rds")
# summary(cap_ncyc)
# RsquareAdj(cap_ncyc)

# Get sample scores (site scores)
site_scores_ncyc <- scores(cap_ncyc, display = "sites") %>% as.data.frame()
site_scores_ncyc$tree_code <- rownames(site_scores_ncyc)

saveRDS(site_scores_ncyc, "rds-objects/ncyc-bray-cap-scores.rds")

### Envfit analysis====
# Prep data
# Sanity check that rows are in the same order
chem_std <- chem_std[order(match(rownames(chem_std), rownames(site_scores_ncyc))), ] # Re-order rows
table(rownames(chem_std) == rownames(site_scores_ncyc)) # Sanity check that rows are in the same order # TRUE

fit <- envfit(cap_ncyc, chem_std, permutations = 999)

### Prepare env data for plotting====
# Join metadata to scores
plot_data_ncyc <- meta %>% 
  rownames_to_column(var = "tree_code") %>% 
  left_join(site_scores_ncyc, by = "tree_code") %>% 
  column_to_rownames(var = "tree_code")

saveRDS(plot_data_ncyc, file = "rds-objects/ncyc-bray-ord-meta-df.rds")

# Apply factor to allow arrow lengths to make full use of plot area (refer to linked stackoverflow thread for more explanation)
arrow_factor <- ordiArrowMul(fit)
fit_scores <- as.data.frame(vegan::scores(fit, display = "vectors"))*arrow_factor
fit_scores <- cbind(fit_scores, envvar = rownames(fit_scores), Pvalues = fit$vectors$pvals, r_square = fit$vectors$r)
fit_scores <- subset(fit_scores, Pvalues < 0.05)

# Rename physicochemical variables
fit_scores <- fit_scores %>% 
  mutate(envname = recode(
    envvar,
    "electrical_conductivity" = "Electrical conductivity", 
    "total_carbon" = "Total carbon",
    "total_hydrogen" = "Total hydrogen",
    "total_nitrogen" = "Total nitrogen",
    "carbon_to_nitrogen" = "C:N ratio",
    "bulk_density" = "Bulk density",
    "moisture_factor" = "Soil moisture",
    "water_holding_capacity" = "Water holding capacity",
    "elevation_mean" = "Elevation",
    "longitude" = "Longitude",
    "latitude" = "Latitude",
    "x_coordinates" = "Easting",
    "y_coordinates" = "Northing"
  ))

saveRDS(fit_scores, file = "rds-objects/ncyc-bray-cap-envfit-scores.rds")

#****************
# CAZy data----
#****************
### Constrained ordination====
cap_cazy <- capscale(cazy ~ site/location, data = meta)

saveRDS(cap_cazy, "rds-objects/cazy-bray-cap-ordination.rds")

# cap_cazy <- readRDS("constrained-ordination/cazy-cap-ordination.rds")
# summary(cap_cazy)
# RsquareAdj(cap_cazy)

# Get sample scores (site scores)
site_scores_cazy <- scores(cap_cazy, display = "sites") %>% as.data.frame()
site_scores_cazy$tree_code <- rownames(site_scores_cazy)

saveRDS(site_scores_cazy, "rds-objects/cazy-bray-cap-scores.rds")

### Envfit analysis====
# Prep data
# Sanity check that rows are in the same order
chem_std <- chem_std[order(match(rownames(chem_std), rownames(site_scores_cazy))), ] # Re-order rows
table(rownames(chem_std) == rownames(site_scores_cazy)) # Sanity check that rows are in the same order # TRUE

fit <- envfit(cap_cazy, chem_std, permutations = 999)

### Prepare env data for plotting====
# Join metadata to scores
plot_data_cazy <- meta %>% 
  rownames_to_column(var = "tree_code") %>% 
  left_join(site_scores_cazy, by = "tree_code") %>% 
  column_to_rownames(var = "tree_code")

saveRDS(plot_data_cazy, file = "rds-objects/cazy-bray-ord-meta-df.rds")

# Apply factor to allow arrow lengths to make full use of plot area (refer to linked stackoverflow thread for more explanation)
arrow_factor <- ordiArrowMul(fit)
fit_scores <- as.data.frame(vegan::scores(fit, display = "vectors"))*arrow_factor
fit_scores <- cbind(fit_scores, envvar = rownames(fit_scores), Pvalues = fit$vectors$pvals, r_square = fit$vectors$r)
fit_scores <- subset(fit_scores, Pvalues < 0.05)

# Rename physicochemical variables
fit_scores <- fit_scores %>% 
  mutate(envname = recode(
    envvar,
    "electrical_conductivity" = "Electrical conductivity", 
    "total_carbon" = "Total carbon",
    "total_hydrogen" = "Total hydrogen",
    "total_nitrogen" = "Total nitrogen",
    "carbon_to_nitrogen" = "C:N ratio",
    "bulk_density" = "Bulk density",
    "moisture_factor" = "Soil moisture",
    "water_holding_capacity" = "Water holding capacity",
    "elevation_mean" = "Elevation",
    "longitude" = "Longitude",
    "latitude" = "Latitude",
    "x_coordinates" = "Easting",
    "y_coordinates" = "Northing"
  ))

saveRDS(fit_scores, file = "rds-objects/cazy-bray-cap-envfit-scores.rds")

## End of Script ##