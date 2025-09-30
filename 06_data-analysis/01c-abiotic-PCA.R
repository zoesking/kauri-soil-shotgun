## ---------------------------
##
## Script name: abiotic_PCA.R
##
## Purpose of script: PCA of soil physicochemical data
##
## Author: Zoe King
##
## Date Created: 2025-07-25
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
library(ggrepel)
library(ggpp)
library(ggnewscale) #ggnewscale: spend 400% more time tweaking your ggplot!

#***********
# Data----
#***********
# Load in metadata
meta <- read_csv("../../Master-data-file-shotgun.csv") %>% 
  filter(type == "Tree") %>% 
  filter(tree_code != "UC-D2") %>% 
  mutate(plot_num = recode(
    location,
    "Infected-Cascades" = "Plot 1",
    "Uinfected-Cascades" = "Plot 2",
    "Infected-Huia" = "Plot 1",
    "Uninfected-Huia" = "Plot 2",
    "Infected-Piha" = "Plot 1",
    "Uninfected-Piha" = "Plot 2"
  ))

# Select physicochemical data
physico <- meta %>% 
  select(tree_code, pH, electrical_conductivity, total_carbon, total_hydrogen, total_nitrogen, carbon_to_nitrogen, bulk_density, moisture_factor, water_holding_capacity) %>% 
  column_to_rownames(var = "tree_code")

#***********
# PCA----
#***********
# Run PCA with centering and scaling (gives Z-scores)
pca_result <- prcomp(physico, center = TRUE, scale. = TRUE)

# Scores (sample projections)
scores <- as.data.frame(pca_result$x)

# Loadings (variable contributions)
loadings <- as.data.frame(pca_result$rotation)

# Variance explained  
var_exp <- summary(pca_result)$importance[2, 1:2] * 100

# saveRDS(var_exp, file = "abiotic-data/abiotic-pca-var-exp.rds")

# Scaling factor to stretch arrows (adjust if needed)
arrow_scale <- 5

loadings <- loadings %>%
  mutate(Variable = rownames(loadings),
         x = PC1 * arrow_scale,
         y = PC2 * arrow_scale) %>% 
  mutate(envname = recode(
  Variable,
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

saveRDS(loadings, file = "abiotic-data/abiotic-loadings.rds")

# Combine scores and metadata
scores <- scores %>%
  mutate(Sample = meta$tree_code,
         Site = meta$site,
         Plot = meta$plot_num,
         Health_status = meta$tree_health)

# saveRDS(scores, file = "abiotic-data/abiotic-pca-scores.rds")

# Generate data for centroids
centroid_dead <- scores %>% 
  select(PC1, PC2, Health_status) %>% 
  filter(Health_status == "dead") %>% 
  select(-Health_status)
centroid_dead <- colMeans(centroid_dead)

centroid_defol <- scores %>% 
  select(PC1, PC2, Health_status) %>% 
  filter(Health_status == "defoliated") %>% 
  select(-Health_status)
centroid_defol <- colMeans(centroid_defol)

centroid_heal <- scores %>% 
  select(PC1, PC2, Health_status) %>% 
  filter(Health_status == "healthy") %>% 
  select(-Health_status)
centroid_heal <- colMeans(centroid_heal)

centroid_th <- data.frame(centroid_dead, centroid_defol, centroid_heal) %>% 
  t() %>% 
  data.frame(Health_status = c("dead", "defoliated", "healthy"))

# saveRDS(centroid_th, "abiotic-data/abiotic-tree-health-centroid-df.rds")

centroid_cas <- scores %>% 
  select(PC1, PC2, Site) %>% 
  filter(Site == "Cascades") %>% 
  select(-Site)
centroid_cas <- colMeans(centroid_cas)

centroid_piha <- scores %>% 
  select(PC1, PC2, Site) %>% 
  filter(Site == "Piha") %>% 
  select(-Site)
centroid_piha <- colMeans(centroid_piha)

centroid_huia <- scores %>% 
  select(PC1, PC2, Site) %>% 
  filter(Site == "Huia") %>% 
  select(-Site)
centroid_huia <- colMeans(centroid_huia)

centroid_site <- data.frame(centroid_cas, centroid_piha, centroid_huia) %>% 
  t() %>% 
  data.frame(Site = c("Cascades", "Piha", "Huia"))

# saveRDS(centroid_site, "abiotic-data/abiotic-site-centroid-df.rds")

#***********
# Plot----
#***********
# Quick plot to view
# ggplot(data = scores, aes(x = PC1, y = PC2)) +
#   geom_point(aes(color = Site, shape = Health_status), size = 3) +
#   geom_segment(data = loadings, aes(x = 0, y = 0, xend = x, yend = y),
#                arrow = arrow(length = unit(0.25, "cm")), color = "gray30") +
#   geom_text(data = loadings, aes(x = x, y = y, label = Variable), 
#             hjust = 0.5, vjust = 1.2, size = 4, color = "black") +
#   labs(x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
#        y = paste0("PC2 (", round(var_exp[2], 1), "%)")) +
#   theme_classic() +
#   coord_fixed(ratio = 1) +
#   theme(text = element_text(size = 12))

# Final plot with metadata
ggplot(data = scores, aes(x = PC1, y = PC2)) + # Increase values on each axis to fill the plot more
  geom_point(
    aes(
      shape = Site, # Set shape to relate to site
      fill = Health_status, # Set fill (colour inside shape) to relate to tree health
      colour = ifelse(Plot == "Plot 1", "black", "NA") # Set colour (outside/border colour) to be "black" for Plot 1 and empty if Plot 2
    ),
    size = 4, # Increase size of points
    stroke = 1, # Increase border thickness of points
    alpha = 0.6 # Decrease opacity of point
  ) +
  scale_shape_manual( # Manually set the shapes for the site variable
    name = "Site", 
    values = c(21, 22, 24), # Shapes must be between options 21-25 to allow both fill and colour options
    guide = "none" # Remove legend
  ) + 
  scale_colour_manual(
    name = "Plot", # Name on legend
    na.value = NA, # Set colour for NA values to appear blank
    values = c("black", "NA"),
    guide = "none" # Remove legend
  ) +
  scale_fill_manual( # Manually set the colours of fill variable (tree health)
    name = "Tree Health Status", # Change name that appears on legend
    labels = c("Dead", "Defoliated", "Healthy"), # Set labels to appear on legend
    values = c("#FB8072", "#BEBADA", "#B3DE69")
  ) +
  new_scale_colour() + # Reset scales to use (easier to add new layers)
  geom_point(data = centroid_th, aes(x = PC1, y = PC2, colour = Health_status),
             size = 8,
             stroke = 2,
             alpha = 0.9, 
             shape = 4, 
             show.legend = F # Remove legend as this layer shows just the centroids of the same variables above
  ) +
  scale_colour_manual(
    name = "Tree Health Status",
    labels = c("Dead", "Defoliated", "Healthy"),
    values = c("#F94C39", "#6B60A9", "#83B828")
  ) +
  geom_point(data = centroid_site, aes(x = PC1, PC2, shape = Site),
             size = 8,
             stroke = 2) +
  guides(fill = guide_legend(override.aes = list( # override aesthetics of the legend
    shape = 21,
    alpha = 1,
    colour = "NA",
    size = 5
  )),
  shape = guide_legend(override.aes = list(size = 5, stroke = 1))) +
  guides(fill = guide_legend(override.aes = list( # override aesthetics of the legend
    shape = 21,
    alpha = 1,
    colour = "NA",
    size = 5
  )),
  shape = guide_legend(override.aes = list(size = 5, stroke = 1))) +
  labs(x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
       y = paste0("PC2 (", round(var_exp[2], 1), "%)")) +
  theme_classic(base_size = 20) +
  coord_fixed() + ## need aspect ratio of 1! (to add arrows)
  geom_segment( # Add arrows
    data = loadings,
    aes(
      x = 0,
      xend = x,
      y = 0,
      yend = y
    ),
    arrow = arrow(length = unit(0.25, "cm")),
    colour = "grey50"
  ) +
  geom_text_repel( # Add text at tips of arrows
    data = loadings,
    aes(x = x, y = y, label = envname),
    size = 6,
    position = position_nudge_center(0.01, 0.025, 0.02, 0.025, direction = "radial") # Nudge text so it doesn't overlap with arrow tips
  ) +
  theme(legend.spacing.x = unit(1.0, "cm"), legend.key.spacing.y = unit(0.3, "cm")) + # Increase spacing of legend text
  theme(axis.text = element_text(colour = "black"))

ggsave(height=10, width=20, dpi=600, plot=last_plot(), filename = "abiotic-PCA-centroids.pdf")

## End of script ##