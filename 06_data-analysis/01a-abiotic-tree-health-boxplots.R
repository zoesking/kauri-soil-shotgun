## ---------------------------
##
## Script name: abiotic-tree-health-boxplots.R
##
## Purpose of script: Generate boxplots and run Wilcoxon tests to 
## assess the differences between soil physicochemical properties
## between healthy vs unhealthy trees.
##
## Author: Zoe King
##
## Date Created: 2025-08-28
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
library(ggpubr)

#***********
# Data----
#***********
# Load in metadata
meta <- read_csv("Master-data-file-shotgun.csv") %>% 
  filter(type == "Tree") %>% 
  filter(tree_code != "UC-D2") %>% 
  mutate(health_status = case_when(
    tree_health == "healthy" ~ "Healthy",
    tree_health == "defoliated" ~ "Unhealthy",
    tree_health == "dead" ~ "Unhealthy"
  ))

# Extract health status info
health_status <- meta$health_status

# Select physicochemical data and health status columns
physico <- meta %>% 
  dplyr::select(tree_code, health_status, pH, electrical_conductivity, total_hydrogen, carbon_to_nitrogen, bulk_density,  water_holding_capacity, total_carbon, total_nitrogen, moisture_factor) %>% # remove total_carbon, total_nitrogen, moisture_factor (significant results)
  column_to_rownames(var = "tree_code")

#*******************************
# Prepare data for plotting----
#*******************************
# Pretty Y-axis labels with units
var_labels <- c(
  total_hydrogen = " Mean total hydrogen (%)",
  total_carbon = "Mean total carbon (%)",
  total_nitrogen  = "Mean total nitrogen (%)",
  carbon_to_nitrogen = "Mean C:N ratio",
  pH  = "Mean pH",
  electrical_conductivity = " Mean Electrical conductivity \n (\u00B5S/cm)", # μS/cm
  bulk_density = "Bulk density (g/cm\u00B3)",  # g/cm³
  moisture_factor = "Soil moisture (%)",
  water_holding_capacity = "Water holding capacity (%)"
)

# Specific order for facets
facet_order <- names(var_labels)

# Reshape to long format
soil_long <- physico %>% 
  pivot_longer(-c(health_status),
               names_to = "variable", values_to = "value") %>% 
  mutate(health_status = factor(health_status, levels = c("Healthy", "Unhealthy")),
         variable      = factor(variable, levels = facet_order),
         var_label = forcats::fct_relabel(variable, ~ var_labels[.x]))

# saveRDS(soil_long, file = "gene-catalog-data/final-R-code/rds-objects/abiotic-univariate-data.rds")

#***********
# Plot----
#***********
# Significant comparisons
soil_long %>% 
  filter(variable %in% c("total_carbon", "total_nitrogen", "moisture_factor")) %>% 
  ggplot(aes(x = health_status, y = value)) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 2, aes(color = health_status)) +
  geom_boxplot(outlier.shape = NA, fill = NA) +
  scale_colour_manual(
    name = "Tree Health Status", # Change name that appears on legend
    labels = c("Healthy", "Unhealthy"), # Set labels to appear on legend
    values = c("#B3DE69", "#FB8072")
  ) +
  facet_wrap(~ var_label, scales = "free_y", strip.position = "left") +
  labs(x = "Tree health", y = NULL) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.signif", # stars only ("*", "**", etc.)
    #hide.ns = TRUE,
    comparisons = list(c("Healthy", "Unhealthy")),  # define the 2-level comparison
    bracket.size = 0.5, tip.length = 0.02) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
  theme_bw(base_size = 20) +
  theme(panel.grid = element_blank(),
        strip.placement    = "outside",
        strip.background = element_blank(),
        strip.text.y.left  = element_text(angle = 90, hjust = 0.5),  # make strips look like y-axis labels
        legend.position = "none",
        axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"))

ggsave(height=8, width=16, dpi=600, plot=last_plot(), filename = "abiotic-tree-health-boxplots-sig-diff.pdf")

# Non-significant comparisons
soil_long %>% 
  filter(!variable %in% c("total_carbon", "total_nitrogen", "moisture_factor")) %>% 
  ggplot(aes(x = health_status, y = value)) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 2, aes(color = health_status)) +
  geom_boxplot(outlier.shape = NA, fill = NA) +
  scale_colour_manual(
    name = "Tree Health Status", # Change name that appears on legend
    labels = c("Healthy", "Unhealthy"), # Set labels to appear on legend
    values = c("#B3DE69", "#FB8072")
  ) +
  facet_wrap(~ var_label, scales = "free_y", strip.position = "left") +
  labs(x = "Tree health", y = NULL) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.signif", # stars only ("*", "**", etc.)
    comparisons = list(c("Healthy", "Unhealthy")),  # define the 2-level comparison
    bracket.size = 0.5, tip.length = 0.02) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
  theme_bw(base_size = 20) +
  theme(panel.grid = element_blank(),
        strip.placement    = "outside",
        strip.background = element_blank(),
        strip.text.y.left  = element_text(angle = 90, hjust = 0.5),  # make strips look like y-axis labels
        legend.position = "none",
        axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"))

# Save plot
ggsave(height=8, width=16, dpi=600, plot=last_plot(), filename = "abiotic-tree-health-boxplots-not-sig.pdf")

## End of Script ##