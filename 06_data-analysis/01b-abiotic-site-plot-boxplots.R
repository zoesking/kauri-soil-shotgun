## ---------------------------
##
## Script name: abiotic_site_plot_boxplots.R
##
## Purpose of script: Generate boxplots and run Kruskal-Wallis tests 
## to assess the differences between soil physicochemical properties
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
library(rstatix)
library(ggpubr)

#***********
# Data----
#***********
# Load in metadata
meta <- read_csv("Master-data-file-shotgun.csv") %>% 
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

# Select physicochemical data and health status columns
physico <- meta %>% 
  select(tree_code, site, plot_num, pH, electrical_conductivity, total_carbon, total_hydrogen, total_nitrogen, carbon_to_nitrogen, bulk_density, moisture_factor, water_holding_capacity) %>% 
  column_to_rownames(var = "tree_code")

#**************************
# Prepare Y-axis labels----
#**************************
# Pretty Y-axis labels with units
var_labels <- c(
  total_hydrogen = "Mean total hydrogen (%)",
  total_carbon = "Mean total carbon (%)",
  total_nitrogen  = "Mean total nitrogen (%)",
  carbon_to_nitrogen = "Mean C:N ratio",
  pH  = "Mean pH",
  electrical_conductivity = "Mean Electrical conductivity \n (\u00B5S/cm)", # μS/cm
  bulk_density = "Bulk density (g/cm\u00B3)",  # g/cm³
  moisture_factor = "Soil moisture (%)",
  water_holding_capacity = "Water holding capacity (%)"
)

# Specific order for facets
facet_order <- names(var_labels)

#************************************
# Prepare site data for plotting-----
#************************************
# Convert to long format
soil_long <- physico %>%
  mutate(Site = factor(site, levels = c("Cascades","Piha","Huia"))) %>%
  pivot_longer(cols = all_of(facet_order), names_to = "variable", values_to = "value") %>%
  mutate(
    variable  = factor(variable, levels = facet_order),
    var_label = factor(var_labels[as.character(variable)], levels = unname(var_labels))
  ) %>%
  drop_na(value)

### Kruskal-Wallis and Dunn test====
# Kruskal–Wallis across the 3 sites (per variable)
kw_tbl <- soil_long %>%
  group_by(var_label) %>%
  kruskal_test(value ~ Site) %>%
  add_significance()

# Dunn post-hoc pairwise (BH)
site_pairs <- list(c("Cascades","Piha"), c("Cascades","Huia"), c("Piha","Huia"))

dunn_tbl <- soil_long %>%
  group_by(var_label) %>%
  dunn_test(value ~ Site, p.adjust.method = "BH") %>%
  ungroup() %>%
  add_significance("p.adj")

### Prepare bracket positions====
# Base y for each facet (top of data range + small margin)
y_base <- soil_long %>%
  group_by(var_label) %>%
  summarise(
    ymax   = max(value, na.rm = TRUE),
    yrange = diff(range(value, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(y.pos.base = ymax + 0.05 * ifelse(yrange > 0, yrange, 1))

# Order pairs so they stack predictably (top to bottom)
pair_order <- tibble(
  group1 = c("Cascades","Cascades","Piha"),
  group2 = c("Piha","Huia","Huia"),
  pair_rank = 1:3
)

### Create plotting dataframe====
dunn_anno <- dunn_tbl %>%
  left_join(pair_order, by = c("group1","group2")) %>%
  group_by(var_label) %>%
  arrange(pair_rank, .by_group = TRUE) %>%
  ungroup() %>%
  left_join(y_base, by = "var_label") %>%
  mutate(
    y.position = y.pos.base + (pmax(yrange, 1) * 0.08) * (pair_rank - 1),
    label = p.adj.signif   # <-- use stars instead of numbers
  ) %>%
  select(var_label, group1, group2, y.position, label)

# Create unique colour palette
site_cols <- c(
  "Cascades" = "#1b9e77",
  "Piha"     = "#d95f02",
  "Huia"     = "#7570b3"
)

#**************************
# Plot site boxplots----
#**************************
ggplot(soil_long, aes(x = Site, y = value)) +
  geom_boxplot(width = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.6, aes(colour = Site)) +
  scale_colour_manual(values = site_cols, name = "Site") +
  facet_wrap(~ var_label, scales = "free_y", strip.position = "left", ncol = 3) +
  labs(x = "Site", y = NULL) +
  theme_bw(base_size = 20) +
  theme(
    panel.grid = element_blank(),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 90, hjust = 0.5),
    legend.position = "none",
    axis.text = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black")
  ) +
  ggpubr::stat_pvalue_manual(
    dunn_anno,
    label = NULL,
    xmin  = "group1",
    xmax  = "group2",
    y.position = "y.position",
    tip.length = 0.01,
    bracket.size = 0.3
  )

# Save plot
ggsave(height=10, width=16, dpi=600, plot=last_plot(), filename = "abiotic-site-boxplots.pdf")

#**********************************************
# Prepare plot/location data for plotting-----
#**********************************************
# Add combined location information
soil_long_sp <- soil_long %>%
  mutate(Site_Plot = paste(site, plot_num, sep = "_"))
# Set levels
siteplot_levels <- c("Cascades_Plot 1","Cascades_Plot 2",
                     "Piha_Plot 1","Piha_Plot 2",
                     "Huia_Plot 1","Huia_Plot 2")
# Convert to long format
soil_long_sp <- soil_long_sp %>%
  mutate(
    Site       = factor(Site, levels = c("Cascades","Piha","Huia")),
    Plot       = factor(plot_num, levels = c("Plot 1","Plot 2")),
    Site_Plot  = factor(Site_Plot, levels = siteplot_levels),
    # ensure variable and var_label exist
    variable   = factor(variable, levels = facet_order),
    var_label  = factor(var_labels[as.character(variable)], levels = unname(var_labels))
  ) %>%
  select(c(Site_Plot, variable, value, var_label)) %>% 
  drop_na(value)

### Kruskal-Wallis and Dunn test====
# Kruskal–Wallis across the 6 plots, per variable
kw_tbl_sp <- soil_long_sp %>%
  group_by(var_label) %>%
  kruskal_test(value ~ Site_Plot) %>%
  add_significance()

# Dunn post-hoc pairwise across the 6 plots (BH)
dunn_tbl_sp <- soil_long_sp %>%
  group_by(var_label) %>%
  dunn_test(value ~ Site_Plot, p.adjust.method = "BH") %>%
  ungroup() %>%
  add_significance("p.adj")

# Keep only significant comparisons to reduce clutter
dunn_tbl_sp <- dunn_tbl_sp %>% filter(p.adj < 0.05)

### Prepare bracket positions====
# Base y for each facet (top of data range + small margin)
y_base_sp <- soil_long_sp %>%
  group_by(var_label) %>%
  summarise(
    ymax   = max(value, na.rm = TRUE),
    yrange = diff(range(value, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(y.pos.base = ymax + 0.05 * ifelse(yrange > 0, yrange, 1))

### Create plotting dataframe====
# Order comparisons (bottom-up stacking by increasing p-value within each facet)
dunn_anno <- dunn_tbl_sp %>%
  group_by(var_label) %>%
  arrange(p.adj, .by_group = TRUE) %>%
  mutate(row_id = row_number()) %>%
  ungroup() %>%
  left_join(y_base_sp, by = "var_label") %>%
  mutate(
    y.position = y.pos.base + (pmax(yrange, 1) * 0.08) * (row_id - 1),
    label = p.adj.signif  # use stars only: "ns", "*", "**", "***", "****"
  ) %>%
  # Keep only the columns ggpubr needs
  select(var_label, group1, group2, y.position, label)

# Unique colours for plotting
sp_plot_cols <- c(
  "Cascades_Plot 1" = "#1b9e77",
  "Cascades_Plot 2" = "#d95f02",
  "Piha_Plot 1"     = "#7570b3",
  "Piha_Plot 2"     = "#e7298a",
  "Huia_Plot 1"     = "#66a61e",
  "Huia_Plot 2"     = "#e6ab02"
)

#************************************
# Plot plot/location boxplots----
#************************************
ggplot(soil_long_sp, aes(x = Site_Plot, y = value)) +
  geom_boxplot(width = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.6, show.legend = FALSE, aes(color = Site_Plot)) +
  scale_color_manual(values = sp_plot_cols, name = "Plot") +
  facet_wrap(~ var_label, scales = "free_y", strip.position = "left", ncol = 3) +
  labs(x = "Plot", y = NULL) +
  theme_bw(base_size = 20) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 90, hjust = 0.5),
    axis.text = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black")
  ) +
  ggpubr::stat_pvalue_manual(
    dunn_anno,
    label = "label",
    xmin  = "group1", # must exactly match factor levels of Site_Plot
    xmax  = "group2",
    y.position = "y.position",
    tip.length = 0.01,
    bracket.size = 0.3
  )

# Save plot
ggsave(height=13, width=16, dpi=600, plot=last_plot(), filename = "abiotic-plot-boxplots.pdf")

## End of Script ##