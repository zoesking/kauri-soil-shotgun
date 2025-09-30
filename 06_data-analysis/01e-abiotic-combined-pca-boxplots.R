## ---------------------------
##
## Script name: abiotic-combined-pca-boxplots.R
##
## Purpose of script: Combine plots of PCA (with centroids)
## and boxplots of significant univariate analysis.
##
## Author: Zoe King
##
## Date Created: 2025-09-15
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
library(patchwork)
library(ggpubr)

### Functions====
# Helper function to add labels to plots
# https://stackoverflow.com/questions/71346978/how-to-position-plot-tag-and-plot-title-side-by-side-on-the-same-line
labs2 <- function(title = NULL, tag = NULL, ...) {
  if (!is.null(tag)) {
    tag <- paste0(tag)
    title <- paste(tag, title)
  }
  labs(title = title, ...)
}

#***********
# Data----
#***********
# Objects for PCA
scores <- readRDS("abiotic-data/abiotic-pca-scores.rds")
var_exp <- readRDS("abiotic-data/abiotic-pca-var-exp.rds")
loadings <- readRDS("abiotic-data/abiotic-loadings.rds")
centroid_th <- readRDS("abiotic-data/abiotic-tree-health-centroid-df.rds")
centroid_site <- readRDS("abiotic-data/abiotic-site-centroid-df.rds")

# Objects for boxplots
soil_long <- readRDS("abiotic-data/abiotic-univariate-data.rds") %>% 
  filter(variable %in% c("total_carbon", "total_nitrogen", "moisture_factor"))

#****************
# PCA Plot----
#****************
# PCA
p1 <- ggplot(data = scores, aes(x = PC1, y = PC2)) + # Increase values on each axis to fill the plot more
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
  theme(legend.spacing.x = unit(1.0, "cm"), legend.key.spacing.y = unit(0.3, "cm"), legend.position = "bottom", legend.title.position = "top", legend.title = element_text(hjust = 0.5)) + # Increase spacing of legend text
  theme(axis.text = element_text(colour = "black")) +
  labs2(tag = "A.") +
  theme(plot.title.position = 'plot', plot.title = ggtext::element_markdown(hjust = 0, vjust = 1, margin = margin(t = 0, r = 0, b = 1, l = 0, unit = 'mm')))

# Quick look
p1

#****************
# Box plots----
#****************
p2 <- ggplot(soil_long, aes(x = health_status, y = value)) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 2, aes(color = health_status)) +
  geom_boxplot(outlier.shape = NA, fill = NA) +
  scale_colour_manual(
    name = "Tree Health Status", # Change name that appears on legend
    labels = c("Healthy", "Unhealthy"), # Set labels to appear on legend
    values = c("#B3DE69", "#FB8072")
  ) +
  facet_wrap(~ var_label, scales = "free_y", strip.position = "left") +
  labs(x = NULL, y = NULL) +
  #stat_compare_means(method = "wilcox.test", label = "paste(..p.format.., ..p.signif..)", hjust = 0.75, vjust = 0.4) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.signif",            # stars only ("*", "**", etc.)
    #hide.ns = TRUE,                # hide "ns" if not significant
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
        axis.ticks = element_line(colour = "black")) +
  labs2(tag = "B.") +
  theme(plot.title.position = 'plot', plot.title = ggtext::element_markdown(hjust = 0, vjust = 1, margin = margin(t = 0, r = 0, b = 1, l = 0, unit = 'mm')))

# Quick look
p2

#*********************
# Combine plots----
#*********************
# p1 + inset_element(p2, left = 0.05, bottom = 0.25, right = 0.28, top = 0.98, align_to = "full")

p1 /
  p2 + plot_layout(heights = c(2, 1))

ggsave(height=12, width=16, dpi=600, plot=last_plot(), filename = "abiotic-PCA-centroids-boxplots.pdf")

## End of Script ##