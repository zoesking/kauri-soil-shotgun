## ---------------------------
##
## Script name: gene-richness-boxplots.R
##
## Purpose of script: Plot richness boxplots 
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
egg <- readRDS("gene-catalog-data/R-scripts/eggnog-richness-df.rds")
ncyc <- readRDS("gene-catalog-data/R-scripts/ncyc-richness-df.rds")
cazy <- readRDS("gene-catalog-data/R-scripts/cazy-richness-df.rds")

#****************
# eggNOG plot----
#****************
ggplot(egg, aes(x = tree_health, y = richness)) +
  geom_jitter(aes(color = tree_health), size = 4, alpha = 0.6) +
  geom_boxplot(outlier.color = NA, fill = NA, linewidth = 1) +
  scale_color_manual( # Manually set the colours of fill variable (tree health)
    name = "Tree Health Status", # Change name that appears on legend
    labels = c("Healthy", "Defoliated", "Dead"), # Set labels to appear on legend
    values = c("#B3DE69", "#BEBADA", "#FB8072")
  ) +
  labs(x = "Tree Health", y = "Observed Richness") +
  theme_classic(base_size = 20) +
  theme(legend.position = "none") +
  theme(axis.text = element_text(colour = "black"), axis.ticks = element_line(colour = "black")) +
  labs2(tag = "B.") +
  theme(plot.title.position = 'plot', plot.title = ggtext::element_markdown(hjust = 0, vjust = 1, margin = margin(t = 1, r = 0, b = 1, l = 0, unit = 'mm'), size = 20))

ggsave(height=7, width=7, dpi=600, plot=last_plot(), filename = "gene-catalog-data/R-scripts/eggnog-richness-boxplot.pdf")

#****************
# NCyc plot----
#****************
ggplot(ncyc, aes(x = tree_health, y = richness)) +
  geom_jitter(aes(color = tree_health), size = 4, alpha = 0.6) +
  geom_boxplot(outlier.color = NA, fill = NA, linewidth = 1) +
  scale_color_manual( # Manually set the colours of fill variable (tree health)
    name = "Tree Health Status", # Change name that appears on legend
    labels = c("Healthy", "Defoliated", "Dead"), # Set labels to appear on legend
    values = c("#B3DE69", "#BEBADA", "#FB8072")
  ) +
  labs(x = "Tree Health", y = "Observed Richness") +
  theme_classic(base_size = 20) +
  theme(legend.position = "none") +
  theme(axis.text = element_text(colour = "black"), axis.ticks = element_line(colour = "black")) +
  labs2(tag = "B.") +
  theme(plot.title.position = 'plot', plot.title = ggtext::element_markdown(hjust = 0, vjust = 1, margin = margin(t = 1, r = 0, b = 1, l = 0, unit = 'mm'), size = 20))

ggsave(height=7, width=7, dpi=600, plot=last_plot(), filename = "gene-catalog-data/R-scripts/ncyc-richness-boxplot.pdf")

#****************
# CAZy plot----
#****************
ggplot(cazy, aes(x = tree_health, y = richness)) +
  geom_jitter(aes(color = tree_health), size = 4, alpha = 0.6) +
  geom_boxplot(outlier.color = NA, fill = NA, linewidth = 1) +
  scale_color_manual( # Manually set the colours of fill variable (tree health)
    name = "Tree Health Status", # Change name that appears on legend
    labels = c("Healthy", "Defoliated", "Dead"), # Set labels to appear on legend
    values = c("#B3DE69", "#BEBADA", "#FB8072")
  ) +
  labs(x = "Tree Health", y = "Observed Richness") +
  theme_classic(base_size = 20) +
  theme(legend.position = "none") +
  theme(axis.text = element_text(colour = "black"), axis.ticks = element_line(colour = "black")) +
  labs2(tag = "B.") +
  theme(plot.title.position = 'plot', plot.title = ggtext::element_markdown(hjust = 0, vjust = 1, margin = margin(t = 1, r = 0, b = 1, l = 0, unit = 'mm'), size = 20))

ggsave(height=7, width=7, dpi=600, plot=last_plot(), filename = "gene-catalog-data/R-scripts/cazy-richness-boxplot.pdf")


## End of Script ##