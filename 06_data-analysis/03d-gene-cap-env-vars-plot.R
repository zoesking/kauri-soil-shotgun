## ---------------------------
##
## Script name: gene-cap-env-vars-plot.R
##
## Purpose of script: Plot constrained ordination
## with env vars overlaid as vectors (from envfit)
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
library(ggrepel)
library(ggnewscale) #ggnewscale: spend 400% more time tweaking your ggplot!
library(patchwork)
library(ggtext)
library(ggpp)

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
ord_df_egg <- readRDS("rds-objects/eggnog-bray-ord-meta-df.rds")
fit_scores_egg <- readRDS("rds-objects/eggnog-bray-cap-envfit-scores.rds")
values_egg <- readRDS("rds-objects/eggnog-bray-cap-scores.rds")
cap_egg <- readRDS("rds-objects/eggnog-bray-cap-ordination.rds")

ord_df_ncyc <- readRDS("rds-objects/ncyc-bray-ord-meta-df.rds")
fit_scores_ncyc <- readRDS("rds-objects/ncyc-bray-cap-envfit-scores.rds")
values_ncyc <- readRDS("rds-objects/ncyc-bray-cap-scores.rds")
cap_ncyc <- readRDS("rds-objects/ncyc-bray-cap-ordination.rds")

ord_df_cazy <- readRDS("rds-objects/cazy-bray-ord-meta-df.rds")
fit_scores_cazy <- readRDS("rds-objects/cazy-bray-cap-envfit-scores.rds")
values_cazy <- readRDS("rds-objects/cazy-bray-cap-scores.rds")
cap_cazy <- readRDS("rds-objects/cazy-bray-cap-ordination.rds")

#****************
# eggNOG plot----
#****************
ggplot(ord_df_egg, aes(x = CAP1, y = CAP2)) +
  geom_point(
    aes(
      shape = site, # Set shape to relate to site
      fill = tree_health, # Set fill (colour inside shape) to relate to tree health
      colour = ifelse(plot_num == "Plot 1", "black", "NA") # Set colour (outside/border colour) to be "black" for Plot 1 and empty if Plot 2
    ),
    size = 4, # Increase size of points
    stroke = 1, # Increase border thickness of points
    alpha = 0.6 # Decrease opacity of point
  ) +
  scale_shape_manual( # Manually set the shapes for the site variable
    name = "Site", 
    values = c(21, 22, 24), # Shapes must be between options 21-25 to allow both fill and colour options
    # guide = "none" # Remove legend
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
    values = c("#FB8072","#BEBADA", "#B3DE69")
  )  +
  labs(x = paste0("CAP1 (", round(cap_egg$CCA$eig[1] / sum(cap_egg$CCA$eig) * 100, 1), "%)"),
       y = paste0("CAP2 (", round(cap_egg$CCA$eig[2] / sum(cap_egg$CCA$eig) * 100, 1), "%)")) +
  theme_classic(base_size = 20) +
  theme(axis.text = element_text(colour = "black"), axis.ticks = element_line(color = "black")) +
  guides(fill = guide_legend(override.aes = list( # override aesthetics of the legend
    shape = 21,
    alpha = 1,
    colour = "NA",
    size = 5
  )),
  shape = guide_legend(override.aes = list(size = 5, stroke = 1))) +
  coord_fixed() + ## need aspect ratio of 1! (to add arrows)
  geom_segment( # Add arrows
    data = fit_scores_egg,
    aes(
      x = 0,
      xend = CAP1,
      y = 0,
      yend = CAP2
    ),
    arrow = arrow(length = unit(0.25, "cm")),
    colour = "grey50"
  ) +
  geom_text_repel( # Add text at tips of arrows
    data = fit_scores_egg,
    aes(x = CAP1, y = CAP2, label = envname),
    size = 6,
    position = ggpp::position_nudge_center(0.01, 0.025, 0.02, 0.025, direction = "radial") # Nudge text so it doesn't overlap with arrow tips
  ) +
  labs2(tag = "B.") +
  theme(plot.title.position = 'plot', plot.title = ggtext::element_markdown(hjust = 0, vjust = 1, margin = margin(t = 0, r = 0, b = 1, l = 0, unit = 'mm')))

ggsave(height=8, width=12, dpi=600, plot=last_plot(), filename = "../eggnog-bray-gene-cap-env.pdf")

#****************
# NCyc plot----
#****************
ggplot(ord_df_ncyc, aes(x = CAP1, y = CAP2)) +
  geom_point(
    aes(
      shape = site, # Set shape to relate to site
      fill = tree_health, # Set fill (colour inside shape) to relate to tree health
      colour = ifelse(plot_num == "Plot 1", "black", "NA") # Set colour (outside/border colour) to be "black" for Plot 1 and empty if Plot 2
    ),
    size = 4, # Increase size of points
    stroke = 1, # Increase border thickness of points
    alpha = 0.6 # Decrease opacity of point
  ) +
  scale_shape_manual( # Manually set the shapes for the site variable
    name = "Site", 
    values = c(21, 22, 24), # Shapes must be between options 21-25 to allow both fill and colour options
    # guide = "none" # Remove legend
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
    values = c("#FB8072","#BEBADA", "#B3DE69")
  )  +
  labs(x = paste0("CAP1 (", round(cap_ncyc$CCA$eig[1] / sum(cap_ncyc$CCA$eig) * 100, 1), "%)"),
       y = paste0("CAP2 (", round(cap_ncyc$CCA$eig[2] / sum(cap_ncyc$CCA$eig) * 100, 1), "%)")) +
  theme_classic(base_size = 20) +
  theme(axis.text = element_text(colour = "black"), axis.ticks = element_line(color = "black")) +
  guides(fill = guide_legend(override.aes = list( # override aesthetics of the legend
    shape = 21,
    alpha = 1,
    colour = "NA",
    size = 5
  )),
  shape = guide_legend(override.aes = list(size = 5, stroke = 1))) +
  coord_fixed() + ## need aspect ratio of 1! (to add arrows)
  geom_segment( # Add arrows
    data = fit_scores_ncyc,
    aes(
      x = 0,
      xend = CAP1,
      y = 0,
      yend = CAP2
    ),
    arrow = arrow(length = unit(0.25, "cm")),
    colour = "grey50"
  ) +
  geom_text_repel( # Add text at tips of arrows
    data = fit_scores_ncyc,
    aes(x = CAP1, y = CAP2, label = envname),
    size = 6,
    position = ggpp::position_nudge_center(0.01, 0.025, 0.02, 0.025, direction = "radial") # Nudge text so it doesn't overlap with arrow tips
  ) +
  labs2(tag = "B.") +
  theme(plot.title.position = 'plot', plot.title = ggtext::element_markdown(hjust = 0, vjust = 1, margin = margin(t = 0, r = 0, b = 1, l = 0, unit = 'mm')))

ggsave(height=8, width=12, dpi=600, plot=last_plot(), filename = "../ncyc-bray-gene-cap-env.pdf")

#****************
# CAZy plot----
#****************
ggplot(ord_df_cazy, aes(x = CAP1, y = CAP2)) +
  geom_point(
    aes(
      shape = site, # Set shape to relate to site
      fill = tree_health, # Set fill (colour inside shape) to relate to tree health
      colour = ifelse(plot_num == "Plot 1", "black", "NA") # Set colour (outside/border colour) to be "black" for Plot 1 and empty if Plot 2
    ),
    size = 4, # Increase size of points
    stroke = 1, # Increase border thickness of points
    alpha = 0.6 # Decrease opacity of point
  ) +
  scale_shape_manual( # Manually set the shapes for the site variable
    name = "Site", 
    values = c(21, 22, 24), # Shapes must be between options 21-25 to allow both fill and colour options
    # guide = "none" # Remove legend
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
    values = c("#FB8072","#BEBADA", "#B3DE69")
  )  +
  labs(x = paste0("CAP1 (", round(cap_cazy$CCA$eig[1] / sum(cap_cazy$CCA$eig) * 100, 1), "%)"),
       y = paste0("CAP2 (", round(cap_cazy$CCA$eig[2] / sum(cap_cazy$CCA$eig) * 100, 1), "%)")) +
  theme_classic(base_size = 20) +
  theme(axis.text = element_text(colour = "black"), axis.ticks = element_line(color = "black")) +
  guides(fill = guide_legend(override.aes = list( # override aesthetics of the legend
    shape = 21,
    alpha = 1,
    colour = "NA",
    size = 5
  )),
  shape = guide_legend(override.aes = list(size = 5, stroke = 1))) +
  coord_fixed() + ## need aspect ratio of 1! (to add arrows)
  geom_segment( # Add arrows
    data = fit_scores_cazy,
    aes(
      x = 0,
      xend = CAP1,
      y = 0,
      yend = CAP2
    ),
    arrow = arrow(length = unit(0.25, "cm")),
    colour = "grey50"
  ) +
  geom_text_repel( # Add text at tips of arrows
    data = fit_scores_cazy,
    aes(x = CAP1, y = CAP2, label = envname),
    size = 6,
    position = ggpp::position_nudge_center(0.01, 0.025, 0.02, 0.025, direction = "radial") # Nudge text so it doesn't overlap with arrow tips
  ) +
  labs2(tag = "B.") +
  theme(plot.title.position = 'plot', plot.title = ggtext::element_markdown(hjust = 0, vjust = 1, margin = margin(t = 0, r = 0, b = 1, l = 0, unit = 'mm')))

ggsave(height=8, width=12, dpi=600, plot=last_plot(), filename = "../cazy-bray-gene-cap-env.pdf")

## End of script ##