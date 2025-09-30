## ---------------------------
##
## Script name: gene-presence-venn.R
##
## Purpose of script:
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
## Final figure adjustments made in Illustrator
##
## ---------------------------

#****************
# Packages----
#****************
library(tidyverse)
library(ggvenn)

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
egg <- readRDS("gene-catalog-data/gene-stats/eggnog-gene-sets-venn.rds")
ncyc <- readRDS("gene-catalog-data/gene-stats/ncyc-gene-sets-venn.rds")
cazy <- readRDS("gene-catalog-data/gene-stats/cazy-gene-sets-venn.rds")

#****************
# eggNOG Venn----
#****************
# Rename
names(egg) <- c("Healthy\n(n = 21)", "Defoliated\n(n = 69)", "Dead\n(n = 5)")

# Plot
ggvenn(
  egg, 
  show_percentage = TRUE,
  digits = 3,
  stroke_size = 0.5,
  set_name_size = 5,
  text_size = 3,
  fill_color = c("#B3DE69", "#BEBADA", "#FB8072"),
  fill_alpha = 0.4
) +
  scale_x_continuous(expand = expansion(mult = .15)) +
  scale_y_continuous(expand = expansion(mult = .15)) +
  labs2(tag = "A.") +
  theme(plot.title.position = 'plot', plot.title = ggtext::element_markdown(hjust = 0, vjust = 1, margin = margin(t = 1, r = 1, b = 1, l = 0, unit = 'mm'), size = 20))

ggsave(height=6, width=6, dpi=600, plot=last_plot(), filename = "gene-catalog-data/R-scripts/eggnog-gene-pres-venn.pdf")

#****************
# NCyc Venn----
#****************
# Rename
names(ncyc) <- c("Healthy\n(n = 21)", "Defoliated\n(n = 69)", "Dead\n(n = 5)")

# Plot
ggvenn(
  ncyc, 
  show_percentage = TRUE,
  digits = 3,
  stroke_size = 0.5,
  set_name_size = 5,
  text_size = 3,
  fill_color = c("#B3DE69", "#BEBADA", "#FB8072"),
  fill_alpha = 0.4
) +
  scale_x_continuous(expand = expansion(mult = .15)) +
  scale_y_continuous(expand = expansion(mult = .15)) +
  labs2(tag = "A.") +
  theme(plot.title.position = 'plot', plot.title = ggtext::element_markdown(hjust = 0, vjust = 1, margin = margin(t = 1, r = 1, b = 1, l = 0, unit = 'mm'), size = 20))

ggsave(height=6, width=6, dpi=600, plot=last_plot(), filename = "gene-catalog-data/R-scripts/ncyc-gene-pres-venn.pdf")

#****************
# CAZy Venn----
#****************
# Rename
names(cazy) <- c("Healthy\n(n = 21)", "Defoliated\n(n = 69)", "Dead\n(n = 5)")

# Plot
ggvenn(
  cazy, 
  show_percentage = TRUE,
  digits = 3,
  stroke_size = 0.5,
  set_name_size = 5,
  text_size = 3,
  fill_color = c("#B3DE69", "#BEBADA", "#FB8072"),
  fill_alpha = 0.4
) +
  scale_x_continuous(expand = expansion(mult = .15)) +
  scale_y_continuous(expand = expansion(mult = .15)) +
  labs2(tag = "A.") +
  theme(plot.title.position = 'plot', plot.title = ggtext::element_markdown(hjust = 0, vjust = 1, margin = margin(t = 1, r = 1, b = 1, l = 0, unit = 'mm'), size = 20))

ggsave(height=6, width=6, dpi=600, plot=last_plot(), filename = "gene-catalog-data/R-scripts/cazy-gene-pres-venn.pdf")

## End of Script ##