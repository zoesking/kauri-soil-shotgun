## ---------------------------
##
## Script name: gene-composition-pcoa.R
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
# Distance matrices
egg <- readRDS("rds-objects/eggnog-bc-distance-matrix.rds")
ncyc <- readRDS("rds-objects/ncyc-bc-distance-matrix.rds")
cazy <- readRDS("rds-objects/cazy-bc-distance-matrix.rds")
# Metadata
meta <- read_csv("../../Master-data-file-shotgun.csv") %>% 
  filter(type == "Tree") %>% 
  filter(tree_code != "UC-D2") %>% 
  dplyr::select(tree_code, canopy_score, tree_health, location, site) %>% 
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

#****************
# eggNOG data----
#****************
### Ordination====
pcoa_egg <- stats::cmdscale(egg, eig=T, add=T)

# Sample coordinates for the first two axes
coords_egg <- as.data.frame(pcoa_egg$points[, 1:2, drop = FALSE])
colnames(coords_egg) <- c("PCo1", "PCo2")
coords_egg$sample <- rownames(coords_egg)

# % variance explained for axis labels (use only positive eigenvalues)
eig_egg <- pcoa_egg$eig
pos_eig_egg <- eig_egg[eig_egg > 0]
var_explained_egg <- eig_egg[1:2] / sum(pos_eig_egg)

# Join metadata
coords_egg <- merge(coords_egg, tibble::rownames_to_column(as.data.frame(meta), "sample"),
                     by = "sample", all.x = TRUE)

### Plot====
ggplot(coords_egg, aes(x = PCo1, y = PCo2)) +
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
  ) +
  labs(x = sprintf("PCoA1 (%.1f%%)", 100 * var_explained_egg[1]), 
       y = sprintf("PCoA2 (%.1f%%)", 100 * var_explained_egg[2])) +
  theme_classic(base_size = 20) +
  theme(axis.text = element_text(colour = "black"), axis.ticks = element_line(colour = "black")) +
  coord_fixed() +
  guides(fill = guide_legend(override.aes = list( # override aesthetics of the legend
    shape = 21,
    alpha = 1,
    colour = "NA",
    size = 5
  )),
  shape = guide_legend(override.aes = list(size = 5, stroke = 1))) +
  labs2(tag = "C.") +
  theme(plot.title.position = 'plot', plot.title = ggtext::element_markdown(hjust = 0, vjust = 1, margin = margin(t = 0, r = 0, b = 1.5, l = 0, unit = 'mm')))

ggsave(height=6, width=10, dpi=600, plot=last_plot(), filename = "../egg-bray-gene-pcoa.pdf")

#****************
# NCyc data----
#****************
### Ordination====
pcoa_ncyc <- stats::cmdscale(ncyc, eig=T, add=T)

# Sample coordinates for the first two axes
coords_ncyc <- as.data.frame(pcoa_ncyc$points[, 1:2, drop = FALSE])
colnames(coords_ncyc) <- c("PCo1", "PCo2")
coords_ncyc$sample <- rownames(coords_ncyc)

# % variance explained for axis labels (use only positive eigenvalues)
eig_ncyc <- pcoa_ncyc$eig
pos_eig_ncyc <- eig_ncyc[eig_ncyc > 0]
var_explained_ncyc <- eig_ncyc[1:2] / sum(pos_eig_ncyc)

# Join metadata
coords_ncyc <- merge(coords_ncyc, tibble::rownames_to_column(as.data.frame(meta), "sample"),
                by = "sample", all.x = TRUE)

### Plot====
ggplot(coords_ncyc, aes(x = PCo1, y = PCo2)) +
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
  ) +
  labs(x = sprintf("PCoA1 (%.1f%%)", 100 * var_explained_ncyc[1]), 
       y = sprintf("PCoA2 (%.1f%%)", 100 * var_explained_ncyc[2])) +
  theme_classic(base_size = 20) +
  theme(axis.text = element_text(colour = "black"), axis.ticks = element_line(colour = "black")) +
  coord_fixed() +
  guides(fill = guide_legend(override.aes = list( # override aesthetics of the legend
    shape = 21,
    alpha = 1,
    colour = "NA",
    size = 5
  )),
  shape = guide_legend(override.aes = list(size = 5, stroke = 1))) +
  labs2(tag = "C.") +
  theme(plot.title.position = 'plot', plot.title = ggtext::element_markdown(hjust = 0, vjust = 1, margin = margin(t = 0, r = 0, b = 1.5, l = 0, unit = 'mm')))

ggsave(height=6, width=10, dpi=600, plot=last_plot(), filename = "../ncyc-bray-gene-pcoa.pdf")

#****************
# CAZy data----
#****************
### Ordination====
pcoa_cazy <- stats::cmdscale(cazy, eig=T, add=T)

# Sample coordinates for the first two axes
coords_cazy <- as.data.frame(pcoa_cazy$points[, 1:2, drop = FALSE])
colnames(coords_cazy) <- c("PCo1", "PCo2")
coords_cazy$sample <- rownames(coords_cazy)

# % variance explained for axis labels (use only positive eigenvalues)
eig_cazy <- pcoa_cazy$eig
pos_eig_cazy <- eig_cazy[eig_cazy > 0]
var_explained_cazy <- eig_cazy[1:2] / sum(pos_eig_cazy)

# Join metadata
coords_cazy <- merge(coords_cazy, tibble::rownames_to_column(as.data.frame(meta), "sample"),
                    by = "sample", all.x = TRUE)

### Plot====
ggplot(coords_cazy, aes(x = PCo1, y = PCo2)) +
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
    values = c(21, 22, 24) # Shapes must be between options 21-25 to allow both fill and colour options
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
  ) +
  labs(x = sprintf("PCoA1 (%.1f%%)", 100 * var_explained_cazy[1]), 
       y = sprintf("PCoA2 (%.1f%%)", 100 * var_explained_cazy[2])) +
  theme_classic(base_size = 20) +
  theme(axis.text = element_text(colour = "black"), axis.ticks = element_line(colour = "black")) +
  coord_fixed() +
  guides(fill = guide_legend(override.aes = list( # override aesthetics of the legend
    shape = 21,
    alpha = 1,
    colour = NA,
    size = 5
  )),
  shape = guide_legend(override.aes = list(size = 5, stroke = 1, colour = "black"))) +
  labs2(tag = "C.") +
  theme(plot.title.position = 'plot', plot.title = ggtext::element_markdown(hjust = 0, vjust = 1, margin = margin(t = 0, r = 0, b = 1.5, l = 0, unit = 'mm')))

ggsave(height=6, width=10, dpi=600, plot=last_plot(), filename = "../cazy-bray-gene-pcoa.pdf")

## End of Script ##