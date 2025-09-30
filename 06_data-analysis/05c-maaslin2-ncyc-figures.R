## ---------------------------
##
## Script name: maaslin2-ncyc-figures.R
##
## Purpose of script:
##
## Author: Zoe King
##
## Date Created: 2025-08-18
##
## Email: zoe.s.king@autuni.ac.nz OR zasking@gmail.com
##
## ---------------------------
##
## Notes:
## 'coef' is the effect size (positive or negative association)
## 'qval' is the FDR-adjusted p-value
##
## ---------------------------

#****************
# Packages----
#****************
library(tidyverse)
library(patchwork)

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
# Load in all maaslin2 results
maas <- read_tsv("../R-scripts/maaslin2_ncyc_output/all_results.tsv")
# Load in maaslin2 results (significant results only)
maas_passed <- read_tsv("../R-scripts/maaslin2_ncyc_output/significant_results.tsv") %>% 
  rename("query" = "feature")
# Annotated genes information
anno <- read_csv("../annotation-tables/ncyc-mapped-output.csv") %>% 
  select(query, gene)
# Higher function information
func <- read_tsv("../ncyc-gene-mapping.tsv")

#*********************
# Volcano plot----
#*********************
### Tidy data====
# Add a column to color points based on significance
maas_col <- maas %>%
  mutate(significant = ifelse(qval < 0.05, "Yes", "No"))

# Add significance and direction info
maas_col <- maas %>%
  mutate(
    direction = case_when(
      qval < 0.05 & coef > 0 ~ "Higher in Unhealthy",
      qval < 0.05 & coef < 0 ~ "Higher in Healthy",
      TRUE ~ "Not significant"
    )
  )

### Plot====
p1 <- ggplot(maas_col, aes(x = coef, y = -log10(qval), color = direction)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = c(
    "Higher in Healthy" = "blue",
    "Higher in Unhealthy" = "red",
    "Not significant" = "grey"
  )) +
  guides(colour = guide_legend(override.aes = list(size = 4))) +
  labs(
    x = "Effect Size (coefficient)",
    y = "-log10(q-value)",
    color = "Direction"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  theme_classic(base_size = 20) +
  theme(axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black")) +
  labs2(tag = "A.") +
  theme(plot.title.position = 'plot', plot.title = ggtext::element_markdown(hjust = 0, vjust = 1, size = 20, margin = margin(t = 0, r = 0, b = 0.8, l = 0, unit = 'mm')))

# Quick look
# p1

# Save volcano plot
# ggsave(height=8, width=10, dpi=300, plot=p1, filename = "ncyc-maaslin2-volcano.png")

#**************************
# Bar plot (functions)----
#**************************
# Join passed ss genes with annotation info
merged_maas <- maas_passed %>% 
  left_join(anno, by = "query") 

# Generate mean coefficient values and stadard deviation of coefficient
# sd generated if more than one gene passed ss with the same functional annotation
top_annotations <- merged_maas %>%
  group_by(gene) %>%
  summarise(
    mean_coef = mean(coef, na.rm = TRUE),
    sd_coef = sd(coef, na.rm = TRUE),
    n = n(),
    se_coef = sd_coef / sqrt(n)
  ) %>%
  arrange(desc(abs(mean_coef)))

# Add higher level annotation information
func_mass <- top_annotations %>% 
  left_join(func, by = "gene")

# Collapse data
collapsed <- func_mass %>%
  group_by(gene) %>%
  summarise(
    mean_coef = first(mean_coef),
    se_coef = first(se_coef),
    pathways = paste(sort(unique(pathway)), collapse = " + ")
  )

# write.table(x = collapsed, file = "maaslin2-ncyc-gene-summary-df.tsv", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

### Plot====
p2 <- ggplot(collapsed, aes(x = reorder(gene, mean_coef), y = mean_coef, fill = pathways)) +
  geom_col() +
  geom_errorbar(aes(ymin = mean_coef - se_coef, ymax = mean_coef + se_coef),
                width = 0.2) +
  coord_flip() +
  scale_fill_brewer(palette = "Set3", name = "Pathway") +
  labs(
    x = "Annotated Gene",
    y = "Mean Effect Size"
  ) +
  theme_classic(base_size = 20) +
  theme(axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black")) +
  labs2(tag = "B.") +
  theme(plot.title.position = 'plot', plot.title = ggtext::element_markdown(hjust = 0, vjust = 1, size = 20, margin = margin(t = 0, r = 0, b = 0.8, l = 0, unit = 'mm')))

# Quick look
# p2

# Save bar plot
# ggsave(height=8, width=10, dpi=300, plot=last_plot(), filename = "ncyc-maaslin2-function-bar-plot.pdf")

#*********************
# Combine plots----
#*********************
(p1 + p2) + plot_layout(guides = "collect") & theme(legend.position = "bottom", legend.title.position = "top", legend.title = element_text(hjust = 0.5)) &
  guides(fill = guide_legend(nrow = 4),
         colour = guide_legend(nrow = 2, override.aes = list(size = 4)))

ggsave(height=14, width=17, dpi=300, plot=last_plot(), filename = "../ncyc-maaslin2-all-plots.png")


## End of Script ##