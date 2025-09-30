## ---------------------------
##
## Script name: maaslin2-eggnog-figures.R
##
## Purpose of script:
##
## Author: Zoe King
##
## Date Created: 2025-08-23
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
library(vroom)
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
maas <- vroom::vroom("tpm-tables/maaslin2_eggNOG_output/all_results.tsv")
# Load in maaslin2 results (significant results only)
maas_passed <- read_tsv("tpm-tables/maaslin2_eggNOG_output/significant_results.tsv") %>% 
  rename("query" = "feature")
# Annotated genes information
anno <- vroom::vroom("annotation-tables/eggnog-cds-annotations.tsv") %>% 
  rename("query" = `#query`)

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
# ggsave(height=8, width=10, dpi=300, plot=p1, filename = "eggnog-maaslin2-volcano.png")

#**************************
# Bar plot (functions)----
#**************************
# Join passed ss genes with annotation info
merged_maas <- maas_passed %>% 
  left_join(anno, by = "query") 

# Quick test subset of entry with multiple cog annotations
# test <- merged_maas %>% filter(query == "k141_692732_1")

# Summarise by OG name
# Helper for NULL/NA-safe fallback
`%or%` <- function(x, y) ifelse(is.na(x) | x == "" | x == "-", y, x)

cog_summary <- merged_maas %>%
  # Only extract COG IDs that are explicitly at the root level: ... "COG0125@1|root"
  mutate(
    eggNOG_OGs = !!sym(names(.)[which(names(.) %in% c("eggNOG_OGs", "eggnog_ogs"))][1] %||% "eggNOG_OGs"),
    # grab only IDs followed by @1|root (positive lookahead)
    COG_ids    = str_extract_all(eggNOG_OGs %||% "", "COG\\d{4}(?=@1\\|root)"),
    # de-duplicate within each gene so weights aren't inflated
    COG_ids    = map(COG_ids, unique),
    n_id       = lengths(COG_ids),
    
    # label preference: Description -> Preferred_name -> COG_id
    COG_name_raw = (Description %or% Preferred_name),
    COG_name_raw = ifelse(is.na(COG_name_raw) | COG_name_raw == "" | COG_name_raw == "-",
                          NA_character_, COG_name_raw)
  ) %>%
  filter(n_id > 0) %>%
  unnest_longer(COG_ids, values_to = "COG_id") %>%
  mutate(
    wcoef    = coef / n_id,
    COG_name = coalesce(COG_name_raw, COG_id)
  ) %>%
  # keep functional class if present
  group_by(query, COG_id, COG_name, COG_category, .add = TRUE) %>%
  summarise(gene_coef = sum(wcoef, na.rm = TRUE), .groups = "drop") %>%
  group_by(COG_id, COG_name, COG_category) %>%
  summarise(
    mean_coef = mean(gene_coef, na.rm = TRUE),
    sd_coef   = sd(gene_coef, na.rm = TRUE),
    n         = n(),
    se_coef   = sd_coef / sqrt(n),
    .groups   = "drop"
  ) %>%
  arrange(desc(abs(mean_coef)))

# Filter output to only show 
cog_filtered <- cog_summary %>%
  filter(mean_coef > 0.4 | mean_coef < -1)

write.table(x = cog_summary, file = "maaslin2-eggnog-cog-summary-df.tsv", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

### Plot====
p2 <- ggplot(cog_filtered, aes(x = reorder(COG_id, mean_coef), y = mean_coef, fill = COG_category)) +
  geom_col() +
  geom_errorbar(aes(ymin = mean_coef - se_coef, ymax = mean_coef + se_coef),
                width = 0.2) +
  coord_flip() +
  scale_fill_brewer(palette = "Set3", name = "COG Category") +
  labs(
    x = "COG",
    y = "Mean Effect Size"
  ) +
  theme_classic(base_size = 20) +
  theme(axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black")) +
  labs2(tag = "B.") +
  theme(plot.title.position = 'plot', plot.title = ggtext::element_markdown(hjust = 0, vjust = 1, size = 20, margin = margin(t = 0, r = 0, b = 0.8, l = 0, unit = 'mm')))

# Quick look
# p2

# Save bar plot
# ggsave(height=8, width=10, dpi=300, plot=last_plot(), filename = "eggnog-maaslin2-function-bar-plot.pdf")

#*********************
# Combine plots----
#*********************
# Reccommended to not check plot before saving
p3 <- (p1 + p2) + plot_layout(guides = "collect") & theme(legend.position = "bottom", legend.title.position = "top", legend.title = element_text(hjust = 0.5)) &
  guides(fill = guide_legend(nrow = 2),
         colour = guide_legend(nrow = 2, override.aes = list(size = 4)))

ggsave(height=12, width=17, dpi=300, plot=p3, filename = "eggNOG-maaslin2-all-plots.png")

## End of Script ##