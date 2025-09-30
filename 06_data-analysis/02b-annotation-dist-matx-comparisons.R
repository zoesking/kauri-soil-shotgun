## ---------------------------
##
## Script name: annotation-dist-matx-comparisons.R
##
## Purpose of script:
##
## Author: Zoe King
##
## Date Created: 2025-09-16
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
### Functions====

#***********
# Data----
#***********
egg <- readRDS("gene-catalog-data/R-scripts/beta-diversity/eggnog-bc-distance-matrix.rds")
ncyc <- readRDS("gene-catalog-data/R-scripts/beta-diversity/ncyc-bc-distance-matrix.rds")
cazy <- readRDS("gene-catalog-data/R-scripts/beta-diversity/cazy-bc-distance-matrix.rds")

set.seed(413)

#**************************
# Spearman correlation----
#**************************
# Helper function to extract lower triangle as a vector
lt <- function(d) as.vector(as.matrix(d)[lower.tri(as.matrix(d))])

eggNOG_vs_NCyc <- cor.test(lt(egg), lt(ncyc), method = "spearman")
eggNOG_vs_CAZy <- cor.test(lt(egg), lt(cazy), method = "spearman")
NCyc_vs_CAZy <- cor.test(lt(ncyc), lt(cazy), method = "spearman")

#****************
# Mantel test----
#****************
mantel_e_n <- mantel(egg, ncyc, method = "spearman", permutations = 9999)
mantel_e_c <- mantel(egg, cazy, method = "spearman", permutations = 9999)
mantel_n_c <- mantel(ncyc, cazy, method = "spearman", permutations = 9999)

#****************
# Procrustes-----
#****************
pe_n <- protest(egg, ncyc, permutations = 9999)
pe_c <- protest(egg, cazy, permutations = 9999)
pn_c <- protest(ncyc, cazy, permutations = 9999)


#************************************
# Distance-distance scatterplot----
#************************************
mk_df <- function(d1, d2, lab){
  data.frame(x = lt(d1), y = lt(d2), pair = lab)
}
plot_df <- bind_rows(
  mk_df(egg,  ncyc, "eggNOG vs NCyc"),
  mk_df(egg,  cazy, "eggNOG vs CAZy"),
  mk_df(ncyc, cazy, "NCyc vs CAZy")
)

ggplot(plot_df, aes(x, y)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  facet_wrap(~ pair, scales = "free") +
  labs(x = "Distance (method A)", y = "Distance (method B)") +
  theme_classic(base_size = 18) +
  theme(axis.ticks = element_line(colour = "black"), axis.text = element_text(colour = "black"))

## End of Script ##