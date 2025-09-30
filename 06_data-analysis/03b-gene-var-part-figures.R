## ---------------------------
##
## Script name: gene-var-part-figures.R
##
## Purpose of script: Plot varition partitioning analysis
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
library(vegan)

### Functions====
fig_label <- function(text, region="figure", pos="topleft", cex=NULL, ...) {
  
  region <- match.arg(region, c("figure", "plot", "device"))
  pos <- match.arg(pos, c("topleft", "top", "topright", 
                          "left", "center", "right", 
                          "bottomleft", "bottom", "bottomright"))
  
  if(region %in% c("figure", "device")) {
    ds <- dev.size("in")
    # xy coordinates of device corners in user coordinates
    x <- grconvertX(c(0, ds[1]), from="in", to="user")
    y <- grconvertY(c(0, ds[2]), from="in", to="user")
    
    # fragment of the device we use to plot
    if(region == "figure") {
      # account for the fragment of the device that 
      # the figure is using
      fig <- par("fig")
      dx <- (x[2] - x[1])
      dy <- (y[2] - y[1])
      x <- x[1] + dx * fig[1:2]
      y <- y[1] + dy * fig[3:4]
    } 
  }
  
  # much simpler if in plotting region
  if(region == "plot") {
    u <- par("usr")
    x <- u[1:2]
    y <- u[3:4]
  }
  
  sw <- strwidth(text, cex=cex) * 60/100
  sh <- strheight(text, cex=cex) * 60/100
  
  x1 <- switch(pos,
               topleft     =x[1] + sw, 
               left        =x[1] + sw,
               bottomleft  =x[1] + sw,
               top         =(x[1] + x[2])/2,
               center      =(x[1] + x[2])/2,
               bottom      =(x[1] + x[2])/2,
               topright    =x[2] - sw,
               right       =x[2] - sw,
               bottomright =x[2] - sw)
  
  y1 <- switch(pos,
               topleft     =y[2] - sh,
               top         =y[2] - sh,
               topright    =y[2] - sh,
               left        =(y[1] + y[2])/2,
               center      =(y[1] + y[2])/2,
               right       =(y[1] + y[2])/2,
               bottomleft  =y[1] + sh,
               bottom      =y[1] + sh,
               bottomright =y[1] + sh)
  
  old.par <- par(xpd=NA)
  on.exit(par(old.par))
  
  text(x1, y1, text, cex=cex, ...)
  return(invisible(c(x,y)))
}

#***********
# Data----
#***********
varven_egg <- readRDS("rds-objects/eggnog-bray-varpart.rds") 
varven_ncyc <- readRDS("rds-objects/ncyc-bray-varpart.rds")
varven_cazy <- readRDS("rds-objects/cazy-bray-varpart.rds")

#****************
# eggNOG plot----
#****************
pdf(file = "../eggnog-bray-var-part-env.pdf", height = 8, width = 12.5)
par(
  mai = c(1, 0.1, 0.1, 0.1), 
  bty = "n", # Remove box around plot
  xpd = NA, # Disable clipping of text 
  oma=c(0,0,1,0)
)

plot(varven_egg, bg = c("forestgreen","skyblue", "orange"), Xnames = NA, cex = 2)
text(-0.85, 0.68, "Physicochemical \n properties", cex = 2)
text(1.7, 0.7, "Site and Plot", cex = 2)
text(0.5, -1.65, "Tree health", cex = 2)
fig_label("A.", region = "figure", pos = "topleft", cex = 2.5) # Label may not appear when view in plot pane but should appear when saved
dev.off()

#****************
# NCyc plot----
#****************
pdf(file = "../ncyc-bray-var-part-env.pdf", height = 8, width = 12.5)
par(
  mai = c(1, 0.1, 0.1, 0.1), 
  bty = "n", # Remove box around plot
  xpd = NA, # Disable clipping of text 
  oma=c(0,0,1,0)
)

plot(varven_ncyc, bg = c("forestgreen","skyblue", "orange"), Xnames = NA, cex = 2)
text(-0.85, 0.68, "Physicochemical \n properties", cex = 2)
text(1.7, 0.7, "Site and Plot", cex = 2)
text(0.5, -1.65, "Tree health", cex = 2)
fig_label("A.", region = "figure", pos = "topleft", cex = 2.5)
dev.off()

#****************
# CAZy plot----
#****************
pdf(file = "../cazy-bray-var-part-env.pdf", height = 8, width = 12.5)
par(
  mai = c(1, 0.1, 0.1, 0.1), 
  bty = "n", # Remove box around plot
  xpd = NA, # Disable clipping of text 
  oma=c(0,0,1,0)
)

plot(varven_cazy, bg = c("forestgreen","skyblue", "orange"), Xnames = NA, cex = 2)
text(-0.85, 0.68, "Physicochemical \n properties", cex = 2)
text(1.7, 0.7, "Site and Plot", cex = 2)
text(0.5, -1.65, "Tree health", cex = 2)
fig_label("A.", region = "figure", pos = "topleft", cex = 2.5)
dev.off()

## End of Script ##