library(readr)
library(dplyr)
library(dartR)
library(adegenet)
source("default_theme.r")

# Load from "gl_f.RData"
# Note: If gl_f.RData does not exist, run create_gl_f.r
load("gl_f.RData")


# Combined filtered data with unfiltered Long Bommie SNPS
# Load the unfiltered data
load("gl.RData")

# Grab Long Bommie from the original, unfiltered Genlight object (gl)
long_bommie <- gl.keep.pop(gl, pop.list = c("Opal Long Bommie"))

# Keep only the loci that are also present in filtered Genlight object (gl_f)
# Otherwise we can't recombine this data
long_bommie <- gl.keep.loc(long_bommie, locNames(gl_f))

# Combine all the filtered data with the unfiltered Long Bommie data
gl_f_with_long_bommie <- rbind(gl_f, long_bommie)

gl_f_with_long_bommie <- gl.compliance.check(gl_f_with_long_bommie) # get back missing metadata that rbind lost

gl_f_with_long_bommie <- gl.recalc.metrics(gl_f_with_long_bommie) # re-calculate locus metadata


# PCA

# Switch between with/without long_bommie
pca_data <- gl_f
# pca_data <- gl_f_with_long_bommie
# # Get rid of the 64 individual outlier
# pca_data <- gl.drop.ind(pca_data, ind.list = c("64"))

# levels(pca_data$pop) = c(
#   "L1 Flat", "L1 Pale" ,"L2 Slope", "Low Isles slope pale",
#   "L2 Flat", "O2 Crest", "O2 Slope", "O3 Flat", "O3 Slope" ,
#   "T1 Slope", "T2 Slope", "O1 Slope"
# )

pca_data <- gl.rename.pop(pca_data, old = "Low Isles flat", new = "L1 Flat")
pca_data <- gl.rename.pop(pca_data, old = "Low Isles flat pale", new = "L1 Pale")
pca_data <- gl.rename.pop(pca_data, old = "Low Isles Woody flat", new = "L2 Flat")
pca_data <- gl.rename.pop(pca_data, old = "Low Isles slope", new = "L2 Slope")
pca_data <- gl.rename.pop(pca_data, old = "Tongue Turtle Bay", new = "T2 Slope")
pca_data <- gl.rename.pop(pca_data, old = "Tongue Blue Hole", new = "T1 Slope")
pca_data <- gl.rename.pop(pca_data, old = "Opal Long Bommie", new = "O1 Slope")
pca_data <- gl.rename.pop(pca_data, old = "Opal Mojo flat", new = "O2 Crest")
pca_data <- gl.rename.pop(pca_data, old = "Opal Mojo slope", new = "O2 Slope")
pca_data <- gl.rename.pop(pca_data, old = "Opal Rayban flat", new = "O3 Flat")
pca_data <- gl.rename.pop(pca_data, old = "Opal Rayban slope", new = "O3 Slope")

pcoa <- gl.pcoa(pca_data, verbose = 3)

# gl.pcoa.plot won't return the plot for us to work with
# We can tell it to save it in a temporary directory though
gl.pcoa.plot(pcoa, pca_data, pop.labels = "pop", ellipse = FALSE)
plot <- ggplot2::last_plot()

# Remove the horizontal and vertical lines
plot$layers[[4]] <- NULL
plot$layers[[3]] <- NULL

# output the plot
plot +
  # stat_ellipse(type = "norm", level = 0.95) +
  xlim(-3, 9) + ylim(-5, 5) +
  default_theme

ggsave("Genetic Analysis/plots/pca.png", width = 6, height = 5)
ggsave("Genetic Analysis/plots/pca_longbommie.png", width = 6, height = 5)
