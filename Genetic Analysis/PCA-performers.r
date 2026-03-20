library(readr)
library(dplyr)
library(dartR)

# Load from "gl_f.RData"
# Note: If gl_f.RData does not exist, run create_gl_f.r
load("gl_f.RData")

# Load phenotypic performance data
sample_to_numeric_ids <- read_csv("Genetic Analysis/ind_metrics_site.csv") %>%
  select(id, "Sample ID") %>%
  rename(SampleID = "Sample ID")


ed50s <- read_csv("Phenotype Analysis/3. Thermal Curves/outputs/FqFm_T3/ed50.csv") %>%
  left_join(sample_to_numeric_ids, by = "SampleID")

# View(sample_to_numeric_ids %>% left_join(read_csv("Phenotype Analysis/3. Thermal Curves/outputs/FqFm_T3/ed50.csv"), by = "SampleID"))

top20 <- quantile(ed50s$ed50, 0.8)
bottom20 <- quantile(ed50s$ed50, 0.2)

ed50s <- ed50s %>%
  mutate(id = as.numeric(id), pop = case_when(
    ed50 > top20 ~ "top",
    ed50 < bottom20 ~ "bottom",
    TRUE ~ "middle"
  )) %>%
  select(id, pop) %>%
  filter(!is.na(id))

pop(gl_f) <- ed50s$pop[match(indNames(gl_f), ed50s$id)]

# PCA
library(adegenet)
source("default_theme.r")

pcoa <- gl.pcoa(gl_f, verbose = 3)

gl.pcoa.plot(pcoa, gl_f, pop.labels = "none", ellipse = FALSE)
plot <- ggplot2::last_plot()

# Remove the horizontal and vertical lines
plot$layers[[3]] <- NULL
plot$layers[[2]] <- NULL


# output the plot
plot +
  # stat_ellipse(type = "norm", level = 0.95) +
  scale_color_manual(
    name = "ED50 Performance",
    values = c("top" = "#F84848", "middle" = "#aeaeae", "bottom" = "#66BEDA"),
    labels = c("top" = "Top 20%", "bottom" = "Bottom 20%"),
    breaks = c("top", "bottom"),
    na.value = "#aeaeae"
  ) +
  xlim(-3, 8) + ylim(-5, 5) +
  default_theme +
  theme(
    legend.position = c(0.95, 0.95),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", color = NA),
    legend.text = element_text(size = rel(1.3)),
    legend.title = element_text(size = rel(1.3))
  )

ggsave("Genetic Analysis/plots/pca-performers.png", width = 6, height = 5)
