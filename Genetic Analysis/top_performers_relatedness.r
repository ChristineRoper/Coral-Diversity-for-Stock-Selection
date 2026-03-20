# Top Performers Relatedness Analysis
# Question: Are the top ED50 performers more closely related to each other than expected by chance?

library(dartR)
library(dplyr)
library(ggplot2)

source("default_theme.r")

# =============================================================================
# 1. LOAD AND PREPARE DATA
# =============================================================================

# Load from "gl_f.RData"
# Note: If gl_f.RData does not exist, run create_gl_f.r
load("gl_f.RData")

# Impute missing data (needed for distance calculation)
gl_f2 <- gl.impute(gl_f, method = "neighbour", verbose = 5)
gl_f2 <- gl.recalc.metrics(gl_f2)

# =============================================================================
# 2. CALCULATE INDIVIDUAL-LEVEL GENETIC DISTANCE
# =============================================================================

# Calculate pairwise genetic distance between all individuals
ind_dist <- gl.dist.ind(gl_f2, method = "euclidean")
ind_dist_matrix <- as.matrix(ind_dist)

# Get individual IDs from genlight object
genetic_ids <- indNames(gl_f2)

# =============================================================================
# 3. LOAD ED50 DATA AND MATCH TO GENETIC IDS
# =============================================================================

# Load ED50 data
ed50_data <- read.csv("Phenotype Analysis/3. Thermal Curves/outputs/FqFm_T3/ed50.csv")

# Load ID mapping
id_mapping <- read.csv("Genetic Analysis/ind_metrics_site.csv")

# The genetic data uses numeric IDs (1, 2, 3...)
# The ED50 data uses Sample IDs like "BH 1", "LB 2", etc.
# Need to match them via the mapping file

# Clean up Sample IDs for matching (remove extra spaces, standardize)
id_mapping$Sample.ID <- trimws(id_mapping$Sample.ID)
ed50_data$SampleID <- trimws(ed50_data$SampleID)

# Create lookup: genetic ID -> Sample ID -> ED50
ed50_lookup <- id_mapping %>%
  select(id, Sample.ID) %>%
  left_join(ed50_data %>% select(SampleID, ed50, Group),
    by = c("Sample.ID" = "SampleID")
  )

# Filter to only individuals that are in our filtered genetic dataset
ed50_matched <- ed50_lookup %>%
  filter(as.character(id) %in% genetic_ids) %>%
  filter(!is.na(ed50))

cat("\n=== DATA MATCHING SUMMARY ===\n")
cat("Individuals in genetic data:", length(genetic_ids), "\n")
cat("Individuals with ED50 matched:", nrow(ed50_matched), "\n")

# =============================================================================
# 4. DEFINE TOP PERFORMERS
# =============================================================================

# Calculate quantiles
quantiles <- quantile(ed50_matched$ed50, probs = c(0.50, 0.80, 0.90))
cat("\nED50 quantiles:\n")
print(quantiles)

# Define top performers as top 20% (above 80th percentile)
top_threshold <- quantiles["80%"]
top_performers <- ed50_matched %>%
  filter(ed50 >= top_threshold)

cat("\nTop performers (ED50 >=", round(top_threshold, 2), "°C):", nrow(top_performers), "individuals\n")
cat("ED50 range of top performers:", round(min(top_performers$ed50), 2), "-", round(max(top_performers$ed50), 2), "°C\n")

# =============================================================================
# 5. CALCULATE WITHIN-GROUP GENETIC DISTANCES
# =============================================================================

# Function to get mean pairwise distance for a group of individuals
get_mean_within_distance <- function(ind_ids, dist_matrix) {
  if (length(ind_ids) < 2) {
    return(NA)
  }

  # Subset the distance matrix
  sub_matrix <- dist_matrix[as.character(ind_ids), as.character(ind_ids)]

  # Get upper triangle (exclude diagonal and duplicates)
  upper_tri <- sub_matrix[upper.tri(sub_matrix)]

  return(mean(upper_tri, na.rm = TRUE))
}

# Calculate observed mean distance for top performers
top_ids <- as.character(top_performers$id)
observed_mean_dist <- get_mean_within_distance(top_ids, ind_dist_matrix)

cat("\n=== OBSERVED DISTANCES ===\n")
cat("Mean genetic distance among top 20% performers:", round(observed_mean_dist, 4), "\n")

# =============================================================================
# 6. PERMUTATION TEST
# =============================================================================

# How likely is it to observe this (or smaller) distance by chance?
# Randomly sample same number of individuals many times

set.seed(42)
n_permutations <- 10000
n_top <- length(top_ids)
all_ids <- as.character(ed50_matched$id)

# Generate null distribution
null_distances <- numeric(n_permutations)

for (i in 1:n_permutations) {
  random_ids <- sample(all_ids, n_top, replace = FALSE)
  null_distances[i] <- get_mean_within_distance(random_ids, ind_dist_matrix)
}

# Calculate p-value (two-tailed: are top performers MORE or LESS related than random?)
null_mean <- mean(null_distances)
p_value <- mean(abs(null_distances - null_mean) >= abs(observed_mean_dist - null_mean))

cat("\n=== PERMUTATION TEST RESULTS (n =", n_permutations, ") ===\n")
cat("Observed mean distance (top 20%):", round(observed_mean_dist, 4), "\n")
cat("Null distribution mean:", round(mean(null_distances), 4), "\n")
cat("Null distribution SD:", round(sd(null_distances), 4), "\n")
cat("P-value (two-tailed, testing if top performers are more or less related):", round(p_value, 4), "\n")

if (p_value < 0.05) {
  cat("\n*** SIGNIFICANT: Top performers are more genetically related than expected by chance ***\n")
} else {
  cat("\n*** NOT SIGNIFICANT: Top performers are not more related than random individuals ***\n")
}

# =============================================================================
# 7. VISUALIZE RESULTS
# =============================================================================

# Plot null distribution with observed value
null_df <- data.frame(distance = null_distances)

p1 <- ggplot(null_df, aes(x = distance)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7, color = "white") +
  geom_vline(xintercept = observed_mean_dist, color = "red", linewidth = 1.2, linetype = "dashed") +
  annotate("text",
    x = observed_mean_dist, y = Inf,
    label = paste("Observed\np =", round(p_value, 3)),
    vjust = 2, hjust = -0.1, color = "red", size = 4
  ) +
  labs(
    title = "Are Top ED50 Performers More Genetically Related?",
    subtitle = paste("Top 20% (n =", n_top, ") vs random samples of same size"),
    x = "Mean Pairwise Genetic Distance (Euclidean)",
    y = "Frequency"
  ) +
  default_theme

ggsave("Genetic Analysis/outputs/top_performers_relatedness.png", p1,
  width = 8, height = 6, dpi = 300
)
# ggsave("Genetic Analysis/outputs/top_performers_relatedness.svg", p1,
#        width = 8, height = 6)

cat("\nPlot saved to: Genetic Analysis/outputs/top_performers_relatedness.png\n")

# =============================================================================
# 8. ADDITIONAL ANALYSIS: By Reef/Population
# =============================================================================

# Check if top performers cluster within specific populations
cat("\n=== POPULATION BREAKDOWN OF TOP PERFORMERS ===\n")
pop_breakdown <- top_performers %>%
  group_by(Group) %>%
  summarise(
    n_top = n(),
    mean_ed50 = mean(ed50),
    .groups = "drop"
  ) %>%
  arrange(desc(n_top))

# Compare to overall population sizes
overall_pop <- ed50_matched %>%
  group_by(Group) %>%
  summarise(n_total = n(), .groups = "drop")

pop_comparison <- left_join(overall_pop, pop_breakdown, by = "Group") %>%
  mutate(
    n_top = ifelse(is.na(n_top), 0, n_top),
    pct_top = round(n_top / n_total * 100, 1)
  ) %>%
  arrange(desc(pct_top))

print(pop_comparison)

# =============================================================================
# 9. SENSITIVITY ANALYSIS: Different thresholds
# =============================================================================

cat("\n=== SENSITIVITY ANALYSIS: Different ED50 thresholds ===\n")

thresholds <- c(0.50, 0.80, 0.90)
sensitivity_results <- data.frame(
  threshold = character(),
  n_individuals = integer(),
  observed_dist = numeric(),
  null_mean = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for (thresh in thresholds) {
  cutoff <- quantile(ed50_matched$ed50, probs = thresh)
  top_ids_sens <- as.character(ed50_matched$id[ed50_matched$ed50 >= cutoff])

  if (length(top_ids_sens) >= 2) {
    obs_dist <- get_mean_within_distance(top_ids_sens, ind_dist_matrix)

    # Quick permutation test
    null_sens <- replicate(1000, {
      random_ids <- sample(all_ids, length(top_ids_sens), replace = FALSE)
      get_mean_within_distance(random_ids, ind_dist_matrix)
    })

    p_val <- mean(null_sens <= obs_dist)

    sensitivity_results <- rbind(sensitivity_results, data.frame(
      threshold = paste0("Top ", (1 - thresh) * 100, "%"),
      n_individuals = length(top_ids_sens),
      observed_dist = round(obs_dist, 4),
      null_mean = round(mean(null_sens), 4),
      p_value = round(p_val, 4)
    ))
  }
}

print(sensitivity_results)

# =============================================================================
# 10. ED50 vs Genetic Distance Correlation (all individuals)
# =============================================================================

cat("\n=== ED50 DIFFERENCES vs GENETIC DISTANCE ===\n")

# Create all pairwise comparisons
pairs_df <- expand.grid(
  id1 = ed50_matched$id,
  id2 = ed50_matched$id
) %>%
  filter(id1 < id2) %>% # Keep only unique pairs
  left_join(ed50_matched %>% select(id, ed50_1 = ed50), by = c("id1" = "id")) %>%
  left_join(ed50_matched %>% select(id, ed50_2 = ed50), by = c("id2" = "id"))

# Add genetic distances
pairs_df$genetic_dist <- mapply(function(i, j) {
  ind_dist_matrix[as.character(i), as.character(j)]
}, pairs_df$id1, pairs_df$id2)

# Calculate ED50 difference
pairs_df$ed50_diff <- abs(pairs_df$ed50_1 - pairs_df$ed50_2)

# Correlation test
cor_test <- cor.test(pairs_df$genetic_dist, pairs_df$ed50_diff, method = "spearman")
cat("Spearman correlation (genetic distance vs ED50 difference):\n")
cat("  rho =", round(cor_test$estimate, 4), "\n")
cat("  p-value =", format.pval(cor_test$p.value, digits = 4), "\n")

# If significant and positive: genetically similar individuals have similar ED50s
# If not significant: ED50 is independent of genetic relatedness

# Plot
p2 <- ggplot(pairs_df, aes(x = genetic_dist, y = ed50_diff)) +
  geom_point(alpha = 0.1, size = 0.5) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  labs(
    title = "Genetic Distance vs ED50 Phenotypic Difference",
    subtitle = paste(
      "Spearman rho =", round(cor_test$estimate, 3),
      ", p =", format.pval(cor_test$p.value, digits = 3)
    ),
    x = "Pairwise Genetic Distance (Euclidean)",
    y = "ED50 Difference (°C)"
  ) +
  default_theme

ggsave("Genetic Analysis/outputs/ed50_vs_genetic_distance.png", p2,
  width = 8, height = 6, dpi = 300
)
# ggsave("Genetic Analysis/outputs/ed50_vs_genetic_distance.svg", p2,
#        width = 8, height = 6)

cat("\nPlot saved to: Genetic Analysis/outputs/ed50_vs_genetic_distance.png\n")

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n")
cat("=" %>% rep(60) %>% paste(collapse = ""), "\n")
cat("ANALYSIS COMPLETE\n")
cat("=" %>% rep(60) %>% paste(collapse = ""), "\n")
cat("\nKey outputs:\n")
cat("1. Permutation test: Do top ED50 performers cluster genetically?\n")
cat("2. Sensitivity analysis across different thresholds\n")
cat("3. Correlation: Do genetically similar individuals have similar ED50?\n")
cat("\nFigures saved to Genetic Analysis/outputs/\n")
