# Top Performers Genetic Metrics Analysis
# Question: Do top ED50 performers have higher heterozygosity or allelic richness?
# Extension of top_performers_relatedness.r to test additional genetic metrics

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

# =============================================================================
# 2. CALCULATE INDIVIDUAL-LEVEL HETEROZYGOSITY
# =============================================================================

cat("\n=== CALCULATING INDIVIDUAL HETEROZYGOSITY ===\n")

# Calculate observed heterozygosity for each individual
# This is the proportion of heterozygous loci per individual
# Extract genotype matrix (0 = homozygous ref, 1 = het, 2 = homozygous alt, NA = missing)
geno_matrix <- as.matrix(gl_f)

# Calculate proportion of heterozygous sites (value = 1) for each individual
ind_het <- rowSums(geno_matrix == 1, na.rm = TRUE) / rowSums(!is.na(geno_matrix))

# Extract heterozygosity values and create dataframe
het_df <- data.frame(
  id = as.character(indNames(gl_f)),
  heterozygosity = as.numeric(ind_het)
)

cat("Calculated heterozygosity for", nrow(het_df), "individuals\n")
cat("Mean heterozygosity:", round(mean(het_df$heterozygosity, na.rm = TRUE), 4), "\n")
cat(
  "Range:", round(min(het_df$heterozygosity, na.rm = TRUE), 4), "-",
  round(max(het_df$heterozygosity, na.rm = TRUE), 4), "\n"
)

# =============================================================================
# 3. CALCULATE INDIVIDUAL-LEVEL ALLELIC RICHNESS
# =============================================================================

cat("\n=== CALCULATING INDIVIDUAL ALLELIC RICHNESS ===\n")

# Allelic richness per individual = number of unique alleles carried
# For biallelic SNPs, each locus can contribute 0, 1, or 2 alleles
# We'll calculate the proportion of loci where the individual carries at least one copy of each allele

# Calculate allelic richness as proportion of heterozygous sites
# (sites where individual carries both alleles)
# This is identical to heterozygosity for biallelic markers
# For a more meaningful metric, we'll calculate the mean number of non-reference alleles
ind_allelic_richness <- rowMeans(geno_matrix, na.rm = TRUE)

allelic_df <- data.frame(
  id = as.character(rownames(geno_matrix)),
  allelic_richness = ind_allelic_richness
)

cat("Calculated allelic richness for", nrow(allelic_df), "individuals\n")
cat("Mean allelic richness:", round(mean(allelic_df$allelic_richness, na.rm = TRUE), 4), "\n")
cat(
  "Range:", round(min(allelic_df$allelic_richness, na.rm = TRUE), 4), "-",
  round(max(allelic_df$allelic_richness, na.rm = TRUE), 4), "\n"
)

# =============================================================================
# 4. LOAD ED50 DATA AND MATCH TO GENETIC IDS
# =============================================================================

# Load ED50 data
ed50_data <- read.csv("Phenotype Analysis/3. Thermal Curves/outputs/FqFm_T3/ed50.csv")

# Load ID mapping
id_mapping <- read.csv("Genetic Analysis/ind_metrics_site.csv")

# Clean up Sample IDs for matching
id_mapping$Sample.ID <- trimws(id_mapping$Sample.ID)
ed50_data$SampleID <- trimws(ed50_data$SampleID)

# Create combined dataset with all metrics
combined_data <- id_mapping %>%
  select(id, Sample.ID, pop) %>%
  mutate(id = as.character(id)) %>%
  left_join(ed50_data %>% select(SampleID, ed50, Group),
    by = c("Sample.ID" = "SampleID")
  ) %>%
  left_join(het_df, by = "id") %>%
  left_join(allelic_df, by = "id") %>%
  filter(id %in% as.character(genetic_ids)) %>%
  filter(!is.na(ed50))

cat("\n=== DATA MATCHING SUMMARY ===\n")
cat("Individuals in genetic data:", length(genetic_ids), "\n")
cat("Individuals with ED50 matched:", nrow(combined_data), "\n")

# =============================================================================
# 5. DEFINE TOP PERFORMERS
# =============================================================================

# Calculate quantiles
quantiles <- quantile(combined_data$ed50, probs = c(0.50, 0.80, 0.90))
cat("\nED50 quantiles:\n")
print(quantiles)

# Define top performers as top 20% (above 80th percentile)
top_threshold <- quantiles["80%"]
combined_data$performance_group <- ifelse(
  combined_data$ed50 >= top_threshold,
  "Top 20%",
  "Bottom 80%"
)

top_performers <- combined_data %>%
  filter(ed50 >= top_threshold)

cat("\nTop performers (ED50 >=", round(top_threshold, 2), "°C):", nrow(top_performers), "individuals\n")

# =============================================================================
# 6. TEST: HETEROZYGOSITY
# =============================================================================

cat("\n")
cat("=" %>% rep(60) %>% paste(collapse = ""), "\n")
cat("HETEROZYGOSITY ANALYSIS\n")
cat("=" %>% rep(60) %>% paste(collapse = ""), "\n")

# Compare mean heterozygosity between top performers and others
het_top <- combined_data %>% filter(performance_group == "Top 20%")
het_bottom <- combined_data %>% filter(performance_group == "Bottom 80%")

cat("\nObserved mean heterozygosity:\n")
cat("  Top 20%:", round(mean(het_top$heterozygosity, na.rm = TRUE), 4), "\n")
cat("  Bottom 80%:", round(mean(het_bottom$heterozygosity, na.rm = TRUE), 4), "\n")

# Permutation test: Are top performers more heterozygous than expected?
set.seed(42)
n_permutations <- 10000
n_top <- nrow(top_performers)

observed_het_mean <- mean(het_top$heterozygosity, na.rm = TRUE)

null_het_means <- replicate(n_permutations, {
  random_sample <- sample_n(combined_data, n_top)
  mean(random_sample$heterozygosity, na.rm = TRUE)
})

# Two-tailed test: are top performers different from random?
p_value_het <- mean(abs(null_het_means - mean(combined_data$heterozygosity, na.rm = TRUE)) >=
  abs(observed_het_mean - mean(combined_data$heterozygosity, na.rm = TRUE)))

cat("\nPermutation test results (n =", n_permutations, "):\n")
cat("  Observed mean (top 20%):", round(observed_het_mean, 4), "\n")
cat("  Null distribution mean:", round(mean(null_het_means), 4), "\n")
cat("  Null distribution SD:", round(sd(null_het_means), 4), "\n")
cat("  P-value (two-tailed):", round(p_value_het, 4), "\n")

if (p_value_het < 0.05) {
  cat("\n*** SIGNIFICANT: Top performers have different heterozygosity than expected ***\n")
} else {
  cat("\n*** NOT SIGNIFICANT: Top performers' heterozygosity is not different from random ***\n")
}

# Visualize
het_null_df <- data.frame(heterozygosity = null_het_means)

p_het <- ggplot(het_null_df, aes(x = heterozygosity)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7, color = "white") +
  geom_vline(xintercept = observed_het_mean, color = "red", linewidth = 1.2, linetype = "dashed") +
  annotate("text",
    x = observed_het_mean, y = Inf,
    label = paste("Observed\np =", round(p_value_het, 3)),
    vjust = 2, hjust = -0.1, color = "red", size = 4
  ) +
  labs(
    title = "Do Top ED50 Performers Have Different Heterozygosity?",
    subtitle = paste("Top 20% (n =", n_top, ") vs random samples of same size"),
    x = "Mean Observed Heterozygosity",
    y = "Frequency"
  ) +
  default_theme

ggsave("Genetic Analysis/outputs/top_performers_heterozygosity.png", p_het,
  width = 8, height = 6, dpi = 300
)

# =============================================================================
# 7. TEST: ALLELIC RICHNESS
# =============================================================================

cat("\n")
cat("=" %>% rep(60) %>% paste(collapse = ""), "\n")
cat("ALLELIC RICHNESS ANALYSIS\n")
cat("=" %>% rep(60) %>% paste(collapse = ""), "\n")

# Compare mean allelic richness between top performers and others
ar_top <- combined_data %>% filter(performance_group == "Top 20%")
ar_bottom <- combined_data %>% filter(performance_group == "Bottom 80%")

cat("\nObserved mean allelic richness:\n")
cat("  Top 20%:", round(mean(ar_top$allelic_richness, na.rm = TRUE), 4), "\n")
cat("  Bottom 80%:", round(mean(ar_bottom$allelic_richness, na.rm = TRUE), 4), "\n")

# Permutation test
observed_ar_mean <- mean(ar_top$allelic_richness, na.rm = TRUE)

null_ar_means <- replicate(n_permutations, {
  random_sample <- sample_n(combined_data, n_top)
  mean(random_sample$allelic_richness, na.rm = TRUE)
})

# Two-tailed test
p_value_ar <- mean(abs(null_ar_means - mean(combined_data$allelic_richness, na.rm = TRUE)) >=
  abs(observed_ar_mean - mean(combined_data$allelic_richness, na.rm = TRUE)))

cat("\nPermutation test results (n =", n_permutations, "):\n")
cat("  Observed mean (top 20%):", round(observed_ar_mean, 4), "\n")
cat("  Null distribution mean:", round(mean(null_ar_means), 4), "\n")
cat("  Null distribution SD:", round(sd(null_ar_means), 4), "\n")
cat("  P-value (two-tailed):", round(p_value_ar, 4), "\n")

if (p_value_ar < 0.05) {
  cat("\n*** SIGNIFICANT: Top performers have different allelic richness than expected ***\n")
} else {
  cat("\n*** NOT SIGNIFICANT: Top performers' allelic richness is not different from random ***\n")
}

# Visualize
ar_null_df <- data.frame(allelic_richness = null_ar_means)

p_ar <- ggplot(ar_null_df, aes(x = allelic_richness)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7, color = "white") +
  geom_vline(xintercept = observed_ar_mean, color = "red", linewidth = 1.2, linetype = "dashed") +
  annotate("text",
    x = observed_ar_mean, y = Inf,
    label = paste("Observed\np =", round(p_value_ar, 3)),
    vjust = 2, hjust = -0.1, color = "red", size = 4
  ) +
  labs(
    title = "Do Top ED50 Performers Have Different Allelic Richness?",
    subtitle = paste("Top 20% (n =", n_top, ") vs random samples of same size"),
    x = "Mean Allelic Richness (Mean Alt Alleles per Locus)",
    y = "Frequency"
  ) +
  default_theme

ggsave("Genetic Analysis/outputs/top_performers_allelic_richness.png", p_ar,
  width = 8, height = 6, dpi = 300
)

# =============================================================================
# 8. CORRELATION ANALYSES
# =============================================================================

cat("\n")
cat("=" %>% rep(60) %>% paste(collapse = ""), "\n")
cat("CORRELATION ANALYSES\n")
cat("=" %>% rep(60) %>% paste(collapse = ""), "\n")

# ED50 vs Heterozygosity
cor_het <- cor.test(combined_data$ed50, combined_data$heterozygosity,
  method = "spearman", exact = FALSE
)
cat("\nED50 vs Heterozygosity (Spearman correlation):\n")
cat("  rho =", round(cor_het$estimate, 4), "\n")
cat("  p-value =", format.pval(cor_het$p.value, digits = 4), "\n")

# ED50 vs Allelic Richness
cor_ar <- cor.test(combined_data$ed50, combined_data$allelic_richness,
  method = "spearman", exact = FALSE
)
cat("\nED50 vs Allelic Richness (Spearman correlation):\n")
cat("  rho =", round(cor_ar$estimate, 4), "\n")
cat("  p-value =", format.pval(cor_ar$p.value, digits = 4), "\n")

# Heterozygosity vs Allelic Richness (sanity check - should be highly correlated)
cor_het_ar <- cor.test(combined_data$heterozygosity, combined_data$allelic_richness,
  method = "spearman", exact = FALSE
)
cat("\nHeterozygosity vs Allelic Richness (Spearman correlation):\n")
cat("  rho =", round(cor_het_ar$estimate, 4), "\n")
cat("  p-value =", format.pval(cor_het_ar$p.value, digits = 4), "\n")

# Plot ED50 vs Heterozygosity
p_cor_het <- ggplot(combined_data, aes(x = heterozygosity, y = ed50)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  labs(
    title = "ED50 vs Individual Heterozygosity",
    subtitle = paste(
      "Spearman rho =", round(cor_het$estimate, 3),
      ", p =", format.pval(cor_het$p.value, digits = 3)
    ),
    x = "Observed Heterozygosity",
    y = "ED50 (°C)"
  ) +
  default_theme

ggsave("Genetic Analysis/outputs/ed50_vs_heterozygosity.png", p_cor_het,
  width = 8, height = 6, dpi = 300
)

# Plot ED50 vs Allelic Richness
p_cor_ar <- ggplot(combined_data, aes(x = allelic_richness, y = ed50)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  labs(
    title = "ED50 vs Individual Allelic Richness",
    subtitle = paste(
      "Spearman rho =", round(cor_ar$estimate, 3),
      ", p =", format.pval(cor_ar$p.value, digits = 3)
    ),
    x = "Allelic Richness (Mean Alt Alleles per Locus)",
    y = "ED50 (°C)"
  ) +
  default_theme

ggsave("Genetic Analysis/outputs/ed50_vs_allelic_richness.png", p_cor_ar,
  width = 8, height = 6, dpi = 300
)

# =============================================================================
# 9. SENSITIVITY ANALYSIS: Different ED50 thresholds
# =============================================================================

cat("\n")
cat("=" %>% rep(60) %>% paste(collapse = ""), "\n")
cat("SENSITIVITY ANALYSIS: Different ED50 thresholds\n")
cat("=" %>% rep(60) %>% paste(collapse = ""), "\n")

thresholds <- c(0.50, 0.80, 0.90)
sensitivity_results <- data.frame(
  threshold = character(),
  n_individuals = integer(),
  mean_het = numeric(),
  p_het = numeric(),
  mean_ar = numeric(),
  p_ar = numeric(),
  stringsAsFactors = FALSE
)

for (thresh in thresholds) {
  cutoff <- quantile(combined_data$ed50, probs = thresh)
  top_group <- combined_data %>% filter(ed50 >= cutoff)

  if (nrow(top_group) >= 2) {
    # Heterozygosity test
    obs_het <- mean(top_group$heterozygosity, na.rm = TRUE)
    null_het <- replicate(1000, {
      mean(sample_n(combined_data, nrow(top_group))$heterozygosity, na.rm = TRUE)
    })
    p_het <- mean(abs(null_het - mean(combined_data$heterozygosity, na.rm = TRUE)) >=
      abs(obs_het - mean(combined_data$heterozygosity, na.rm = TRUE)))

    # Allelic richness test
    obs_ar <- mean(top_group$allelic_richness, na.rm = TRUE)
    null_ar <- replicate(1000, {
      mean(sample_n(combined_data, nrow(top_group))$allelic_richness, na.rm = TRUE)
    })
    p_ar <- mean(abs(null_ar - mean(combined_data$allelic_richness, na.rm = TRUE)) >=
      abs(obs_ar - mean(combined_data$allelic_richness, na.rm = TRUE)))

    sensitivity_results <- rbind(sensitivity_results, data.frame(
      threshold = paste0("Top ", (1 - thresh) * 100, "%"),
      n_individuals = nrow(top_group),
      mean_het = round(obs_het, 4),
      p_het = round(p_het, 4),
      mean_ar = round(obs_ar, 4),
      p_ar = round(p_ar, 4)
    ))
  }
}

print(sensitivity_results)

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n")
cat("=" %>% rep(60) %>% paste(collapse = ""), "\n")
cat("ANALYSIS COMPLETE\n")
cat("=" %>% rep(60) %>% paste(collapse = ""), "\n")
cat("\nKey outputs:\n")
cat("1. Heterozygosity analysis: Do top ED50 performers have higher heterozygosity?\n")
cat("2. Allelic richness analysis: Do top ED50 performers have higher allelic richness?\n")
cat("3. Correlation tests: Relationship between ED50 and genetic metrics\n")
cat("4. Sensitivity analysis across different thresholds\n")
cat("\nFigures saved to Genetic Analysis/outputs/\n")
cat("  - top_performers_heterozygosity.png\n")
cat("  - top_performers_allelic_richness.png\n")
cat("  - ed50_vs_heterozygosity.png\n")
cat("  - ed50_vs_allelic_richness.png\n")
