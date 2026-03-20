# Top Performers Group-Level Heterozygosity Analysis
# Question: Do top ED50 performers have higher group-level heterozygosity than expected?
# Group-level Ho: proportion of heterozygous genotypes per locus, averaged across loci
# Group-level He: expected heterozygosity from allele frequencies (1 - sum(p^2)), averaged across loci

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
# 2. BUILD GENOTYPE MATRIX AND ID INDEX
# =============================================================================

cat("\n=== BUILDING GENOTYPE MATRIX ===\n")

# Genotype matrix: rows = individuals, cols = loci
# Values: 0 = homozygous ref, 1 = het, 2 = homozygous alt, NA = missing
geno_matrix <- as.matrix(gl_f)
rownames(geno_matrix) <- indNames(gl_f)

cat("Genotype matrix dimensions:", nrow(geno_matrix), "individuals x", ncol(geno_matrix), "loci\n")

# =============================================================================
# 3. HELPER FUNCTIONS FOR GROUP-LEVEL HETEROZYGOSITY
# =============================================================================

# Observed heterozygosity for a group of individuals (subset of geno_matrix rows)
# Ho at each locus = proportion of non-missing genotypes that are heterozygous (value == 1)
# Returns mean Ho across loci
group_ho <- function(geno_sub) {
  locus_ho <- colSums(geno_sub == 1, na.rm = TRUE) / colSums(!is.na(geno_sub))
  mean(locus_ho, na.rm = TRUE)
}

# Expected heterozygosity for a group of individuals
# He at each locus = 1 - (p^2 + q^2) where p = freq of ref allele, q = freq of alt allele
# Allele counts from genotype codes: 0 -> 2 ref alleles, 1 -> 1 ref + 1 alt, 2 -> 2 alt alleles
# Returns mean He across loci
group_he <- function(geno_sub) {
  n_alleles <- 2 * colSums(!is.na(geno_sub)) # total alleles per locus
  alt_alleles <- colSums(geno_sub, na.rm = TRUE) # alt allele count (0->0, 1->1, 2->2)
  ref_alleles <- n_alleles - alt_alleles
  q <- alt_alleles / n_alleles # alt allele frequency
  p <- ref_alleles / n_alleles # ref allele frequency
  locus_he <- 1 - (p^2 + q^2)
  mean(locus_he, na.rm = TRUE)
}

# =============================================================================
# 4. LOAD ED50 DATA AND MATCH TO GENETIC IDS
# =============================================================================

# Get individual IDs from genlight object
genetic_ids <- indNames(gl_f)

ed50_data <- read.csv("Phenotype Analysis/3. Thermal Curves/outputs/FqFm_T3/ed50.csv")
id_mapping <- read.csv("Genetic Analysis/ind_metrics_site.csv")

id_mapping$Sample.ID <- trimws(id_mapping$Sample.ID)
ed50_data$SampleID <- trimws(ed50_data$SampleID)

combined_data <- id_mapping %>%
  select(id, Sample.ID, pop) %>%
  mutate(id = as.character(id)) %>%
  left_join(ed50_data %>% select(SampleID, ed50, Group),
    by = c("Sample.ID" = "SampleID")
  ) %>%
  filter(id %in% as.character(genetic_ids)) %>%
  filter(!is.na(ed50))

cat("\n=== DATA MATCHING SUMMARY ===\n")
cat("Individuals in genetic data:", length(genetic_ids), "\n")
cat("Individuals with ED50 matched:", nrow(combined_data), "\n")

# Map combined_data rows to geno_matrix rows
combined_data$geno_row <- match(combined_data$id, rownames(geno_matrix))
stopifnot("Some individuals not found in geno_matrix" = !any(is.na(combined_data$geno_row)))

# =============================================================================
# 5. DEFINE TOP PERFORMERS
# =============================================================================

quantiles <- quantile(combined_data$ed50, probs = c(0.50, 0.80, 0.90))
cat("\nED50 quantiles:\n")
print(quantiles)

top_threshold <- quantiles["80%"]
combined_data$performance_group <- ifelse(
  combined_data$ed50 >= top_threshold,
  "Top 20%",
  "Bottom 80%"
)

top_rows <- combined_data$geno_row[combined_data$performance_group == "Top 20%"]
all_rows <- combined_data$geno_row

n_top <- length(top_rows)
cat("\nTop performers (ED50 >=", round(top_threshold, 2), "°C):", n_top, "individuals\n")

# =============================================================================
# 6. OBSERVED HETEROZYGOSITY (Ho)
# =============================================================================

cat("\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("OBSERVED HETEROZYGOSITY (Ho)\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

observed_ho <- group_ho(geno_matrix[top_rows, ])
overall_ho <- group_ho(geno_matrix[all_rows, ])

cat("\nGroup-level Ho:\n")
cat("  Top 20%:", round(observed_ho, 4), "\n")
cat("  All individuals:", round(overall_ho, 4), "\n")

set.seed(42)
n_permutations <- 10000

null_ho <- replicate(n_permutations, {
  idx <- sample(all_rows, n_top)
  group_ho(geno_matrix[idx, ])
})

p_value_ho <- mean(abs(null_ho - overall_ho) >= abs(observed_ho - overall_ho))

cat("\nPermutation test results (n =", n_permutations, "):\n")
cat("  Observed Ho (top 20%):", round(observed_ho, 4), "\n")
cat("  Null distribution mean:", round(mean(null_ho), 4), "\n")
cat("  Null distribution SD:", round(sd(null_ho), 4), "\n")
cat("  P-value (two-tailed):", round(p_value_ho, 4), "\n")

if (p_value_ho < 0.05) {
  cat("\n*** SIGNIFICANT: Top performers have different Ho than expected ***\n")
} else {
  cat("\n*** NOT SIGNIFICANT: Top performers' Ho is not different from random ***\n")
}

ho_null_df <- data.frame(ho = null_ho)

p_ho <- ggplot(ho_null_df, aes(x = ho)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7, color = "white") +
  geom_vline(xintercept = observed_ho, color = "red", linewidth = 1.2, linetype = "dashed") +
  annotate("text",
    x = observed_ho, y = Inf,
    label = paste("Observed\np =", round(p_value_ho, 3)),
    vjust = 2, hjust = -0.1, color = "red", size = 4
  ) +
  labs(
    title = "Do Top ED50 Performers Have Different Observed Heterozygosity?",
    subtitle = paste("Top 20% (n =", n_top, ") vs random samples of same size"),
    x = "Group-Level Observed Heterozygosity (Ho)",
    y = "Frequency"
  ) +
  default_theme

ggsave("Genetic Analysis/outputs/top_performers_Ho.png", p_ho,
  width = 8, height = 6, dpi = 300
)

# =============================================================================
# 7. EXPECTED HETEROZYGOSITY (He)
# =============================================================================

cat("\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("EXPECTED HETEROZYGOSITY (He)\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

observed_he <- group_he(geno_matrix[top_rows, ])
overall_he <- group_he(geno_matrix[all_rows, ])

cat("\nGroup-level He:\n")
cat("  Top 20%:", round(observed_he, 4), "\n")
cat("  All individuals:", round(overall_he, 4), "\n")

null_he <- replicate(n_permutations, {
  idx <- sample(all_rows, n_top)
  group_he(geno_matrix[idx, ])
})

p_value_he <- mean(abs(null_he - overall_he) >= abs(observed_he - overall_he))

cat("\nPermutation test results (n =", n_permutations, "):\n")
cat("  Observed He (top 20%):", round(observed_he, 4), "\n")
cat("  Null distribution mean:", round(mean(null_he), 4), "\n")
cat("  Null distribution SD:", round(sd(null_he), 4), "\n")
cat("  P-value (two-tailed):", round(p_value_he, 4), "\n")

if (p_value_he < 0.05) {
  cat("\n*** SIGNIFICANT: Top performers have different He than expected ***\n")
} else {
  cat("\n*** NOT SIGNIFICANT: Top performers' He is not different from random ***\n")
}

he_null_df <- data.frame(he = null_he)

p_he <- ggplot(he_null_df, aes(x = he)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7, color = "white") +
  geom_vline(xintercept = observed_he, color = "red", linewidth = 1.2, linetype = "dashed") +
  annotate("text",
    x = observed_he, y = Inf,
    label = paste("Observed\np =", round(p_value_he, 3)),
    vjust = 2, hjust = -0.1, color = "red", size = 4
  ) +
  labs(
    title = "Do Top ED50 Performers Have Different Expected Heterozygosity?",
    subtitle = paste("Top 20% (n =", n_top, ") vs random samples of same size"),
    x = "Group-Level Expected Heterozygosity (He)",
    y = "Frequency"
  ) +
  default_theme

ggsave("Genetic Analysis/outputs/top_performers_He.png", p_he,
  width = 8, height = 6, dpi = 300
)

# =============================================================================
# 8. SENSITIVITY ANALYSIS: Different ED50 thresholds
# =============================================================================

cat("\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("SENSITIVITY ANALYSIS: Different ED50 thresholds\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

thresholds <- c(0.50, 0.80, 0.90)
sensitivity_results <- data.frame(
  threshold = character(),
  n_individuals = integer(),
  ho = numeric(),
  p_ho = numeric(),
  he = numeric(),
  p_he = numeric(),
  stringsAsFactors = FALSE
)

for (thresh in thresholds) {
  cutoff <- quantile(combined_data$ed50, probs = thresh)
  top_group <- combined_data %>% filter(ed50 >= cutoff)
  top_idx <- top_group$geno_row

  if (length(top_idx) >= 2) {
    obs_ho <- group_ho(geno_matrix[top_idx, ])
    null_ho_s <- replicate(1000, group_ho(geno_matrix[sample(all_rows, length(top_idx)), ]))
    p_ho_s <- mean(abs(null_ho_s - overall_ho) >= abs(obs_ho - overall_ho))

    obs_he <- group_he(geno_matrix[top_idx, ])
    null_he_s <- replicate(1000, group_he(geno_matrix[sample(all_rows, length(top_idx)), ]))
    p_he_s <- mean(abs(null_he_s - overall_he) >= abs(obs_he - overall_he))

    sensitivity_results <- rbind(sensitivity_results, data.frame(
      threshold = paste0("Top ", (1 - thresh) * 100, "%"),
      n_individuals = length(top_idx),
      ho = round(obs_ho, 4),
      p_ho = round(p_ho_s, 4),
      he = round(obs_he, 4),
      p_he = round(p_he_s, 4)
    ))
  }
}

print(sensitivity_results)

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("ANALYSIS COMPLETE\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("\nKey outputs:\n")
cat("1. Observed heterozygosity (Ho): proportion of het genotypes in the group, averaged across loci\n")
cat("2. Expected heterozygosity (He): 1 - sum(p^2) from allele frequencies, averaged across loci\n")
cat("3. Sensitivity analysis across different ED50 thresholds\n")
cat("\nFigures saved to Genetic Analysis/outputs/\n")
cat("  - top_performers_Ho.png\n")
cat("  - top_performers_He.png\n")
