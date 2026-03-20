# Top Performers Group-Level Allelic Richness Analysis
# Question: Do top ED50 performers have higher group-level allelic richness than expected?
# Group-level AR: rarefied allelic richness per locus, averaged across loci
# Rarefaction standardises for group size so top-20% vs random groups are comparable

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
# 3. HELPER FUNCTION FOR GROUP-LEVEL ALLELIC RICHNESS (rarefied)
# =============================================================================

# Rarefied allelic richness for a group of individuals.
#
# For each locus we count ref and alt allele copies among non-missing genotypes:
#   ref count = 2*(n_hom_ref) + n_het
#   alt count = 2*(n_hom_alt) + n_het
# where n = number of diploid individuals with data at that locus (2n alleles total).
#
# Rarefied AR at a locus is the expected number of distinct alleles (max 2) in a
# random draw of 2*g alleles, where g = rarefaction_n (the minimum number of
# diploid individuals with data across all groups being compared).
#
# E[AR] = sum over alleles a of: 1 - C(N - n_a, 2g) / C(N, 2g)
#   N  = total alleles at locus (2 * non-missing individuals)
#   n_a = count of allele a
#   g  = rarefaction sample size in diploid individuals
#
# Returns mean rarefied AR across loci (loci with N < 2*g are excluded).

rarefy_ar_locus <- function(n_ref, n_alt, g) {
  N <- n_ref + n_alt
  if (N < 2 * g) return(NA_real_)

  # log-scale combinatorial helper: log C(n, k)
  log_choose_safe <- function(n, k) {
    if (k < 0 || k > n) return(-Inf)
    lchoose(n, k)
  }

  log_denom <- log_choose_safe(N, 2 * g)

  # Probability that allele a is ABSENT in the draw = C(N - n_a, 2g) / C(N, 2g)
  p_absent_ref <- exp(log_choose_safe(N - n_ref, 2 * g) - log_denom)
  p_absent_alt <- exp(log_choose_safe(N - n_alt, 2 * g) - log_denom)

  # Expected number of distinct alleles = sum of P(allele present)
  (1 - p_absent_ref) + (1 - p_absent_alt)
}

group_ar <- function(geno_sub, g) {
  n_nonmissing <- colSums(!is.na(geno_sub))
  n_alt <- colSums(geno_sub, na.rm = TRUE) # 0->0, 1->1, 2->2 alt alleles
  n_ref <- 2 * n_nonmissing - n_alt

  locus_ar <- mapply(rarefy_ar_locus, n_ref, n_alt, MoreArgs = list(g = g))
  mean(locus_ar, na.rm = TRUE)
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
# 6. SET RAREFACTION SIZE
# =============================================================================
# Rarefy to the top-20% group size so that random samples drawn in permutations
# always match the observed group size -- comparisons are automatically fair.

rarefaction_n <- n_top
cat("\nRarefaction size (g):", rarefaction_n, "diploid individuals (", 2 * rarefaction_n, "alleles )\n")

# =============================================================================
# 7. ALLELIC RICHNESS (AR)
# =============================================================================

cat("\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("ALLELIC RICHNESS (AR, rarefied)\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

observed_ar <- group_ar(geno_matrix[top_rows, ], g = rarefaction_n)
overall_ar  <- group_ar(geno_matrix[all_rows, ],  g = rarefaction_n)

cat("\nGroup-level AR (rarefied to g =", rarefaction_n, "):\n")
cat("  Top 20%:", round(observed_ar, 4), "\n")
cat("  All individuals:", round(overall_ar, 4), "\n")

set.seed(42)
n_permutations <- 10000

null_ar <- replicate(n_permutations, {
  idx <- sample(all_rows, n_top)
  group_ar(geno_matrix[idx, ], g = rarefaction_n)
})

p_value_ar <- mean(abs(null_ar - overall_ar) >= abs(observed_ar - overall_ar))

cat("\nPermutation test results (n =", n_permutations, "):\n")
cat("  Observed AR (top 20%):", round(observed_ar, 4), "\n")
cat("  Null distribution mean:", round(mean(null_ar), 4), "\n")
cat("  Null distribution SD:", round(sd(null_ar), 4), "\n")
cat("  P-value (two-tailed):", round(p_value_ar, 4), "\n")

if (p_value_ar < 0.05) {
  cat("\n*** SIGNIFICANT: Top performers have different AR than expected ***\n")
} else {
  cat("\n*** NOT SIGNIFICANT: Top performers' AR is not different from random ***\n")
}

ar_null_df <- data.frame(ar = null_ar)

p_ar <- ggplot(ar_null_df, aes(x = ar)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7, color = "white") +
  geom_vline(xintercept = observed_ar, color = "red", linewidth = 1.2, linetype = "dashed") +
  annotate("text",
    x = observed_ar, y = Inf,
    label = paste("Observed\np =", round(p_value_ar, 3)),
    vjust = 2, hjust = -0.1, color = "red", size = 4
  ) +
  labs(
    title = "Do Top ED50 Performers Have Different Allelic Richness?",
    subtitle = paste0("Top 20% (n = ", n_top, ") vs random samples of same size, rarefied to g = ", rarefaction_n),
    x = "Group-Level Allelic Richness (rarefied AR)",
    y = "Frequency"
  ) +
  default_theme

ggsave("Genetic Analysis/outputs/top_performers_AR.png", p_ar,
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
  ar = numeric(),
  p_ar = numeric(),
  stringsAsFactors = FALSE
)

for (thresh in thresholds) {
  cutoff <- quantile(combined_data$ed50, probs = thresh)
  top_group <- combined_data %>% filter(ed50 >= cutoff)
  top_idx <- top_group$geno_row
  g_s <- length(top_idx)

  if (g_s >= 2) {
    obs_ar <- group_ar(geno_matrix[top_idx, ], g = g_s)
    overall_ar_s <- group_ar(geno_matrix[all_rows, ], g = g_s)
    null_ar_s <- replicate(1000, group_ar(geno_matrix[sample(all_rows, g_s), ], g = g_s))
    p_ar_s <- mean(abs(null_ar_s - overall_ar_s) >= abs(obs_ar - overall_ar_s))

    sensitivity_results <- rbind(sensitivity_results, data.frame(
      threshold = paste0("Top ", (1 - thresh) * 100, "%"),
      n_individuals = g_s,
      ar = round(obs_ar, 4),
      p_ar = round(p_ar_s, 4)
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
cat("1. Allelic richness (AR): expected number of distinct alleles per locus in a random\n")
cat("   draw of 2*g alleles (rarefied), averaged across loci. Max value = 2 (biallelic SNPs).\n")
cat("2. Sensitivity analysis across different ED50 thresholds\n")
cat("\nFigures saved to Genetic Analysis/outputs/\n")
cat("  - top_performers_AR.png\n")
