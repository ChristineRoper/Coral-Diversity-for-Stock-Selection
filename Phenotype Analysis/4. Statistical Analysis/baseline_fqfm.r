input <- "Phenotype Analysis/2. Rapid Light Curves/outputs/FqFm/rlc_metrics.csv"
sample_list <- "Phenotype Analysis/1. Fluorcam data/Sample list.csv"


library(dplyr)

sample_list <- read_csv(sample_list) %>%
  # Get SampleID and Group from sample_list
  dplyr::select(c(SampleID, Group))

baseline_data <- read_csv(input) %>%
  filter(
    # Filter leaving only datapoints at baseline temp
    Area == 1,
    # filter mangrove sample as no replication in that group
    SampleID != "M1",
    # remove W5
    SampleID != "W5"
  ) %>%
  left_join(sample_list, by = "SampleID") %>%
  mutate(Group = as.factor(Group))

write.csv(baseline_data, "Phenotype Analysis/4. Statistical Analysis/outputs/Baseline_FqFm.csv")

aggregated <- baseline_data %>%
  group_by(Group) %>%
  summarise(
    Mean = mean(FqFm_max),
    StdErr = sd(FqFm_max) / sqrt(n()),
    Min = min(FqFm_max),
    Max = max(FqFm_max),
    Range = max(FqFm_max) - min(FqFm_max)
  )

write.csv(aggregated, "Phenotype Analysis/4. Statistical Analysis/outputs/Baseline_FqFm.csv")


## statistical tests ##

# Check for Normality
#####################

### ANOVA on mean fq/fm_max ####

# load required libraries
library(stats)
library(FSA)

# Check assumptions for ANOVA
shapiro.test(baseline_data$FqFm_max) # did not pass
bartlett.test(FqFm_max ~ Group, data = baseline_data) # did not pass

# Assumptions for ANOVA not met, using non-parametric alternative
# Kruskal-Wallis Rank Sum Test
kruskal.test(FqFm_max ~ Group, data = baseline_data)

# Post-hoc Dunn's test of multiple comparisons
# with Benjamini-Hochberg adjustment
dunn_result <- dunnTest(FqFm_max ~ Group, data = baseline_data, method = "bh")
View(dunn_result$res)
write.csv(dunn_result$res, "Phenotype Analysis/4. Statistical Analysis/outputs/fq_fm_baseline_posthoc_results.csv")



# leveneTest
############

# density plot to view distribution
plot(density(baseline_data$FqFm_max))
shapiro.test(baseline_data$FqFm_max) # did not pass

# Non-parametric test for equal variances
library(car) # for leveneTest
leveneTest(FqFm_max ~ Group, baseline_data)
# Null Hypothesis: All populations variances are equal


# Manual Levene
###############

# Homebaked leveneTest by running an anova on distance to median

# Calculate means and add them to data
group_medians <- aggregate(FqFm_max ~ Group, data = baseline_data, FUN = median)
colnames(group_medians)[2] <- "group_median"
leveneAovData <- merge(baseline_data, group_medians, by = "Group")
# Add residuals to data: abs(value - median)
leveneAovData$residuals <- abs(leveneAovData$FqFm_max - leveneAovData$group_median)

print(aggregate(residuals ~ Group, data = leveneAovData, FUN = mean))

plot <- leveneAovData %>% ggplot(aes(x = Group, y = residuals, fill=Group)) +
    geom_boxplot() +
    scale_fill_manual(values = c(
        "Low Isles slope" = "red",
        "Tongue Turtle Bay" = "lightblue",
        "Opal Rayban flat" = "lightblue",
        "Tongue Blue Hole" = "lightblue"
    ))
ggsave("box.png", plot)

# anova on residuals
levene <- aov(residuals ~ Group, data = leveneAovData)
summary(levene)

tukey_result <- TukeyHSD(levene)
print(tukey_result)
write.csv(tukey_result$Group, "Phenotype Analysis/4. Statistical Analysis/outputs/Fq_fm_baseline_variance_results.csv")
