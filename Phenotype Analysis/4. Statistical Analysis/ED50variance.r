library(dplyr)

Fq_Fm_data <- read.csv("Phenotype Analysis/3. Thermal Curves/outputs/FqFm_T3/ed50.csv")

# make Group a factor
Fq_Fm_data$Group <- as.factor(Fq_Fm_data$Group)


# Check for Normality
#####################

# density plot to view distribution
plot(density(Fq_Fm_data$ed50))
shapiro.test(Fq_Fm_data$ed50) # did not pass


# leveneTest
############
# Non-parametric test for equal variances
library(car) # for leveneTest
leveneTest(ed50 ~ Group, Fq_Fm_data)
# Null Hypothesis: All populations variances are equal


# Manual Levene
###############

# Homebaked leveneTest by running an anova on distance to median

# Calculate means and add them to data
group_medians <- aggregate(ed50 ~ Group, data = Fq_Fm_data, FUN = median)
colnames(group_medians)[2] <- "group_median"
leveneAovData <- merge(Fq_Fm_data, group_medians, by = "Group")
# Add residuals to data: abs(value - median)
leveneAovData$residuals <- abs(leveneAovData$ed50 - leveneAovData$group_median)

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
write.csv(tukey_result$Group, "Phenotype Analysis/4. Statistical Analysis/ED50_variance_results.csv")
