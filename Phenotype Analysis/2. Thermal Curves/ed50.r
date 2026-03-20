metric <- "FqFm_T3"
fluorcam_data <- "Phenotype Analysis/1. Fluorcam data/outputs/fluorcam_data.csv"
sample_list <- "Phenotype Analysis/1. Fluorcam data/Sample list.csv"
output <- paste0("Phenotype Analysis/3. Thermal Curves/outputs/", metric, "/")
plots <- paste0("Phenotype Analysis/3. Thermal Curves/plots/", metric, "/")

metric_names <- list(
  FqFm_T3 = "Fq'/Fm'",
  FqFv_max = "Initial Fq'/Fm'"
)

library(dplyr)
library(readr)
library(drc)
library(broom)
library(progress)
source("default_theme.r")

data <- NULL

if (metric == "FqFm_T3") {
  fluorcam_data <- read.csv(fluorcam_data)

  for (i in seq_len(nrow(fluorcam_data))) {
    sample <- fluorcam_data[i, ]
    new_row <- data.frame(
      SampleID = sample$SampleID,
      Area = sample$Area,
      Temperature = sample$Temperature,
      Metric = fluorcam_data[i, "T3_Fq_Lss1"] / fluorcam_data[i, "T3_Fm_Lss1"]
    )
    data <- rbind(data, new_row)
  }
}

data$SampleID <- as.factor(data$SampleID)
data$Area <- as.factor(data$Area)

sample_list <- read_csv(sample_list) %>%
  # Drop some columns we don't need now
  dplyr::select(!c(Folder,FirstArea,LastArea,Temp1,Temp2,Temp3,Temp4,Temp5,Temp6))
data <- left_join(data, sample_list, by = "SampleID")

#################
# By Individual #
#################

# samples <- list("W7", "LI F5", "MO S3", "TB 2") # trial representative samples
samples <- levels(as.factor(data$SampleID)) # all samples
output_data <- NULL

message("ED50 by individual")
progress <- progress_bar$new(total = length(samples), format = ":sample [:bar] :percent")
for (sample_name in samples) {
  progress$tick(tokens = list(sample = sample_name))

  sample_data <- filter(data, SampleID == sample_name)

  # Remove well 6 if it's higher than well 5
  # Remove well 5 if it's higher than well 4
  if (!is.na(sample_data[sample_data$Area == 6, "Metric"])) {
    if (sample_data[sample_data$Area == 6, "Metric"] > sample_data[sample_data$Area == 5, "Metric"]) {
      sample_data[sample_data$Area == 6, "Metric"] <- NA
    } else if (sample_data[sample_data$Area == 6, "Metric"] > sample_data[sample_data$Area == 4, "Metric"]) {
      sample_data[sample_data$Area == 6, "Metric"] <- NA
    }
  }
  if (sample_data[sample_data$Area == 5, "Metric"] > sample_data[sample_data$Area == 4, "Metric"]) {
    sample_data[sample_data$Area == 5, "Metric"] <- NA
  }

  # plot data and model fit
  plot <- ggplot(sample_data, aes(Temperature, Metric)) +
    geom_point() +
    default_theme +
    labs(
      x = "Temperature (ºC)", y = metric_names[[metric]],
      title = sample_name
    )

  #######
  # Fit #
  #######
  drc <- tryCatch(
    {
      drm(Metric ~ Temperature,
        data = sample_data,
        fct = LL.3(),
        # limits for Slope, Max and ED50
        upperl = c(300,1,40),
        lowerl = c(0,0,30)
      )
    },
    error = function(cond) {
      message('\n', sample_name)
      message(conditionMessage(cond))
      NULL
    }
  )

  # skip the rest if we couldn't fit
  if (is.null(drc)) {
    suppressMessages(ggsave(paste0(plots, "failed/", sample_name, ".png"), plot))
    next # go on to the next sample
  }

  #################
  # Output Params #
  #################

  slope <- coef(drc)["b:(Intercept)"]
  max <- coef(drc)["d:(Intercept)"]
  ed50 <- coef(drc)["e:(Intercept)"]

  # Identify and notify if the sample's ED50 is below the lowest value
  # In this case, don't include it in output
  if(max / 2 < min(sample_data$Metric, na.rm = TRUE)){
    message('\nED50 lower than last point: ', sample_name)
    next # Move on to the next sample
  }

  output_row <- data.frame(
    SampleID = sample_name,
    max = max,
    ed50 = ed50,
    slope = slope,
    ip = ed50 - max / slope / 2
  )

  output_data <- rbind(output_data, output_row)

  ########
  # Plot #
  ########

  new_data <- data.frame(Temperature = seq(20, 40, 0.2))
  new_data$predicted <- predict(drc, newdata = new_data)


  # plot data and model fit
  plot <- plot +
    # ED50
    geom_segment(aes(x = ed50, xend = ed50, y = -1, yend = max / 2), output_row, linetype = "dashed") +
    # ip
    geom_segment(aes(x = ip, xend = ip, y = -1, yend = max), output_row, linetype = "dashed") +
    # Max
    geom_hline(aes(yintercept = max), output_row, linetype = "dashed") +
    # Slope
    geom_abline(aes(slope = -slope, intercept = slope * ed50 + max / 2), output_row, linetype = "dotted") +
    # Curve
    geom_line(aes(Temperature, predicted), new_data, col = "blue") +
    # ED50 point
    geom_point(aes(ed50, max / 2), output_row, color = "blue", size = 4, pch = 4) +
    coord_cartesian(
      xlim = c(24, 38),
      ylim = c(0, max(sample_data[["Metric"]]))
    )

  # Start a new PNG (much faster than ggsave)
  png(
    filename = paste0(plots, "individual/", sample_name, ".png"),
    width = 1600, height = 1200, units = "px", res = 300
  )
  print(plot) # Create and print your plot
  dev.off() # Output and close the png file
}

# Add in metadata such as Group from sample_list
output_data = output_data %>%
  left_join(sample_list, by = "SampleID")

write.csv(output_data, paste0(output, "ed50.csv"))


# Write mean and SE summary CSV
output_data %>%
  group_by(Group) %>%
  summarise(
    ED50 = mean(ed50),
    stderr = sd(ed50) / sqrt(n()),
    Variance = var(ed50),
    Max = max(ed50),
    Min = min(ed50),
    Range = max(ed50) - min(ed50)
  ) %>%
  write.csv(paste0(output, "ed50_summary.csv"))


#####################
# Log Logistic Plot #
#####################

means_ed50 <- output_data %>%
  group_by(Group) %>%
  summarise(ed50 = mean(ed50), slope = mean(slope), max = mean(max))

mean_points = data %>%
  group_by(Group, Temperature) %>%
  summarise(Metric = mean(Metric))

groups = levels(as.factor(means_ed50$Group))
temperatures <- seq(23, 42, 0.05) # making a sequence from 30 to 42 in steps of 0.5
data.frame(
  Group = rep(groups, each = length(temperatures)), # repeat the group for the sequence of temperatures
  Temperature = rep(temperatures, times = length(groups)) #
) %>%
  left_join(means_ed50, "Group") %>%
  mutate(Metric = max / (1 + (Temperature / ed50)^slope)) %>%
  ggplot(aes(x = Temperature, y = Metric, color = Group)) +
  geom_line() +
  # geom_errorbar(data = means, aes(y = mean, ymin = mean - stderr, ymax = mean + stderr, color = Group), width = 0.1) +
  geom_point(data = mean_points) +
  # geom_point(aes(x = ed50, y = max * 0.50)) +
  # geom_point(aes(x = ed95, y = max * 0.05)) +
  # geom_point(aes(x = ed5, y = max * 0.95)) +
  coord_cartesian(xlim = c(24, 40)) +
  labs(x = "Temperature (ºC)", y = metric_names[[metric]]) +
  default_theme

ggsave(paste0(plots, "ed-50-curves.png"), width = 2000, height = 2000, units = "px")

