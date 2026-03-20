# Reads all Fluorcam files (T1, T2, T3) and creates a data frame with one row per measurement
# Input: Sample List.csv (SampleID, Folder, FirstArea, LastArea, Temp1, Temp2, Temp3, Temp4, Temp5, Temp6, ...)
# If the column "Folder" is blank, that sample is skipped
# columns from 11 onwards are not used but will be included in the output

flourcam_files <- "Phenotype Analysis/1. Fluorcam data/Fluorcam Files" # no trailing slash
sample_list <- "Phenotype Analysis/1. Fluorcam data/Sample list.csv"
output_file <- "Phenotype Analysis/1. Fluorcam data/outputs/fluorcam_data.csv"

# Metrics of interest from each timepoint file
t1_metrics <- c("Fm")
t2_metrics <- c(
    paste0("PAR", seq(1:15)), # Equivelent to writing PAR1, PAR2, ... PAR15
    paste0("ETR", seq(1:15)),
    paste0("QY_Lss", seq(1:15)),
    paste0("NPQ_Lss", seq(1:15)),
    paste0("Fq_Lss", seq(1:15)),
    paste0("Fm_Lss", seq(1:15))
)
t3_metrics <- c("Fq_Lss1", "Fm_Lss1", "Fv/Fm_Lss1", "Fm", "QY_max")

library(dplyr) # for filter
library(progress)

sample_list <- read.csv(sample_list, header = TRUE)
sample_list$SampleID <- as.factor(sample_list$SampleID)
sample_list$Folder <- as.factor(sample_list$Folder)

# The output data starts blank and we will add each row to it
output_data <- NULL


progress <- progress_bar$new(total = length(levels(sample_list$Folder)), format = "[:bar] :percent eta: :eta")

# go through each folder that data is stored in
for (folder in levels(sample_list$Folder)) {
    progress$tick()
    # Skip samples with no folder
    if (folder == "") next

    # Load the TXT files
    t1 <- read.table(paste(flourcam_files, folder, "T1.TXT", sep = "/"), header = TRUE, sep = "\t", skip = 2)
    t2 <- read.table(paste(flourcam_files, folder, "T2.TXT", sep = "/"), header = TRUE, sep = "\t", skip = 2)
    t3 <- read.table(paste(flourcam_files, folder, "T3.TXT", sep = "/"), header = TRUE, sep = "\t", skip = 2)

    # get the samples in this folder
    samples <- filter(sample_list, Folder == folder)

    # for each sample
    for (row_number in seq_len(nrow(samples))) {
        # This sample from sample_list
        sample <- samples[row_number, ]

        # Create a row for this sample
        new_row <- data.frame(SampleID = sample$SampleID)

        for (area_number in sample$FirstArea:sample$LastArea) {
            # columns are named "Area <area_number>"
            column_name <- paste0("Area.", area_number)

            new_row$Area <- area_number - sample$FirstArea + 1

            # temperature columns are named "Temp<temp_name>"
            temp_name <- paste0("Temp", new_row$Area)
            new_row$Temperature <- sample[, temp_name]

            # create a function to add metrics to 'data' that can be used multiple times
            add_metrics <- function(metrics, file, timepoint) {
                # for each metric that we're interested in
                for (metric in metrics) {
                    # Get the row for this metric from the TXT file
                    metric_data <- filter(file, across(1) == metric)

                    if (nrow(metric_data) == 0) {
                        stop(paste("There is no", metric, "in", folder, timepoint))
                    }

                    # set this metric's values in the new row
                    metric_name <- paste(timepoint, metric, sep = "_")
                    new_row[1, metric_name] <<- metric_data[, column_name]
                }
            }

            add_metrics(t1_metrics, t1, "T1")
            add_metrics(t2_metrics, t2, "T2")
            add_metrics(t3_metrics, t3, "T3")

            # Add the new row to 'data'
            output_data <- rbind(output_data, new_row)
        }
    }
}

write.csv(output_data, output_file)
