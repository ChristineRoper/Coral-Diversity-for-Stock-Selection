# Process RLC part of fluorcam_data from extract_fluorcam_data.rmd into a tidy format with one row per datapoint
fluorcam_data <- "Phenotype Analysis/1. Fluorcam data/outputs/fluorcam_data.csv"
output <- "Phenotype Analysis/2. Rapid Light Curves/outputs/fluorcam_RLC_data.csv"

library(progress)

fluorcam_data <- read.csv(fluorcam_data)

# A variable to start adding outputs to
etr_data <- NULL

# Each row in fluorcam_data is a sample
progress <- progress_bar$new(total = nrow(fluorcam_data), format = "[:bar] :percent")
for (i in seq_len(nrow(fluorcam_data))) {
  progress$tick()

  sample <- fluorcam_data[i, ]

  # There are 15 columns for each of PAR, ETR and QY_Lss
  for (m in seq(1:15)) {
    # A new row for our output with sample info plus PAR, ETR and QY_Lss
    new_row <- data.frame(
      SampleID = sample$SampleID,
      Area = sample$Area,
      Temperature = sample$Temperature,
      PAR = fluorcam_data[i, paste0("T2_PAR", m)],
      ETR = fluorcam_data[i, paste0("T2_ETR", m)],
      QY = fluorcam_data[i, paste0("T2_QY_Lss", m)],
      FqFm = fluorcam_data[i, paste0("T2_Fq_Lss", m)] / fluorcam_data[i, paste0("T2_Fm_Lss", m)]
    )
    # Add the row to the output
    etr_data <- rbind(etr_data, new_row)
  }
}

write.csv(etr_data, output)
