# Processes the RLC part of the Fluorcam data with PI curves to extract Pmax, alpha, ek

input <- "Phenotype Analysis/2. Rapid Light Curves/outputs/fluorcam_RLC_data.csv"
output <- "Phenotype Analysis/2. Rapid Light Curves/outputs/FqFm/rlc_metrics.csv"
plots <- "Phenotype Analysis/2. Rapid Light Curves/outputs/FqFm/plots/"

library(readr)
library(ggplot2)
library(nls.multstart)
library(minpack.lm) # for nlsLM
library(broom) # for augment
library(progress)
source("default_theme.r")

rlc_data <- read_csv(input)
rlc_data$SampleID <- as.factor(rlc_data$SampleID)
rlc_data$Area <- as.factor(rlc_data$Area)


eilers_pi <- function(pmax, iopt, a, i) {
  pi <- (pmax * i) / ((pmax / (a * iopt^2)) * i^2 + ((1 - (2 * pmax) / (a * iopt)) * i) + (pmax / a))
  return(pi)
}

eilers_FqFm_pi <- function(FqFm_max, IoptExtra, Ek, I) {
  Iopt <- Ek + IoptExtra
  Pmax <- FqFm_max * Ek
  alpha <- Pmax / Ek
  return(eilers_pi(Pmax, Iopt, alpha, I) / I)
}

sample_ids <- levels(as.factor(rlc_data$SampleID)) # all samples

output_data <- NULL

progress <- progress_bar$new(total = length(sample_ids), format = ":sample [:bar] :percent eta: :eta")
for (sample in sample_ids) {
  progress$tick(tokens = list(sample = sample))

  plot <- ggplot() + default_theme
  fit_curves <- NULL

  sample_data <- dplyr::filter(rlc_data, SampleID == sample)

  # Ensure no FqFm is negative
  sample_data$FqFm <- pmax(sample_data$FqFm, 0)

  for (area in 1:6) {
    area_data <- dplyr::filter(sample_data, Area == area)

    simple_max <- max(area_data$FqFm)

    # Fit
    ###########

    fit <- tryCatch(
      {
        fit <- nls_multstart(
          FqFm ~ eilers_FqFm_pi(FqFm_max, IoptExtra, Ek, PAR),
          data = area_data,
          iter = c(1, 3, 5),
          start_lower = list(FqFm_max = simple_max, IoptExtra = 0, Ek = 20),
          start_upper = list(FqFm_max = simple_max, IoptExtra = 200, Ek = 1000),
          lower = c(0, 0, 0),
          upper = c(0.8, 1700, 1700),
          supp_errors = "Y",
          convergence_count = FALSE
        )
      },
      error = function(e) {
        cat("error in ", sample, " area ", area, " fit 1\n")
        stop(e)
      }
    )


    Ek <- coef(fit)[["Ek"]]
    FqFm_max <- coef(fit)[["FqFm_max"]]
    IoptExtra <- coef(fit)[["IoptExtra"]]
    Iopt <- Ek + IoptExtra

    new_row <- data.frame(
      FqFm_max = FqFm_max,
      Ek = Ek,
      FqFm_ek = eilers_FqFm_pi(FqFm_max, IoptExtra, Ek, Ek),
      Iopt = Iopt,
      AIC = AIC(fit),
      SampleID = sample,
      Area = area
    )

    output_data <- rbind(output_data, new_row)

    plot <- plot +
      geom_segment(
        data = data.frame(FqFm_max = FqFm_max, Ek = Ek, IoptExtra = IoptExtra),
        mapping = aes(x = Ek, y = 0, xend = Ek, yend = eilers_FqFm_pi(FqFm_max, IoptExtra, Ek, Ek)),
        color = "grey", alpha = 0.5
      ) +
      geom_point(
        data = data.frame(FqFm_max = FqFm_max, Ek = Ek, IoptExtra = IoptExtra),
        mapping = aes(x = Ek, y = eilers_FqFm_pi(FqFm_max, IoptExtra, Ek, Ek)),
        pch = 4, size = 2.5
      ) +
      geom_segment(
        data = data.frame(FqFm_max = FqFm_max, Ek = Ek, IoptExtra = IoptExtra),
        mapping = aes(x = 0, y = FqFm_max, xend = Ek, yend = 0),
        color = "grey", linetype = "dotted", alpha = 0.4
      )

    # Save the curve for plotting later
    pars <- seq(0, 1800, 1) # get pars from 0 to 1800 to plot the curve against
    curve <- data.frame(PAR = pars, FqFm = eilers_FqFm_pi(FqFm_max, IoptExtra, Ek, pars), Area = area)
    fit_curves <- rbind(fit_curves, curve)
  }

  fit_curves$Area <- as.factor(fit_curves$Area)


  # plot
  plot <- plot +
    geom_hline(yintercept = 0, color = "grey", alpha = 0.6) +
    # Data
    geom_line(data = sample_data, mapping = aes(x = PAR, y = FqFm, col = Area), alpha = 0.2) +
    geom_point(data = sample_data, mapping = aes(x = PAR, y = FqFm, col = Area)) +
    # Fit 1
    # geom_line(data = sample_data, mapping = aes(x = PAR, y = FqFm_fit, col = Area), linetype = 1) +
    geom_line(data = fit_curves, mapping = aes(x = PAR, y = FqFm, col = Area), linetype = 1) +

    coord_cartesian(ylim = c(0, 1), xlim = c(0, max(sample_data$PAR))) +
    labs(x = "PAR", y = "FqFm") +
    ggtitle(sample)

  suppressMessages(ggsave(paste0(plots, sample, ".png"), plot))
}

write.csv(output_data, output)
