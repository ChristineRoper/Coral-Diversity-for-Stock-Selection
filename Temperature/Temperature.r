library(tibble) # for tibble
library(ggplot2)
library(readr) # for read_csv
library(dplyr)
source("default_theme.r")


sst_data <- read_csv("Temperature/outputs/sst_data.csv")
# reorder
sst_data$site <- factor(sst_data$site, levels = c(
    "Low Isles Reef",
    "Woody Isles Mangrove",
    "Tongue Reef",
    "Opal Reef Long Bommie",
    "Opal Reef Mojo",
    "Opal Reef Rayban"
))
sst_data$date_end <- as.POSIXct(sst_data$date_end)

# Add Opal Reef as the mean of the three Opal Reef sites
opal_reef <- sst_data %>%
    filter(site %in% c("Opal Reef Long Bommie", "Opal Reef Mojo", "Opal Reef Rayban")) %>%
    group_by(date_end) %>%
    summarise(sst = mean(sst, na.rm = TRUE), .groups = "drop") %>%
    mutate(site = "Opal Reef")

sst_data <- bind_rows(filter(sst_data, site %in% c("Low Isles Reef", "Tongue Reef")), opal_reef)
sst_data$site <- factor(sst_data$site, levels = c(
    "Low Isles Reef",
    "Tongue Reef",
    "Opal Reef"
))

# Add monthly means
sst_data <- sst_data %>%
    mutate(month = format(date_end, "%m")) %>%
    group_by(site, month) %>%
    mutate(overall_mean = mean(sst, na.rm = TRUE)) %>%
    ungroup()


plot <- ggplot(
    sst_data,
    aes(x = date_end, y = sst, color = site)
) +
    # geom_line(aes(alpha = ifelse(logger_data_exists, 0.1, 1))) +
    geom_line() +
    scale_alpha_continuous(range = c(0.7, 0.3)) +
    # Satellite mean
    geom_line(aes(y = overall_mean), linetype = "dashed") +
    #
    ylab("Temperature (ºC)") +
    xlab("Date") +
    xlim(as.POSIXct("2022-02-01 00:00:00"), as.POSIXct("2023-01-31 00:00:00")) +
    # xlim(as.POSIXct("2022-01-01 00:00:00"), as.POSIXct("2022-12-31 00:00:00")) +
    # scale_color_manual(values = c("#F84848", "#F6C64C", "#5ebb7f", "#82c4d9", "#9AB7EF", "#AB9AEF")) +
    scale_color_manual(values = c("#F84848", "#5ebb7f", "#AB9AEF")) +
    guides(
        alpha = "none", # Hide the legend for line opacity
        color = guide_legend(title = "Site")
    ) +
    default_theme

print(plot)

ggsave(
    "Temperature/plots/plotMM.png",
    plot,
    width = 8, height = 6
)
ggsave(
    "Temperature/plots/plot2Wrap.png",
    plot + facet_wrap(vars(site), nrow = 4) + theme(legend.position = "none"),
    width = 7, height = 8
)
