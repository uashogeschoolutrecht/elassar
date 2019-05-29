#' # Elisa data analysis - concept
#' # Marc Teunis, Annemarie Stam
#' # Date July 2018
#'
#' Date of last change `r Sys.time()`
#'
#' ## Packages
library(lmtest)
library(sandwich)
library(drc)
library(tidyverse)

#' ## Read data
#' Reads the data for the calibration curve
calibration_df <- readxl::read_excel(
  path = here::here(
    "./data-raw/170217 - TNF ELISA.xlsx"),
    sheet = "calibration_data")
calibration_df

#' Reads the data for the samples (absorbance)
samples_df <- readxl::read_excel(
  path = here::here(
    "./data-raw/170217 - TNF ELISA.xlsx"),
    sheet = "sample_data")
samples_df

#' ## Calibration curves
#' Quick calibration curve
x <- calibration_df$conc
y <- calibration_df$absorbance
plot(y ~ x)

#' Get absorbance range for samples
sample_absorbance_range <- c(min(samples_df$absorbance),
                               max(samples_df$absorbance))
sample_absorbance_range

#' Calibration curve with ggplot2
  calibration_df %>%
    ggplot(aes(x = conc, y = absorbance)) +
      geom_point() +
      geom_smooth(data = calibration_df %>%
                group_by(conc) %>%
                summarise(mean_abs = mean(absorbance)),
              aes(x = conc, y = mean_abs)) +
      scale_x_log10() +
      ylab("Absorbance") +
      xlab("Standard concentration (pg/ml)") +
    geom_hline(yintercept = sample_absorbance_range[1],
               color = "red", linetype = "dashed") +
    geom_hline(yintercept = sample_absorbance_range[2],
               color = "red", linetype = "dashed") +
    theme_bw()


#' ## Model predictions using drc package
#' Define the model with a design formula
model <- drc::drm(absorbance ~ conc,
                  fct = drc:::LL.4(),
                  data = calibration_df)
summary(model)
drc:::plot.drc(model)

#' ## Estimating the sample concentrations via the models
predicted <- ED(model,
               respLev = samples_df$absorbance,
               type = 'absolute') %>%
    as.data.frame()

#' ## Add predicted conc to samples data
  samples_df <- samples_df %>%
    mutate(conc = predicted$Estimate) %>%
    mutate(dilution = fct_recode(as_factor(dilution),
          "1" = "undiluted",
          "5" = "5",
          "10" = "10")) %>%
    print()

#' ## ggplot2 graph for final results
#' Here we plot the calibration curve and the samples
#' to see a graph of all results. The red open circles are
#' the samples, the grey closed circles represent the
#' calibration curve. The purple vertical lines represent the
#' part of the curve that coincides with the absorbance.
#' The grey area around the curve represents the 95% confidence
#' interval of the fitted model.
#'
#' Some samples fall outside this
#' range, so probable the experiment needs repeating with
#' different dilutions of the samples or the amount of
#' measured cytokine is below the detection limit of the
#' ELISA assay.
ggplot(data = calibration_df, aes(x = conc, y = absorbance)) +
  geom_point(size = 4, alpha = 0.6) +
  geom_point(data = samples_df, aes(x = conc,
                                    y = absorbance),
             colour = "red",
             shape = 21,
             size = 5,
             alpha = 0.6) +
  geom_smooth(data = calibration_df %>%
                    group_by(conc) %>%
                    summarise(
                      mean_absorbance = mean(absorbance)
                      ),
                  aes(x = conc, y = mean_absorbance)) +
  scale_x_sqrt() +
  geom_vline(xintercept = max(samples_df$conc, na.rm = TRUE),
             linetype = "dashed",
                size = 1.5, colour = "purple") +
  geom_vline(xintercept = min(samples_df$conc, na.rm = TRUE),
             linetype = "dashed",
             size = 1.5, colour = "purple") +
  theme_bw()

#' ## Exploratory data analysis
calibration_df

samples_df

legend_ord <- samples_df %>%
  dplyr::select(treatment, conc) %>%
  na.omit() %>%
  group_by(treatment) %>%
  summarise(conc = mean(conc)) %>%
  mutate(legend_order = reorder(as_factor(treatment),
                                conc)) %>%
  arrange(conc)

legend <- legend_ord$legend_order

## preprocessing the data for graph -
## summary per experimental group
plot_groups <- samples_df %>%
  dplyr::select(treatment, conc) %>%
  na.omit() %>%
  group_by(treatment) %>%
  summarise(mean_conc = mean(conc)) %>%

## start plot call
  ggplot(aes(x = reorder(as_factor(treatment), mean_conc),
             y = mean_conc)) +
  geom_bar(aes(fill = treatment, color = treatment),
          position = 'dodge', stat = 'identity') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_color_discrete(breaks=legend) +
  scale_fill_discrete(breaks=legend) +
  ylab("Mean concentration") +
  xlab("Treatment")


plot_groups




