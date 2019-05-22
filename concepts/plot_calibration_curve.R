#' # Elisa data analysis - concept
#' # Marc Teunis, Annemarie Stam
#' # Date July 2018
#'
#' Date of last change `r Sys.time()`
#'
#' ## Packages
#'
#'+ packages
library(lmtest)
library(sandwich)
library(drc)
library(tidyverse)
library(tidyverse)

#' ## read data
#"+ data_load
calibration_df <- readxl::read_excel(path = "./data-raw/170217 - TNF ELISA.xlsx",
                                        sheet = "calibration_data")
calibration_df

samples_df <- readxl::read_excel(path = "./data-raw/170217 - TNF ELISA.xlsx",
                                    sheet = "sample_data")
samples_df

## quick calibration curve
x <- calibration_df$conc
y <- calibration_df$absorbance
plot(y ~ x)

## get absorbance range for samples
sample_absorbance_range <- c(min(samples_df$absorbance),
                               max(samples_df$absorbance))
sample_absorbance_range

## calibration curve with ggplot2
  calibration_df %>%
    ggplot(aes(x = conc, y = absorbance)) +
      geom_point() +
      geom_line(data = calibration_df %>%
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

## model predictions using drc package
model <- drc::drm(absorbance ~ conc,
                  fct = drc:::LL.4(),
                  data = calibration_df)
summary(model)
drc:::plot.drc(model)



## heterogeneity in variance; robust standard errors
#robust <- lmtest::coeftest(model,
#                    vcov = sandwich)
#robust

## estimates
#?ED
predicted <- ED(model,
               respLev = samples_df$absorbance,
               type = 'absolute') %>%
    as.data.frame()
  predicted

  # warnings()

## add predicted conc to samples data
  samples_df <- samples_df %>%
    mutate(conc = predicted$Estimate) %>%
    mutate(dilution = fct_recode(as_factor(dilution),
          "1" = "undiluted",
          "5" = "5",
          "10" = "10")) %>%
    mutate(conc_new = conc*(as.numeric(dilution))) %>%
    print()

  samples_df

## predictions using broom::augment
  levels(as_factor(samples_df$dilution))

#

#  ggplot model
model$curve
 # ggplot for final results



ggplot(data = calibration_df, aes(x = conc, y = absorbance)) +
  geom_point(size = 4, alpha = 0.6) +
  geom_point(data = samples_df, aes(x = conc,
                                    y = absorbance),
                colour = "red", shape = 21, size = 5, alpha = 0.6) +



        geom_line(data = calibration_df %>%
                    group_by(conc) %>%
                    summarise(mean_absorbance = mean(absorbance)),
                  aes(x = conc, y = mean_absorbance)) +
  # geom_line(data = samples_df,  aes(x = conc, y = absorbance ),
    #        colour = 'darkgreen', size = 3, alpha = 0.8, linetype = "dashed") +

     scale_x_sqrt() +
     geom_vline(xintercept = 0.31, linetype = "dashed",
                size = 1.5, colour = "purple") +

    geom_vline(xintercept = 5, linetype = "dashed",
             size = 1.5, colour = "purple") +
  theme_bw()
 plot

## exploratory data analysis

## concentrations per group
calibration_df

samples_df

legend_ord <- samples_df %>%
  dplyr::select(treatment, conc_new) %>%
  na.omit() %>%
  group_by(treatment) %>%
  summarise(mean_conc = mean(conc_new)) %>%
  mutate(legend_order = reorder(as_factor(treatment), mean_conc)) %>%
  arrange(mean_conc)

legend <- legend_ord$legend_order


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

 # geom_point(aes(color = treatment), size = 2) +
 # theme_bw()
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_color_discrete(breaks=legend) +
  scale_fill_discrete(breaks=legend) +
  ylab("Mean concentration") +
  xlab("Treatment")


plot_groups




