## create dataset "elisa"
library(tidyverse)

elisa_calibration <- readxl::read_excel(path = "./data-raw/170217 - TNF ELISA (version 1).xlsx",
                   sheet = "calibration_data")

devtools::use_data(elisa_calibration)




elisa_samples <- readxl::read_excel(path = "./data-raw/170217 - TNF ELISA (version 1).xlsx",
                                        sheet = "sample_data")

elisa_samples_tidy <- elisa_samples %>%
  gather(absorbance, X__1, key = "measured", value = "value") %>%
  print()


devtools::use_data(elisa_samples, overwrite = TRUE)
