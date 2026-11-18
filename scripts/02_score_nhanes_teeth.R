
# Dental scoring code is within data-raw/02_dental_scoring.R, which is sourced here.

stopifnot(requireNamespace("here", quietly = TRUE))
stopifnot(requireNamespace("dplyr", quietly = TRUE))
stopifnot(requireNamespace("tidyr", quietly = TRUE))
stopifnot(requireNamespace("stringr", quietly = TRUE))
source(here::here("data-raw", "02_dental_scoring.R"), local = TRUE)
message("Dental scores saved to data/icspack_Tooth_Level_Data.RData")
