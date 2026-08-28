# Data Preparation Script for GII Analysis
# Produces the processed data and descriptive-statistics tables that the
# report reads. Nothing downstream in the pipeline consumes these files -
# 02/03/03b/04 each call load_gii_data() themselves - so this script's job is
# to persist the report's inputs, not to feed the later analysis stages.

library(dplyr)
library(COINr)
source("scripts/R/config.R")
source("scripts/R/functions.R")
source("scripts/R/data_utils.R")

message("Loading GII data...")
gii_data <- load_gii_data()

idata <- gii_data$raw_data
complete_data <- gii_data$complete_data

message("Loaded data with ", nrow(complete_data), " countries and ",
        ncol(complete_data) - 1, " indicators")

write.csv(complete_data, CONFIG$paths$data_complete, row.names = FALSE)
write.csv(gii_data$standardized_data, CONFIG$paths$data_std, row.names = FALSE)

# Imputation summary quoted by the report's appendix
write.csv(gii_data$missingness, CONFIG$paths$imputation, row.names = FALSE)
message("Imputation report saved to: ", CONFIG$paths$imputation)

#' Descriptive statistics for every indicator column of `df`
describe_indicators <- function(df) {
  data.frame(
    Indicator = names(df),
    Mean = colMeans(df, na.rm = TRUE),
    Median = apply(df, 2, median, na.rm = TRUE),
    Min = apply(df, 2, min, na.rm = TRUE),
    Max = apply(df, 2, max, na.rm = TRUE),
    Kurtosis = apply(df, 2, function(x) COINr::kurt(x, na.rm = TRUE)),
    Skewness = apply(df, 2, function(x) COINr::skew(x, na.rm = TRUE))
  )
}

message("Generating descriptive statistics...")
write.csv(describe_indicators(idata[, -1]), CONFIG$paths$stats_original, row.names = FALSE)
write.csv(describe_indicators(complete_data[, -1]), CONFIG$paths$stats_complete, row.names = FALSE)

message("Descriptive statistics tables saved to /data directory.")
message("Data preparation completed.")
