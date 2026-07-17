library(dplyr)
library(stats)

print("Loading prepared GII data for indicator selection analysis...")

# Load standardized data
data_standardized <- read.csv("data/gii_data_standardized.csv", stringsAsFactors = FALSE)
print(paste("Loaded data with", nrow(data_standardized), "countries and", ncol(data_standardized)-1, "indicators"))

# Prepare data for analysis (exclude Country column)
country_names <- data_standardized$Country
indicator_data <- data_standardized[, -1]  # Remove Country column

# Load the metadata to identify pillars
imeta <- read.csv("gii-data/imeta.csv", stringsAsFactors = FALSE)

# Extract indicators (rows where Level == 1)
indicators <- imeta[imeta$Level == 1, ]

# Ensure that all indicators have a positive direction
rev_indicators <- imeta$iCode[imeta$Direction == -1]
indicator_data <- indicator_data |> 
  mutate(across(all_of(rev_indicators), ~ max(.) - .))

# Remove all aggregates from the dataset - pillars, sub-index etc.
indicator_data <- indicator_data |> 
  select(all_of(indicators$iCode))

variable_names <- colnames(indicator_data)
print(paste("Final dataset:", nrow(indicator_data), "countries,", ncol(indicator_data), "indicators"))


select_indicators_by_pca(indicator_data)
