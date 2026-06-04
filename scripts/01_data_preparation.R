# Data Preparation Script for GII Analysis

library(readxl)
library(dplyr)
library(e1071) # For kurtosis and skewness functions

# Load GII data
print("Loading GII data...")

# Load data files
idata <- read_excel("gii-data/idata.xlsx")
imeta <- read_excel("gii-data/imeta.xlsx")
print("Data loaded successfully")

# Examine the structure
print(paste("Original data dimensions:", dim(idata)[1], "rows,", dim(idata)[2], "columns"))

# Identify indicator columns that are at Level 1 from imeta
level1_indicators <- imeta$iCode[imeta$Level == 1]
indicator_columns <- level1_indicators

# Subset the data to only include Level 1 indicators
idata_level1 <- idata[, c("uCode", indicator_columns)]

# Check missing data pattern for Level 1 indicators
missing_counts <- colSums(is.na(idata_level1))
total_missing_before <- sum(missing_counts)
num_vars_missing <- sum(missing_counts > 0)
num_countries_missing <- sum(rowSums(is.na(idata_level1[, -1])) > 0)

print("Missing data summary for Level 1 indicators:")
print(paste("Total missing values:", total_missing_before))
print(paste("Columns with missing values:", num_vars_missing))
print(paste("Countries with missing data:", num_countries_missing))

# Install and load VIM package for kNN imputation if not already installed
if (!require("VIM", quietly = TRUE)) {
  install.packages("VIM", repos = "https://cloud.r-project.org/")
  library(VIM)
} else {
  library(VIM)
}

# Perform kNN imputation instead of case-wise deletion
print("Performing kNN imputation for missing data...")
imputed_data <- kNN(idata_level1, 
  variable = names(idata_level1)[-1],
   k = 5, imp_var = FALSE )

# Check missing data after imputation
missing_counts_after <- colSums(is.na(imputed_data))
total_missing_after <- sum(missing_counts_after)
num_countries_after <- sum(rowSums(is.na(imputed_data[, -1])) > 0)

print("Missing data summary after imputation:")
print(paste("Total missing values:", total_missing_after))
print(paste("Countries with missing data:", num_countries_after))

# Remove rows that still have missing values after imputation (if any)
complete_cases <- complete.cases(imputed_data[, -1])
complete_data <- imputed_data[complete_cases, ]

final_missing <- sum(is.na(complete_data))
print(paste("Countries with complete data after imputation:", nrow(complete_data)))
print(paste("Missing values remaining:", final_missing))

# Create imputation report for documentation
imputation_report <- data.frame(
  Metric = c("Total Missing Values", 
             "Variables with Missing Data", 
             "Countries with Missing Data"),
  Count = c(total_missing_before, num_vars_missing,
            num_countries_missing)
)

write.csv(imputation_report, "data/imputation_report.csv", row.names = FALSE)

# Standardize data (z-scores) for indicators only
print("Standardizing data...")
indicator_data <- complete_data[, -1]  # Exclude uCode column
standardized_data <- as.data.frame(scale(indicator_data))
standardized_data_with_country <- cbind(Country = complete_data$uCode, standardized_data)

# Save cleaned and standardized data
write.csv(complete_data, "data/gii_data_complete_cases.csv", row.names = FALSE)
write.csv(standardized_data_with_country, "data/gii_data_standardized.csv", row.names = FALSE)
write.csv(imputation_report, "data/imputation_report.csv", row.names = FALSE)

# Generate descriptive statistics for original indicator data
print("Generating descriptive statistics for original data...")
original_indicators <- idata_level1[, -1]  # Exclude uCode column
original_stats <- data.frame(
  Indicator = names(original_indicators),
  Mean = colMeans(original_indicators, na.rm = TRUE),
  Median = apply(original_indicators, 2, median, na.rm = TRUE),
  Min = apply(original_indicators, 2, min, na.rm = TRUE),
  Max = apply(original_indicators, 2, max, na.rm = TRUE),
  Kurtosis = apply(original_indicators, 2, function(x) kurtosis(x, na.rm = TRUE)),
  Skewness = apply(original_indicators, 2, function(x) skewness(x, na.rm = TRUE))
)

# Generate descriptive statistics for complete data after imputation and cleaning
print("Generating descriptive statistics for complete data...")
complete_indicators <- complete_data[, -1]  # Exclude uCode column
complete_stats <- data.frame(
  Indicator = names(complete_indicators),
  Mean = colMeans(complete_indicators),
  Median = apply(complete_indicators, 2, median),
  Min = apply(complete_indicators, 2, min),
  Max = apply(complete_indicators, 2, max),
  Kurtosis = apply(complete_indicators, 2, function(x) kurtosis(x)),
  Skewness = apply(complete_indicators, 2, function(x) skewness(x))
)

# Save descriptive statistics tables
write.csv(original_stats, "data/descriptive_stats_original.csv", row.names = FALSE)
write.csv(complete_stats, "data/descriptive_stats_complete.csv", row.names = FALSE)

print("Descriptive statistics tables saved to /data directory.")
print("Data preparation completed.")
print(paste("Final dataset has", nrow(complete_data), "countries and", ncol(complete_data)-1, "indicators."))

