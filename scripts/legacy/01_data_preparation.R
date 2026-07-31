# Data Preparation Script for GII Analysis
# Thin wrapper calling load_gii_data() and standardize_with_direction()

# Load GII data using centralized function
print("Loading GII data...")
gii_data <- load_gii_data()

# Extract components
idata <- gii_data$raw_data
imeta <- gii_data$imeta
complete_data <- gii_data$complete_data
standardized_data_with_country <- gii_data$standardized_data
indicator_data_only <- gii_data$indicator_data_only
country_names <- gii_data$country_names

print(paste("Loaded data with", nrow(complete_data), "countries and",
            ncol(complete_data)-1, "indicators"))

# Apply direction reversal and standardization
print("Applying direction reversal and standardization...")
final_data <- standardize_with_direction(standardized_data_with_country, imeta)
# Save cleaned and standardized data
print("Saving cleaned and standardized data...")
write.csv(complete_data, "data/gii_data_complete_cases.csv", row.names = FALSE)
write.csv(final_data, "data/gii_data_standardized.csv", row.names = FALSE)
# Generate descriptive statistics for original indicator data
print("Generating descriptive statistics for original data...")
original_indicators <- idata[, -1]  # Exclude uCode column
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
print(paste("Final dataset has", nrow(complete_data), "countries and",
            ncol(complete_data)-1, "indicators."))

