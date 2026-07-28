# Exploratory Data Analysis for GII Data

library(readxl)
library(dplyr)

# Load GII data
print("Loading GII data for detailed exploration...")

# Load metadata and data files
imeta <- read_excel("gii-data/imeta.xlsx")
idata <- read_excel("gii-data/idata.xlsx")

print("Data loaded successfully")
# Examine metadata structure
print("=== METADATA STRUCTURE ===")
print(paste("Metadata dimensions:", dim(imeta)[1], "rows,", dim(imeta)[2], "columns"))
print("Metadata column names:")
print(colnames(imeta))
print("First 10 rows of metadata:")
print(head(imeta, 10))

# Examine data structure
print("=== DATA STRUCTURE ===")
print(paste("Data dimensions:", dim(idata)[1], "rows,", dim(idata)[2], "columns"))
print("Data column names:")
print(colnames(idata))
print("First 10 rows of data:")
print(head(idata, 10))
# Check for missing values
print("=== MISSING DATA ANALYSIS ===")
missing_counts <- colSums(is.na(idata))
print("Columns with missing values:")
print(missing_counts[missing_counts > 0])
print(paste("Total missing values:", sum(missing_counts)))

# Identify countries
print("=== COUNTRY INFORMATION ===")
countries <- idata$uCode  # Assuming uCode contains country codes
print(paste("Number of countries/regions:", length(unique(countries))))
print("First 20 countries:")
print(head(unique(countries), 20))

# Identify indicators
indicators <- colnames(idata)[2:ncol(idata)]  # Excluding uCode column
print(paste("Number of indicators:", length(indicators)))
print("First 10 indicators:")
print(head(indicators, 10))

# Check data types
print("=== DATA TYPES ===")
print(sapply(idata, class))

# Summary statistics for numeric columns
print("=== SUMMARY STATISTICS ===")
numeric_data <- idata[, sapply(idata, is.numeric)]
if(ncol(numeric_data) > 0) {
  print("Summary for first 5 numeric indicators:")
  print(summary(numeric_data[, 1:min(5, ncol(numeric_data))]))
}

# Save cleaned data for next steps
print("Saving cleaned data...")
write.csv(idata, "data/gii_data_cleaned.csv", row.names = FALSE)
print("Data exploration completed.")
# Create a list of variables for analysis
variable_list <- data.frame(
  Variable = colnames(idata),
  Type = sapply(idata, class)
)
write.csv(variable_list, "data/variable_list.csv", row.names = FALSE)
print("Variable list saved.")
