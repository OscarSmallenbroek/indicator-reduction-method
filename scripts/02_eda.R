# Exploratory Data Analysis for GII Data
# Produces outputs/pca_components.csv for the report

library(dplyr)
source("scripts/R/config.R")
source("scripts/R/functions.R")
source("scripts/R/data_utils.R")

message("Loading GII data for EDA...")
gii_data <- load_gii_data()
indicator_data <- gii_data$indicator_data_only

message("Loaded data with ", nrow(indicator_data), " countries and ",
        ncol(indicator_data), " indicators")

message("Performing PCA analysis...")
pca_info <- get_pca_components(indicator_data, threshold = CONFIG$pca$variance_threshold)

message("Number of principal components: ", length(pca_info$eigenvalues))

# Number of components needed to reach each variance level
variance_levels <- c(0.8, 0.9, 0.95)
pca_components <- data.frame(
  Variance_Level = paste0(variance_levels * 100, "%"),
  Components_Needed = sapply(
    variance_levels, function(v) which(pca_info$cumulative_variance >= v)[1]
  )
)
print(pca_components)

if (!dir.exists(CONFIG$paths$outputs_dir)) {
  dir.create(CONFIG$paths$outputs_dir, recursive = TRUE)
}
pca_file <- file.path(CONFIG$paths$outputs_dir, "pca_components.csv")
write.csv(pca_components, pca_file, row.names = FALSE)
message("PCA components table saved to ", pca_file)

message("Exploratory data analysis completed.")
