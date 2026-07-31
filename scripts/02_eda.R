# Exploratory Data Analysis for GII Data
# Thin wrapper that produces pca_components.csv for QMD

library(dplyr)

print("Loading GII data for EDA...")

# Load data using centralized function
gii_data <- load_gii_data()
final_data <- gii_data$standardized_data
indicator_data_only <- gii_data$indicator_data_only

print(paste("Loaded data with", nrow(final_data), "countries and", 
            ncol(final_data)-1, "indicators"))

# Perform PCA analysis
print("Performing PCA analysis...")
pca_result <- prcomp(indicator_data_only, center = TRUE, scale. = TRUE)
summary_pca <- summary(pca_result)

# Eigenvalues in decreasing order
eig_values <- pca_result$sdev^2

print("PCA Results:")
print(paste("Number of principal components:", length(eig_values)))
print("Proportion of variance explained by first 10 components:")
print(round(summary_pca$importance[2, 1:min(10, length(eig_values))], 4))

# Calculate cumulative variance explained
cumulative_variance <- cumsum(summary_pca$importance[2, ])
components_80_pct <- which(cumulative_variance >= 0.8)[1]
components_90_pct <- which(cumulative_variance >= 0.9)[1]
components_95_pct <- which(cumulative_variance >= 0.95)[1]

print(paste("Components needed for 80% variance:", components_80_pct))
print(paste("Components needed for 90% variance:", components_90_pct))
print(paste("Components needed for 95% variance:", components_95_pct))

# Create table for components needed
pca_components <- data.frame(
  Variance_Level = c("80%", "90%", "95%"),
  Components_Needed = c(components_80_pct, components_90_pct, components_95_pct)
)

# Export PCA components table
output_dir <- "outputs"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

write.csv(pca_components, file.path(output_dir, "pca_components.csv"), row.names = FALSE)
print("PCA components table saved to outputs/pca_components.csv")

print("Exploratory data analysis completed.")
