# Variable Perspective Analysis for GII Data

library(dplyr)
library(stats)

print("Loading prepared GII data for variable perspective analysis...")

# Load standardized data
data_standardized <- read.csv("data/gii_data_standardized.csv", stringsAsFactors = FALSE)
print(paste("Loaded data with", nrow(data_standardized), "countries and", ncol(data_standardized)-1, "indicators"))

# Prepare data for analysis (exclude Country column)
country_names <- data_standardized$Country
indicator_data <- data_standardized[, -1]  # Remove Country column

print("Performing PCA analysis...")

# Perform PCA on the full dataset
pca_result <- prcomp(indicator_data, center = TRUE, scale. = TRUE)
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
  write.csv(pca_components, "outputs/pca_components.csv", row.names = FALSE)

  print("PCA components table saved to outputs/pca_components.csv")

# Function to calculate weighted r_m for a proposed subset of variables
# This implements Equation 1 from the manuscript:
# r_m = sqrt((sum(i=1 to p)(λ_i * (r_m)_i^2)) / (sum(j=1 to p)λ_j))

# First, we need matrices for multiplication
dm <- as.matrix(indicator_data)
pcm <- as.matrix(pca_result$rotation)
data_pc <- dm %*% pcm

# Function to find the multiple correlation between a PC and a subset of variables
PC_cor <- function(PC, subset_indices, data_pc_matrix, dm_matrix) {
  out <- data_pc_matrix[, PC]
  inputs <- dm_matrix[, subset_indices]
  matrix_for_analysis <- cbind(out, inputs)
  data <- as.data.frame(matrix_for_analysis)
  model <- lm(out ~ ., data = data)
  model_summary <- summary(model)
  mult_r.squared <- model_summary$r.squared
  return(mult_r.squared)
}

# Function to calculate weighted r_m for a proposed subset of variables
r_m <- function(subset_indices, data_pc_matrix, dm_matrix, eig_values_vector) {
  r_sum <- 0
  n_components <- min(length(eig_values_vector), ncol(data_pc_matrix))
  
  for (i in 1:n_components) {
    r_m.i.squared <- PC_cor(i, subset_indices, data_pc_matrix, dm_matrix)
    r_sum <- r_sum + eig_values_vector[i] * r_m.i.squared 
  }
  
  eig_sum = sum(eig_values_vector[1:n_components])
  return(sqrt(r_sum / eig_sum))
}

print("Testing r_m calculation with full dataset...")
# Test with all variables (should be close to 1.0)
all_variables_indices <- 1:ncol(indicator_data)
r_m_all <- r_m(all_variables_indices, data_pc, dm, eig_values)
print(paste("r_m for full dataset:", round(r_m_all, 4)))

# Test with random subsets of different sizes
print("Testing with random subsets...")

set.seed(456)  # For reproducibility
subset_sizes <- c(10, 15, 20, 25, 30)
r_m_results <- data.frame(
  size = subset_sizes,
  mean_r_m = numeric(length(subset_sizes)),
  max_r_m = numeric(length(subset_sizes))
)

for (i in 1:length(subset_sizes)) {
  size <- subset_sizes[i]
  print(paste("Testing subset size:", size))
  
  # Test 20 random subsets of this size (reduced from 50 for faster execution)
  r_m_values <- numeric(20)
  
  for (j in 1:20) {
    # Select random subset of indicators
    selected_indices <- sample(1:ncol(indicator_data), min(size, ncol(indicator_data)))
    r_m_values[j] <- r_m(selected_indices, data_pc, dm, eig_values)
  }
  
  r_m_results$mean_r_m[i] <- mean(r_m_values)
  r_m_results$max_r_m[i] <- max(r_m_values)
  
  print(paste("  Mean r_m:", round(mean(r_m_values), 4)))
  print(paste("  Max r_m:", round(max(r_m_values), 4)))
}

print("Variable subset size comparison results (r_m):")
print(r_m_results)

# Save results
write.csv(r_m_results, "outputs/r_m_subset_size_comparison.csv", row.names = FALSE)

# Identify a good subset using a simpler approach (rather than full greedy)
print("Identifying a good subset of 15 variables...")

# For computational efficiency, we'll use a correlation-based approach
# Select variables that are most representative of the dataset

# Calculate the average absolute correlation of each variable with all others
cor_matrix <- cor(indicator_data)
mean_abs_cor <- apply(cor_matrix, 2, function(x) mean(abs(x)))

# Select top 15 variables with highest mean absolute correlation
top_15_indices <- order(mean_abs_cor, decreasing = TRUE)[1:15]
top_15_r_m <- r_m(top_15_indices, data_pc, dm, eig_values)

print(paste("Top 15 variables by correlation r_m:", round(top_15_r_m, 4)))
print("Selected variables (by index):")
print(top_15_indices)

# Get variable names for the selected indices
variable_names <- colnames(indicator_data)
print("Selected variables (by name):")
selected_variable_names <- variable_names[top_15_indices]
print(selected_variable_names)

# Save the subset information
optimal_subset_info <- data.frame(
  index = top_15_indices,
  variable_name = selected_variable_names
)
write.csv(optimal_subset_info, "outputs/optimal_15_variable_subset.csv", row.names = FALSE)

# Also save just the variable names for easy reference
write.csv(selected_variable_names, "outputs/selected_variables.csv", row.names = FALSE, col.names = FALSE)

print("Variable perspective analysis completed.")