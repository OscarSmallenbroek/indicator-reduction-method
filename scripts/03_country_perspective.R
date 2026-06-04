# Country Perspective Analysis for GII Data

library(dplyr)
library(stats)

print("Loading prepared GII data...")

# Load standardized data
data_standardized <- read.csv("data/gii_data_standardized.csv", stringsAsFactors = FALSE)
print(paste("Loaded data with", nrow(data_standardized), "countries and", ncol(data_standardized)-1, "indicators"))

# Function for finding the similarity/distance matrix
sim_matrix <- function(df) {
  # Note that other methods are possible
  sim_mat <- dist(df, method = "euclidean")
  return(sim_mat)
}

# Prepare data for analysis (exclude Country column)
country_names <- data_standardized$Country
indicator_data <- data_standardized[, -1]  # Remove Country column

print("Computing distance matrix...")
# Compute distance matrix for full dataset
full_distance_matrix <- sim_matrix(indicator_data)
print(paste("Distance matrix computed with", attr(full_distance_matrix, "Size"), "countries"))

# For the country perspective analysis, we need to compare subsets of indicators
# Let's first try with a random subset to establish a baseline

print("Testing with a random subset of indicators...")

# Select random subset of 20 indicators (to demonstrate the method)
set.seed(708557)  # For reproducibility
n_indicators <- ncol(indicator_data)
selected_indices <- sample(1:n_indicators, min(20, n_indicators))
reduced_data <- indicator_data[, selected_indices]

print(paste("Reduced dataset has", ncol(reduced_data), "indicators"))

# Compute distance matrix for reduced dataset
reduced_distance_matrix <- sim_matrix(reduced_data)

# Compare the two distance matrices using correlation
print("Comparing distance matrices...")

# Convert distance objects to vectors for correlation
full_distances <- as.vector(full_distance_matrix)
reduced_distances <- as.vector(reduced_distance_matrix)

# Compute rank correlation (Spearman)
correlation_spearman <- cor(full_distances, reduced_distances, method = "spearman")
correlation_pearson <- cor(full_distances, reduced_distances, method = "pearson")

print(paste("Pearson correlation between distance matrices:", round(correlation_pearson, 4)))
print(paste("Spearman correlation between distance matrices:", round(correlation_spearman, 4)))

# Compare to random subsets to assess statistical significance
print("Assessing statistical significance against random subsets...")

size_test <- function(n_var, full_data, full_distances) {
  # Select random subset of indicators
  selected_indices <- sample(1:ncol(full_data), min(n_var, ncol(full_data)))
  reduced_data <- full_data[, selected_indices]

  # Compute distance matrix for reduced dataset
  reduced_distance_matrix <- sim_matrix(reduced_data)
  reduced_distances <- as.vector(reduced_distance_matrix)

  # Compute correlation
  s <- cor(full_distances, reduced_distances, method = "spearman")
  return(s)
}

# Test different subset sizes
subset_sizes <- c(10, 15, 20, 25, 30)
results <- data.frame(
  size = subset_sizes,
  mean_correlation = numeric(length(subset_sizes)),
  max_correlation = numeric(length(subset_sizes))
)

for (i in 1:length(subset_sizes)) {
  size <- subset_sizes[i]
  print(paste("Testing subset size:", size))

  # Test 100 random subsets of this size
  correlations <- replicate(100, size_test(size, indicator_data, full_distances))

  results$mean_correlation[i] <- mean(correlations)
  results$max_correlation[i] <- max(correlations)

  print(paste("  Mean correlation:", round(mean(correlations), 4)))
  print(paste("  Max correlation:", round(max(correlations), 4)))
}

print("Subset size comparison results:")
print(results)

# Save results
write.csv(results, "outputs/subset_size_comparison.csv", row.names = FALSE)

print("Country perspective analysis completed.")
