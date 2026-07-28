# Test script for PCA-based indicator selection

# Load the functions
source("scripts/select_indicators_by_pca.R")

# Set up test parameters
set.seed(42)
n_obs <- 300
n_factors <- 5
n_indicators <- 20

# Generate synthetic data with known structure
cat("Generating synthetic data...\n")
synthetic_data <- generate_synthetic_data(
  n_obs = n_obs,
  n_factors = n_factors, 
  n_indicators = n_indicators,
  loading_range = c(0.6, 0.9),
  noise_sd = 0.3,
  seed = 42
)

cat("Synthetic data generated successfully.\n")
cat("True loading structure:\n")
print(synthetic_data$true_loadings)
cat("\n")

# Run the indicator selection algorithm
cat("Running indicator selection algorithm...\n")
result <- select_indicators_by_pca(synthetic_data$data)

# Print summary
cat("\n=== Selection Results ===\n")
print_summary(result)

# Verify the selection makes sense
cat("\n=== Verification ===\n")

# Check if selected indicators correspond to the expected structure
true_loadings <- synthetic_data$true_loadings
selected_indices <- result$selected_indicators

cat("Checking if selected indicators match expected structure:\n")
for (i in seq_along(selected_indices)) {
  selected_idx <- selected_indices[i]
  factor_idx <- ceiling(selected_idx / 4)  # Each factor loads on 4 indicators
  cat(sprintf("Selected indicator %d (PC %d) -> True factor: %d\n", 
              selected_idx, result$retained_pc_indices[i], factor_idx))
}

# Additional diagnostics
cat("\n=== Additional Diagnostics ===\n")
cat("Eigenvalues (retained):", result$eigenvalues, "\n")
cat("Total eigenvalues:", result$total_eigenvalues, "\n")
cat("Total indicators:", n_indicators, "\n")
cat("Retained PCs:", result$retained_pcs, "\n")

# Plot explained variance
cat("\nExplained variance plot:\n")
cumulative_variance <- cumsum(result$total_eigenvalues / sum(result$total_eigenvalues))
for (i in seq_along(result$total_eigenvalues)) {
  cat(sprintf("PC %d: %.2f%% cumulative\n", i, cumulative_variance[i] * 100))
}

# Check reconstruction quality
cat("\nReconstruction Quality (R²):", paste(round(result$r_squared, 3), collapse = ", "), "\n")
