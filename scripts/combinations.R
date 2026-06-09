# Calculate unique combinations for index subpillars

# Define the subpillars and their parameters
subpillars <- c("IN.1", "IN.2", "IN.3", "IN.4", "IN.5", "OUT.6", "OUT.7")

# Number of items to select from each subpillar
select_k <- c(2, 3, 2, 3, 4, 4, 3)

# Total number of indicators available in each subpillar
total_n <- c(6, 12, 9, 11, 15, 14, 11)

# Create a data frame for clarity
index_data <- data.frame(
  subpillar = subpillars,
  n_indicators = total_n,
  select_k = select_k
)

cat("Index Subpillar Combination Analysis\n")
cat("=====================================\n\n")
print(index_data)
cat("\n")

# Calculate combinations for each subpillar
combinations <- numeric(length(subpillars))
for (i in seq_along(subpillars)) {
  combinations[i] <- choose(total_n[i], select_k[i])
}

# Display individual combinations
cat("Combinations per subpillar:\n")
for (i in seq_along(subpillars)) {
  cat(sprintf(\