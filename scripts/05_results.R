# Results Synthesis and Reporting

# This script will:
# 1. Compile results from previous analyses
# 2. Generate visualizations
# 3. Create summary tables
# 4. Perform final comparisons

library(ggplot2)
library(dplyr)
library(readr)

# Create outputs directory if it doesn't exist
if (!dir.exists("outputs")) {
  dir.create("outputs")
}

print("Compiling results from GII Index Reduction Analysis...")

# Read in the results
try({
  r_m_results <- read_csv("outputs/r_m_subset_size_comparison.csv")
  print("

Subset size (r_m) results loaded.")
  print(r_m_results)
}, silent = TRUE)

try({
  country_results <- read_csv("outputs/subset_size_comparison.csv")
  print("

Subset size (country perspective) results loaded.")
  print(country_results)
}, silent = TRUE)

try({
  variable_subset_info <- read_csv("outputs/optimal_15_variable_subset.csv")
  print("

Optimal variable subset information loaded.")
  print(variable_subset_info)
}, silent = TRUE)

try({
  recommended_vars <- read_csv("outputs/selected_variables.csv", col_names = FALSE)
  print("

Recommended indicators:")
  print(recommended_vars)
}, silent = TRUE)

# Create visualizations
if (exists("r_m_results")) {
  # Plot for r_m comparison
  r_m_plot <- ggplot(r_m_results, aes(x = size, y = mean_r_m)) +
    geom_point() +
    geom_line() +
    geom_hline(yintercept = 0.85, linetype = "dashed", color = "red") +
    labs(title = "GII Index: Variable Subset Size vs r_m",
         x = "Number of Indicators", y = "r_m value") +
    theme_minimal()

  # Save the plot
  ggsave("outputs/r_m_results_plot.png", plot = r_m_plot, width = 8, height = 6)
  print("r_m comparison plot saved.")
}

