###############################################################################
# Country Perspective Analysis for GII Data
# Produces:
#   - Fig 2: Rank correlation of random subsets vs full similarity matrix
#   - Exhaustive search: two strategies (3 per pillar, 1 per sub-pillar)
#   - SIMPROF + agglomerative clustering: compare full vs subset clustering
# Based on: Example/Matrix_Work.2.R and Example/PCA_Work.2.R
###############################################################################

# ============================================================================
# 0. SETUP
# ============================================================================
# Load required libraries
library(dplyr)
library(cluster)   # for silhouette()

# Source shared functions and config
source("scripts/R/config.R")
source("scripts/R/functions.R")
source("scripts/R/data_utils.R")

# Set seed for reproducibility
set.seed(CONFIG$seed)

# ============================================================================
# 1. LOAD DATA
# ============================================================================
print("Loading GII data...")
gii_data <- load_gii_data()
data_standardized <- gii_data$standardized_data
imeta <- gii_data$imeta
country_names <- gii_data$country_names

# Prepare data for analysis (exclude Country column)
indicator_data <- gii_data$indicator_data_only
n_indicators <- ncol(indicator_data)

print(paste("Loaded data with", nrow(data_standardized), "countries and", n_indicators, "indicators"))

# ============================================================================
# 2. FULL DISTANCE MATRIX (BASELINE)
# ============================================================================
print("Computing distance matrix for full dataset...")
full_distance_matrix <- sim_matrix(indicator_data)
print(paste("Distance matrix computed with", attr(full_distance_matrix, "Size"), "countries"))

# ============================================================================
# 3. FIG 2: RANK CORRELATION OF RANDOM SUBSETS
# ============================================================================
print("\n=== FIG 2: Random Subset Analysis ===")

# Create directory for fig2 best subsets if it doesn't exist
fig2_dir <- file.path(CONFIG$paths$outputs_dir, "fig2_best_subsets")
if (!dir.exists(fig2_dir)) {
  dir.create(fig2_dir, recursive = TRUE)
}

# Loop through k values (10% to 90% of indicators)
fig2_results <- list()

for (k_idx in seq_along(CONFIG$fig2$k_values)) {
  k <- CONFIG$fig2$k_values[k_idx]
  print(paste("\nProcessing k =", k, "(", round(k / n_indicators * 100), "% of indicators)"))
  
  # Store correlations for this k
  scores <- numeric(CONFIG$fig2$n_random_subsets)
  best_score <- -1
  best_subset <- NULL
  
  # Generate random subsets
  for (i in 1:CONFIG$fig2$n_random_subsets) {
    selected_indices <- sample(1:n_indicators, k, replace = FALSE)
    reduced_data <- indicator_data[, selected_indices, drop = FALSE]
    reduced_distance_matrix <- sim_matrix(reduced_data) 
    s <- spearman_dist_cor(full_distance_matrix, reduced_distance_matrix)
    scores[i] <- s
    
    if (s > best_score) {
      best_score <- s
      best_subset <- selected_indices
    }
  }
  
  # Store results for this k
  fig2_results[[as.character(k)]] <- scores
  
  # Save best subset for this k
  best_subset_names <- colnames(indicator_data)[best_subset]
  best_subset_df <- data.frame(
    rank = 1:length(best_subset_names),
    indicator = best_subset_names,
    stringsAsFactors = FALSE
  )
  subset_file <- file.path(fig2_dir, paste0("k_", k, ".csv"))
  write.csv(best_subset_df, subset_file, row.names = FALSE)
  
  print(paste("  Mean correlation:", round(mean(scores), 4)))
  print(paste("  Max correlation:", round(best_score, 4)))
}

# Save overall Fig 2 summary
fig2_summary <- data.frame(
  k = as.integer(names(fig2_results)),
  mean_correlation = sapply(fig2_results, mean),
  median_correlation = sapply(fig2_results, median),
  min_correlation = sapply(fig2_results, min),
  max_correlation = sapply(fig2_results, max),
  stringsAsFactors = FALSE
)

fig2_file <- file.path(CONFIG$paths$outputs_dir, "fig2_rank_correlations.csv")
write.csv(fig2_summary, fig2_file, row.names = FALSE)
print(paste("Fig 2 results saved to:", fig2_file))

# ============================================================================
# 4. EXHAUSTIVE SEARCH (TWO STRATEGIES)
# ============================================================================
print("\n=== EXHAUSTIVE SEARCH ===")

# --- Strategy 1: Best 3 per pillar (Level 3) ---
print("\nRunning Strategy 1: Best 3 variables per pillar...")
strategy1_results <- best_k_within_group(
  k = CONFIG$exhaustive$strategy1_k_per_pillar,
  group_codes = CONFIG$pillars_level3,
  imeta = imeta,
  data = indicator_data,
  level = 3
)

strategy1_vars <- unlist(strsplit(strategy1_results$variables, ","))
strategy1_data <- indicator_data[, strategy1_vars, drop = FALSE]
strategy1_distance <- sim_matrix(strategy1_data)
strategy1_correlation <- spearman_dist_cor(full_distance_matrix, strategy1_distance)

print(paste("Strategy 1 (3 per pillar):", length(strategy1_vars), "variables"))
print(paste("Correlation with full dataset:", round(strategy1_correlation, 4)))

strategy1_output <- strategy1_results %>%
  mutate(correlation = strategy1_correlation) %>%
  select(group, best_rm, variables, correlation)
strategy1_file <- file.path(CONFIG$paths$outputs_dir, "exhaustive_strategy1_pillar.csv")
write.csv(strategy1_output, strategy1_file, row.names = FALSE)
print(paste("Strategy 1 results saved to:", strategy1_file))

# --- Strategy 2: Best 1 per sub-pillar (Level 2) ---
print("\nRunning Strategy 2: Best 1 variable per sub-pillar...")
strategy2_results <- best_k_within_group(
  k = CONFIG$exhaustive$strategy2_k_per_subpillar,
  group_codes = CONFIG$subpillars_level2,
  imeta = imeta,
  data = indicator_data,
  level = 2
)

strategy2_vars <- unlist(strsplit(strategy2_results$variables, ","))
strategy2_data <- indicator_data[, strategy2_vars, drop = FALSE]
strategy2_distance <- sim_matrix(strategy2_data)
strategy2_correlation <- spearman_dist_cor(full_distance_matrix, strategy2_distance)

print(paste("Strategy 2 (1 per sub-pillar):", length(strategy2_vars), "variables"))
print(paste("Correlation with full dataset:", round(strategy2_correlation, 4)))

strategy2_output <- strategy2_results %>%
  mutate(correlation = strategy2_correlation) %>%
  select(group, best_rm, variables, correlation)
strategy2_file <- file.path(CONFIG$paths$outputs_dir, "exhaustive_strategy2_subpillar.csv")
write.csv(strategy2_output, strategy2_file, row.names = FALSE)
print(paste("Strategy 2 results saved to:", strategy2_file))

# Select the better strategy for clustering (higher correlation)
if (strategy1_correlation >= strategy2_correlation) {
  best_strategy_vars <- strategy1_vars
  best_strategy_name <- "Strategy1_3perPillar"
  best_correlation <- strategy1_correlation
} else {
  best_strategy_vars <- strategy2_vars
  best_strategy_name <- "Strategy2_1perSubpillar"
  best_correlation <- strategy2_correlation
}

print(paste("\nSelected best strategy for clustering:", best_strategy_name))
print(paste("Number of variables:", length(best_strategy_vars)))
print(paste("Correlation with full dataset:", round(best_correlation, 4)))

# Save selected variables for use in clustering
selected_vars_file <- file.path(CONFIG$paths$outputs_dir, "selected_variables_country_perspective.csv")
selected_vars_df <- data.frame(
  variable = best_strategy_vars,
  stringsAsFactors = FALSE
)
write.csv(selected_vars_df, selected_vars_file, row.names = FALSE)

# ============================================================================
# 5. SIMPROF + AGGLOMERATIVE CLUSTERING
# ============================================================================
print("\n=== SIMPROF + AGGLOMERATIVE CLUSTERING ===")

# Run hierarchical clustering on full dataset
print("Running clustering on full dataset (all indicators)...")
full_results <- run_simprof(
  indicator_data,
  method.cluster = CONFIG$clustering$method_cluster,
  method.distance = CONFIG$clustering$method_dist
)

# Run hierarchical clustering on subset dataset
print("Running clustering on subset dataset...")
subset_data <- indicator_data[, best_strategy_vars, drop = FALSE]
subset_results <- run_simprof(
  subset_data,
  method.cluster = CONFIG$clustering$method_cluster,
  method.distance = CONFIG$clustering$method_dist
)

# Save clustering results
full_simprof_file <- file.path(CONFIG$paths$outputs_dir, "clustering_full.rds")
saveRDS(full_results, full_simprof_file)
print(paste("Full clustering results saved to:", full_simprof_file))

subset_simprof_file <- file.path(CONFIG$paths$outputs_dir, "clustering_subset.rds")
saveRDS(subset_results, subset_simprof_file)
print(paste("Subset clustering results saved to:", subset_simprof_file))

# Calculate Rand Index
print("Calculating Rand Index...")
rand_idx <- rand_index_from_simprof(full_results, subset_results)
print(paste("Rand Index:", round(rand_idx, 4)))

# Save Rand Index
rand_index_file <- file.path(CONFIG$paths$outputs_dir, "clustering_rand_index.csv")
rand_index_df <- data.frame(
  comparison = "Full_vs_Subset",
  rand_index = rand_idx,
  n_indicators_full = n_indicators,
  n_indicators_subset = length(best_strategy_vars),
  strategy = best_strategy_name,
  stringsAsFactors = FALSE
)
write.csv(rand_index_df, rand_index_file, row.names = FALSE)
print(paste("Rand Index saved to:", rand_index_file))

# Calculate cluster equivalence
print("Calculating cluster equivalence...")
cluster_equiv <- cluster_equivalence(full_results, subset_results)
cluster_equiv_file <- file.path(CONFIG$paths$outputs_dir, "cluster_equivalence.csv")
write.csv(cluster_equiv, cluster_equiv_file, row.names = FALSE)
print(paste("Cluster equivalence saved to:", cluster_equiv_file))

# Print clustering summary
print("\n=== CLUSTERING SUMMARY ===")
print(paste("Full dataset clusters:", full_results$numgroups))
print(paste("Subset dataset clusters:", subset_results$numgroups))
print(paste("Rand Index:", round(rand_idx, 4)))
print(paste("Full silhouette widths:"))
if (!is.null(full_results$silhouette_widths)) {
  print(round(full_results$silhouette_widths, 3))
}
print(paste("Subset silhouette widths:"))
if (!is.null(subset_results$silhouette_widths)) {
  print(round(subset_results$silhouette_widths, 3))
}

print("\n=== Country perspective analysis completed. ===")