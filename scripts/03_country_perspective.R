###############################################################################
# Country Perspective Analysis for GII Data
# Produces:
#   - Fig 2: Rank correlation of random subsets vs full similarity matrix
#   - Exhaustive search: two strategies (3 per pillar, 1 per sub-pillar),
#     selected on the DISTANCE metric (Spearman correlation of distance
#     matrices, after Example/Matrix_Work.2.R). The same allocation rules are
#     searched on the r_m metric in 04_variable_perspective.R.
#   - SIMPROF + agglomerative clustering: compare full vs subset clustering
# Based on: Example/Matrix_Work.2.R and Example/PCA_Work.2.R
###############################################################################

# ============================================================================
# 0. SETUP
# ============================================================================
# Load required libraries
library(dplyr)
#remotes::install_github("douglaswhitaker/clustsig")
library(clustsig)
library(partitionComparison)  # rand_ind() in functions.R needs this attached
# Source shared functions and config
source("scripts/R/config.R")
source("scripts/R/functions.R")
source("scripts/R/data_utils.R")
# Speed patch for clustsig::simprof (the slow step below). Same arguments, same
# RNG stream, same results - only the inner profile computation is rewritten.
source("scripts/R/simprof_fast.R")
use_fast_simprof()

# Set seed for reproducibility
set.seed(CONFIG$seed)

# ============================================================================
# 1. LOAD DATA
# ============================================================================
message("Loading GII data...")
gii_data <- load_gii_data()
data_standardized <- gii_data$standardized_data
imeta <- gii_data$imeta
country_names <- gii_data$country_names

# Prepare data for analysis (exclude Country column)
indicator_data <- gii_data$indicator_data_only
n_indicators <- ncol(indicator_data)

message(paste("Loaded data with", nrow(data_standardized), "countries and", n_indicators, "indicators"))

# ============================================================================
# 2. FULL DISTANCE MATRIX (BASELINE)
# ============================================================================
message("Computing distance matrix for full dataset...")
indicator_matrix <- as.matrix(indicator_data)
full_distance_matrix <- sim_matrix(indicator_matrix)
message(paste("Distance matrix computed with", attr(full_distance_matrix, "Size"), "countries"))

# ============================================================================
# 3. FIG 2: RANK CORRELATION OF RANDOM SUBSETS
# ============================================================================
message("\n=== FIG 2: Random Subset Analysis ===")

# Create directory for fig2 best subsets if it doesn't exist
fig2_dir <- file.path(CONFIG$paths$outputs_dir, "fig2_best_subsets")
if (!dir.exists(fig2_dir)) {
  dir.create(fig2_dir, recursive = TRUE)
}

# Ranks of the full-dataset distance vector: fixed across every draw and every k
full_distance_ranks <- rank(as.vector(full_distance_matrix))

# Loop through k values (10% to 90% of indicators)
fig2_results <- list()

for (k_idx in seq_along(CONFIG$fig2$k_values)) {
  k <- CONFIG$fig2$k_values[k_idx]
  message(paste("\nProcessing k =", k, "(", round(k / n_indicators * 100), "% of indicators)"))
  
  # Store correlations for this k
  scores <- numeric(CONFIG$fig2$n_random_subsets)
  best_score <- -1
  best_subset <- NULL
  
  # Generate random subsets
  for (i in 1:CONFIG$fig2$n_random_subsets) {
    selected_indices <- sample(1:n_indicators, k, replace = FALSE)
    reduced_data <- indicator_matrix[, selected_indices, drop = FALSE]
    reduced_distance_matrix <- sim_matrix(reduced_data)
    # Spearman is Pearson on ranks, and the full-distance ranks never change,
    # so they are computed once outside the loop instead of on all 1000 draws.
    s <- cor(full_distance_ranks, rank(as.vector(reduced_distance_matrix)))
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
  
  message(paste("  Mean correlation:", round(mean(scores), 4)))
  message(paste("  Max correlation:", round(best_score, 4)))
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
message(paste("Fig 2 results saved to:", fig2_file))

# ============================================================================
# 4. EXHAUSTIVE SEARCH (TWO STRATEGIES)
# ============================================================================
message("\n=== EXHAUSTIVE SEARCH ===")

# --- Strategy 1: Best 3 per pillar (Level 3) ---
message("\nRunning Strategy 1: Best 3 variables per pillar...")
strategy1_results <- best_k_within_group(
  k = CONFIG$exhaustive$strategy1_k_per_pillar,
  group_codes = CONFIG$pillars_level3,
  imeta = imeta,
  data = indicator_data,
  level = 3,
  metric = "distance"
)

strategy1_vars <- unlist(strsplit(strategy1_results$variables, ","))
strategy1_data <- indicator_data[, strategy1_vars, drop = FALSE]
strategy1_distance <- sim_matrix(strategy1_data)
strategy1_correlation <- spearman_dist_cor(full_distance_matrix, strategy1_distance)

message(paste("Strategy 1 (3 per pillar):", length(strategy1_vars), "variables"))
message(paste("Correlation with full dataset:", round(strategy1_correlation, 4)))

strategy1_output <- strategy1_results %>%
  mutate(correlation = strategy1_correlation) %>%
  select(group, best_spearman, variables, correlation)
strategy1_file <- file.path(CONFIG$paths$outputs_dir, "exhaustive_strategy1_pillar.csv")
write.csv(strategy1_output, strategy1_file, row.names = FALSE)
message(paste("Strategy 1 results saved to:", strategy1_file))

# --- Strategy 2: Best 1 per sub-pillar (Level 2) ---
message("\nRunning Strategy 2: Best 1 variable per sub-pillar...")
strategy2_results <- best_k_within_group(
  k = CONFIG$exhaustive$strategy2_k_per_subpillar,
  group_codes = CONFIG$subpillars_level2,
  imeta = imeta,
  data = indicator_data,
  level = 2,
  metric = "distance"
)

strategy2_vars <- unlist(strsplit(strategy2_results$variables, ","))
strategy2_data <- indicator_data[, strategy2_vars, drop = FALSE]
strategy2_distance <- sim_matrix(strategy2_data)
strategy2_correlation <- spearman_dist_cor(full_distance_matrix, strategy2_distance)

message(paste("Strategy 2 (1 per sub-pillar):", length(strategy2_vars), "variables"))
message(paste("Correlation with full dataset:", round(strategy2_correlation, 4)))

strategy2_output <- strategy2_results %>%
  mutate(correlation = strategy2_correlation) %>%
  select(group, best_spearman, variables, correlation)
strategy2_file <- file.path(CONFIG$paths$outputs_dir, "exhaustive_strategy2_subpillar.csv")
write.csv(strategy2_output, strategy2_file, row.names = FALSE)
message(paste("Strategy 2 results saved to:", strategy2_file))

# Rank the strategies by rank correlation; the winner is simply the first.
strategies <- list(
  list(name = "Strategy1_3perPillar",    vars = strategy1_vars, correlation = strategy1_correlation),
  list(name = "Strategy2_1perSubpillar", vars = strategy2_vars, correlation = strategy2_correlation)
)
strategies <- strategies[order(sapply(strategies, `[[`, "correlation"), decreasing = TRUE)]
best_strategy <- strategies[[1]]

message("
Selected best strategy for clustering: ", best_strategy$name)
message("Number of variables: ", length(best_strategy$vars))
message("Correlation with full dataset: ", round(best_strategy$correlation, 4))

# Save selected variables for use in clustering
selected_vars_file <- file.path(CONFIG$paths$outputs_dir, "selected_variables_country_perspective.csv")
selected_vars_df <- data.frame(
  variable = best_strategy$vars,
  stringsAsFactors = FALSE
)
write.csv(selected_vars_df, selected_vars_file, row.names = FALSE)

# ============================================================================
# 5. SIMPROF + AGGLOMERATIVE CLUSTERING
# ============================================================================
message("\n=== SIMPROF + AGGLOMERATIVE CLUSTERING ===")

# Each SIMPROF run is seeded explicitly rather than inheriting wherever the
# stream happens to be. Otherwise reusing the cached full-dataset clustering
# would leave the subset runs at a different point in the RNG stream than a
# from-scratch run, and the two would disagree.

# Run hierarchical clustering on full dataset. This is the single most
# expensive step in the pipeline, and it does not depend on which strategy is
# selected above, so a previous run's result is reused when present - the same
# guard 03b_country_perspective_proportional.R uses. Delete
# outputs/clustering_full.rds to force a recomputation.
full_simprof_file <- file.path(CONFIG$paths$outputs_dir, "clustering_full.rds")
if (file.exists(full_simprof_file)) {
  message(paste("Reusing cached full-dataset clustering from:", full_simprof_file))
  full_results <- readRDS(full_simprof_file)
} else {
  message("Running clustering on full dataset (all indicators)...")
  set.seed(CONFIG$seed)
  full_results <- simprof(
    indicator_data,
    num.expected = CONFIG$clustering$simprof_expected,
    num.simulated = CONFIG$clustering$simprof_simulated,
    method.transform = "identity",
    alpha = CONFIG$clustering$alpha,
    method.cluster = CONFIG$clustering$method_cluster,
    method.distance = CONFIG$clustering$method_dist
  )
  saveRDS(full_results, full_simprof_file)
  message(paste("Full clustering results saved to:", full_simprof_file))
}

# Run hierarchical clustering on each strategy's subset dataset and compute
# its Rand Index against the full clustering. `strategies` is already sorted
# best-first, so downstream consumers that assume row 1 is "the" selected
# strategy keep working.
rand_index_rows <- list()

for (strategy in strategies) {
  strategy_name <- strategy$name
  strategy_vars <- strategy$vars

  message(paste("\nRunning clustering on subset dataset:", strategy_name))
  subset_data <- indicator_data[, strategy_vars, drop = FALSE]
  set.seed(CONFIG$seed)
  subset_results <- simprof(
    subset_data,
    num.expected = CONFIG$clustering$simprof_expected,
    num.simulated = CONFIG$clustering$simprof_simulated,
    method.transform = "identity",
    alpha = CONFIG$clustering$alpha,
    method.cluster = CONFIG$clustering$method_cluster,
    method.distance = CONFIG$clustering$method_dist
  )

  subset_simprof_file <- file.path(CONFIG$paths$outputs_dir, paste0("clustering_subset_", strategy_name, ".rds"))
  saveRDS(subset_results, subset_simprof_file)
  message(paste("Subset clustering results saved to:", subset_simprof_file))

  message("Calculating Rand Index...")
  rand_idx <- rand_ind(full_results, subset_results)
  message(paste("Rand Index:", round(rand_idx, 4)))

  rand_index_rows[[strategy_name]] <- data.frame(
    comparison = paste0("Full_vs_", strategy_name),
    rand_index = rand_idx,
    n_indicators_full = n_indicators,
    n_indicators_subset = length(strategy_vars),
    strategy = strategy_name,
    stringsAsFactors = FALSE
  )
}

# Save Rand Index results for all strategies
rand_index_file <- file.path(CONFIG$paths$outputs_dir, "clustering_rand_index.csv")
rand_index_df <- do.call(rbind, rand_index_rows)
write.csv(rand_index_df, rand_index_file, row.names = FALSE)
message(paste("Rand Index results saved to:", rand_index_file))

message("\n=== Country perspective analysis completed. ===")