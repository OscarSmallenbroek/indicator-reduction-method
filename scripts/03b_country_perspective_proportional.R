###############################################################################
# Country Perspective Analysis for GII Data - Proportional Allocation
# Produces:
#   - Exhaustive search: two proportional-allocation strategies
#       1a: budget allocated across pillars (Level 3) by largest-remainder
#           apportionment, instead of a flat 3-per-pillar (Strategy 1)
#       2a: budget allocated across sub-pillars (Level 2) by largest-remainder
#           apportionment, instead of a flat 1-per-sub-pillar (Strategy 2)
#   - SIMPROF + agglomerative clustering: compare full vs 1a/2a-subset clustering
#
# Fig 2 (random-subset rank correlation) is NOT reproduced here - it doesn't
# depend on which exhaustive strategy is used, so it's only computed in
# 03_country_perspective.R. The full-dataset SIMPROF clustering is reused from
# outputs/clustering_full.rds if 03_country_perspective.R already produced it
# (identical computation, and SIMPROF is the slow step).
###############################################################################

# ============================================================================
# 0. SETUP
# ============================================================================
library(dplyr)
library(cluster)   # for silhouette()
library(clustsig)
library(partitionComparison)  # rand_ind() in functions.R needs this attached
source("scripts/R/config.R")
source("scripts/R/functions.R")
source("scripts/R/data_utils.R")

set.seed(CONFIG$seed)

# ============================================================================
# 1. LOAD DATA
# ============================================================================
print("Loading GII data...")
gii_data <- load_gii_data()
data_standardized <- gii_data$standardized_data
imeta <- gii_data$imeta
country_names <- gii_data$country_names

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
# 3. PROPORTIONAL EXHAUSTIVE SEARCH (STRATEGIES 1a AND 2a)
# ============================================================================
print("\n=== PROPORTIONAL EXHAUSTIVE SEARCH ===")
print(paste("Budget:", CONFIG$exhaustive$proportional_target, "indicators (largest-remainder allocation)"))

# --- Strategy 1a: proportional allocation across pillars (Level 3) ---
print("\nRunning Strategy 1a: proportional allocation across pillars...")
allocation_1a <- largest_remainder_allocation(
  target = CONFIG$exhaustive$proportional_target,
  imeta = imeta,
  group_codes = CONFIG$pillars_level3,
  level = 3
)
print("Pillar allocation (1a):")
print(allocation_1a)

strategy1a_results <- best_allocation_within_group(
  allocation = allocation_1a,
  imeta = imeta,
  data = indicator_data,
  level = 3
)

strategy1a_vars <- unlist(strsplit(strategy1a_results$variables, ","))
strategy1a_data <- indicator_data[, strategy1a_vars, drop = FALSE]
strategy1a_distance <- sim_matrix(strategy1a_data)
strategy1a_correlation <- spearman_dist_cor(full_distance_matrix, strategy1a_distance)

print(paste("Strategy 1a (proportional per pillar):", length(strategy1a_vars), "variables"))
print(paste("Correlation with full dataset:", round(strategy1a_correlation, 4)))

strategy1a_output <- strategy1a_results %>%
  mutate(correlation = strategy1a_correlation) %>%
  select(group, allocated, best_rm, variables, correlation)
strategy1a_file <- file.path(CONFIG$paths$outputs_dir, "exhaustive_strategy1a_pillar.csv")
write.csv(strategy1a_output, strategy1a_file, row.names = FALSE)
print(paste("Strategy 1a results saved to:", strategy1a_file))

# --- Strategy 2a: proportional allocation across sub-pillars (Level 2) ---
print("\nRunning Strategy 2a: proportional allocation across sub-pillars...")
allocation_2a <- largest_remainder_allocation(
  target = CONFIG$exhaustive$proportional_target,
  imeta = imeta,
  group_codes = CONFIG$subpillars_level2,
  level = 2
)
print("Sub-pillar allocation (2a):")
print(allocation_2a)
if (any(allocation_2a == 0)) {
  print(paste("Sub-pillars excluded (allocated 0):",
              paste(names(allocation_2a)[allocation_2a == 0], collapse = ", ")))
}

strategy2a_results <- best_allocation_within_group(
  allocation = allocation_2a,
  imeta = imeta,
  data = indicator_data,
  level = 2
)

strategy2a_vars <- unlist(strsplit(strategy2a_results$variables, ","))
strategy2a_data <- indicator_data[, strategy2a_vars, drop = FALSE]
strategy2a_distance <- sim_matrix(strategy2a_data)
strategy2a_correlation <- spearman_dist_cor(full_distance_matrix, strategy2a_distance)

print(paste("Strategy 2a (proportional per sub-pillar):", length(strategy2a_vars), "variables"))
print(paste("Correlation with full dataset:", round(strategy2a_correlation, 4)))

strategy2a_output <- strategy2a_results %>%
  mutate(correlation = strategy2a_correlation) %>%
  select(group, allocated, best_rm, variables, correlation)
strategy2a_file <- file.path(CONFIG$paths$outputs_dir, "exhaustive_strategy2a_subpillar.csv")
write.csv(strategy2a_output, strategy2a_file, row.names = FALSE)
print(paste("Strategy 2a results saved to:", strategy2a_file))

# Select the better strategy for clustering (higher correlation)
if (strategy1a_correlation >= strategy2a_correlation) {
  best_strategy_vars <- strategy1a_vars
  best_strategy_name <- "Strategy1a_ProportionalPillar"
  best_correlation <- strategy1a_correlation
} else {
  best_strategy_vars <- strategy2a_vars
  best_strategy_name <- "Strategy2a_ProportionalSubpillar"
  best_correlation <- strategy2a_correlation
}

print(paste("\nSelected best proportional strategy for clustering:", best_strategy_name))
print(paste("Number of variables:", length(best_strategy_vars)))
print(paste("Correlation with full dataset:", round(best_correlation, 4)))

# Save selected variables for use in clustering
selected_vars_file <- file.path(CONFIG$paths$outputs_dir, "selected_variables_country_perspective_proportional.csv")
selected_vars_df <- data.frame(
  variable = best_strategy_vars,
  stringsAsFactors = FALSE
)
write.csv(selected_vars_df, selected_vars_file, row.names = FALSE)

# ============================================================================
# 4. SIMPROF + AGGLOMERATIVE CLUSTERING
# ============================================================================
print("\n=== SIMPROF + AGGLOMERATIVE CLUSTERING ===")

# Reuse the full-dataset clustering from 03_country_perspective.R if it exists
# (identical computation - SIMPROF is the slow step, no need to repeat it).
full_simprof_file <- file.path(CONFIG$paths$outputs_dir, "clustering_full.rds")
if (file.exists(full_simprof_file)) {
  print(paste("Reusing cached full-dataset clustering from:", full_simprof_file))
  full_results <- readRDS(full_simprof_file)
} else {
  print("No cached full-dataset clustering found - running clustering on full dataset...")
  full_results <- simprof(
    indicator_data,
    num.expected = 1000, num.simulated = 999,
    method.transform = "identity", alpha = 0.05,
    method.cluster = CONFIG$clustering$method_cluster,
    method.distance = CONFIG$clustering$method_dist
  )
  saveRDS(full_results, full_simprof_file)
  print(paste("Full clustering results saved to:", full_simprof_file))
}

# Run hierarchical clustering on the selected proportional subset
print("Running clustering on proportional subset dataset...")
subset_data <- indicator_data[, best_strategy_vars, drop = FALSE]
subset_results <- simprof(
  subset_data,
  num.expected = 1000, num.simulated = 999,
  method.transform = "identity", alpha = 0.05,
  method.cluster = CONFIG$clustering$method_cluster,
  method.distance = CONFIG$clustering$method_dist
)

subset_simprof_file <- file.path(CONFIG$paths$outputs_dir, "clustering_subset_proportional.rds")
saveRDS(subset_results, subset_simprof_file)
print(paste("Proportional subset clustering results saved to:", subset_simprof_file))

# Calculate Rand Index
print("Calculating Rand Index...")
rand_idx <- rand_ind(full_results, subset_results)
print(paste("Rand Index:", round(rand_idx, 4)))

# Save Rand Index
rand_index_file <- file.path(CONFIG$paths$outputs_dir, "clustering_rand_index_proportional.csv")
rand_index_df <- data.frame(
  comparison = "Full_vs_ProportionalSubset",
  rand_index = rand_idx,
  n_indicators_full = n_indicators,
  n_indicators_subset = length(best_strategy_vars),
  strategy = best_strategy_name,
  stringsAsFactors = FALSE
)
write.csv(rand_index_df, rand_index_file, row.names = FALSE)
print(paste("Rand Index saved to:", rand_index_file))

print("\n=== Proportional country perspective analysis completed. ===")
