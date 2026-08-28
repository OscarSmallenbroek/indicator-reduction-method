###############################################################################
# Country Perspective Analysis for GII Data - Proportional Allocation
# Produces:
#   - Exhaustive search: three proportional-allocation strategies
#       1a: budget allocated across pillars (Level 3) by largest-remainder
#           apportionment, instead of a flat 3-per-pillar (Strategy 1)
#       2a: budget allocated across sub-pillars (Level 2) by largest-remainder
#           apportionment, instead of a flat 1-per-sub-pillar (Strategy 2)
#       1b: the 1a rule at the flat strategies' 21-indicator budget, so that
#           flat vs proportional can be compared without budget confounded
#           into it (pillar level only - see config for why)
#   - SIMPROF + agglomerative clustering: compare full vs 1a/2a-subset clustering
#
# Selection uses the DISTANCE metric (Spearman correlation of distance matrices,
# after Example/Matrix_Work.2.R), as in 03. The same proportional allocations
# are searched on the r_m metric in 04_variable_perspective.R.
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
library(clustsig)
library(partitionComparison)  # rand_ind() in functions.R needs this attached
source("scripts/R/config.R")
source("scripts/R/functions.R")
source("scripts/R/data_utils.R")
# Speed patch for clustsig::simprof (the slow step below). Same arguments, same
# RNG stream, same results - only the inner profile computation is rewritten.
source("scripts/R/simprof_fast.R")
use_fast_simprof()

set.seed(CONFIG$seed)

# ============================================================================
# 1. LOAD DATA
# ============================================================================
message("Loading GII data...")
gii_data <- load_gii_data()
data_standardized <- gii_data$standardized_data
imeta <- gii_data$imeta
country_names <- gii_data$country_names

indicator_data <- gii_data$indicator_data_only
n_indicators <- ncol(indicator_data)

message(paste("Loaded data with", nrow(data_standardized), "countries and", n_indicators, "indicators"))

# ============================================================================
# 2. FULL DISTANCE MATRIX (BASELINE)
# ============================================================================
message("Computing distance matrix for full dataset...")
full_distance_matrix <- sim_matrix(indicator_data)
message(paste("Distance matrix computed with", attr(full_distance_matrix, "Size"), "countries"))

# ============================================================================
# 3. PROPORTIONAL EXHAUSTIVE SEARCH (STRATEGIES 1a AND 2a)
# ============================================================================
message("\n=== PROPORTIONAL EXHAUSTIVE SEARCH ===")
message(paste("Budget:", CONFIG$exhaustive$proportional_target, "indicators (largest-remainder allocation)"))

# --- Strategy 1a: proportional allocation across pillars (Level 3) ---
message("\nRunning Strategy 1a: proportional allocation across pillars...")
allocation_1a <- largest_remainder_allocation(
  target = CONFIG$exhaustive$proportional_target,
  imeta = imeta,
  group_codes = CONFIG$pillars_level3,
  level = 3
)
message("Pillar allocation (1a):")
print(allocation_1a)

strategy1a_results <- best_allocation_within_group(
  allocation = allocation_1a,
  imeta = imeta,
  data = indicator_data,
  level = 3,
  metric = "distance"
)

strategy1a_vars <- unlist(strsplit(strategy1a_results$variables, ","))
strategy1a_data <- indicator_data[, strategy1a_vars, drop = FALSE]
strategy1a_distance <- sim_matrix(strategy1a_data)
strategy1a_correlation <- spearman_dist_cor(full_distance_matrix, strategy1a_distance)

message(paste("Strategy 1a (proportional per pillar):", length(strategy1a_vars), "variables"))
message(paste("Correlation with full dataset:", round(strategy1a_correlation, 4)))

strategy1a_output <- strategy1a_results %>%
  mutate(correlation = strategy1a_correlation) %>%
  select(group, allocated, best_spearman, variables, correlation)
strategy1a_file <- file.path(CONFIG$paths$outputs_dir, "exhaustive_strategy1a_pillar.csv")
write.csv(strategy1a_output, strategy1a_file, row.names = FALSE)
message(paste("Strategy 1a results saved to:", strategy1a_file))

# --- Strategy 2a: proportional allocation across sub-pillars (Level 2) ---
message("\nRunning Strategy 2a: proportional allocation across sub-pillars...")
allocation_2a <- largest_remainder_allocation(
  target = CONFIG$exhaustive$proportional_target,
  imeta = imeta,
  group_codes = CONFIG$subpillars_level2,
  level = 2
)
message("Sub-pillar allocation (2a):")
print(allocation_2a)
if (any(allocation_2a == 0)) {
  message(paste("Sub-pillars excluded (allocated 0):",
              paste(names(allocation_2a)[allocation_2a == 0], collapse = ", ")))
}

strategy2a_results <- best_allocation_within_group(
  allocation = allocation_2a,
  imeta = imeta,
  data = indicator_data,
  level = 2,
  metric = "distance"
)

strategy2a_vars <- unlist(strsplit(strategy2a_results$variables, ","))
strategy2a_data <- indicator_data[, strategy2a_vars, drop = FALSE]
strategy2a_distance <- sim_matrix(strategy2a_data)
strategy2a_correlation <- spearman_dist_cor(full_distance_matrix, strategy2a_distance)

message(paste("Strategy 2a (proportional per sub-pillar):", length(strategy2a_vars), "variables"))
message(paste("Correlation with full dataset:", round(strategy2a_correlation, 4)))

strategy2a_output <- strategy2a_results %>%
  mutate(correlation = strategy2a_correlation) %>%
  select(group, allocated, best_spearman, variables, correlation)
strategy2a_file <- file.path(CONFIG$paths$outputs_dir, "exhaustive_strategy2a_subpillar.csv")
write.csv(strategy2a_output, strategy2a_file, row.names = FALSE)
message(paste("Strategy 2a results saved to:", strategy2a_file))

# --- Strategy 1b: proportional allocation across pillars at the FLAT budget ---
# Matched-budget control for the flat-vs-proportional comparison: same rule as
# 1a, same 21-indicator budget as Strategy 1. There is no sub-pillar
# counterpart - at 21 that apportionment is exactly 1 per sub-pillar, i.e.
# Strategy 2 itself (see CONFIG$exhaustive$proportional_target_matched).
message("
Running Strategy 1b: proportional allocation across pillars at the flat budget...")
allocation_1b <- largest_remainder_allocation(
  target = CONFIG$exhaustive$proportional_target_matched,
  imeta = imeta,
  group_codes = CONFIG$pillars_level3,
  level = 3
)
message("Pillar allocation (1b):")
print(allocation_1b)

strategy1b_results <- best_allocation_within_group(
  allocation = allocation_1b,
  imeta = imeta,
  data = indicator_data,
  level = 3,
  metric = "distance"
)

strategy1b_vars <- unlist(strsplit(strategy1b_results$variables, ","))
strategy1b_correlation <- spearman_dist_cor(
  full_distance_matrix, sim_matrix(indicator_data[, strategy1b_vars, drop = FALSE])
)

message(paste("Strategy 1b (proportional per pillar, matched budget):",
              length(strategy1b_vars), "variables"))
message(paste("Correlation with full dataset:", round(strategy1b_correlation, 4)))

strategy1b_output <- strategy1b_results %>%
  mutate(correlation = strategy1b_correlation) %>%
  select(group, allocated, best_spearman, variables, correlation)
strategy1b_file <- file.path(CONFIG$paths$outputs_dir, "exhaustive_strategy1b_pillar.csv")
write.csv(strategy1b_output, strategy1b_file, row.names = FALSE)
message(paste("Strategy 1b results saved to:", strategy1b_file))

# Rank the proportional strategies by rank correlation; the winner is first.
strategies <- list(
  list(name = "Strategy1a_ProportionalPillar",    vars = strategy1a_vars, correlation = strategy1a_correlation),
  list(name = "Strategy2a_ProportionalSubpillar", vars = strategy2a_vars, correlation = strategy2a_correlation)
)
strategies <- strategies[order(sapply(strategies, `[[`, "correlation"), decreasing = TRUE)]
best_strategy <- strategies[[1]]

message("
Selected best proportional strategy for clustering: ", best_strategy$name)
message("Number of variables: ", length(best_strategy$vars))
message("Correlation with full dataset: ", round(best_strategy$correlation, 4))

# Save selected variables for use in clustering
selected_vars_file <- file.path(CONFIG$paths$outputs_dir, "selected_variables_country_perspective_proportional.csv")
selected_vars_df <- data.frame(
  variable = best_strategy$vars,
  stringsAsFactors = FALSE
)
write.csv(selected_vars_df, selected_vars_file, row.names = FALSE)

# ============================================================================
# 4. SIMPROF + AGGLOMERATIVE CLUSTERING
# ============================================================================
message("\n=== SIMPROF + AGGLOMERATIVE CLUSTERING ===")

# Each SIMPROF run is seeded explicitly rather than inheriting wherever the
# stream happens to be. Otherwise reusing the cached full-dataset clustering
# would leave the subset runs at a different point in the RNG stream than a
# from-scratch run, and the two would disagree.

# Reuse the full-dataset clustering from 03_country_perspective.R if it exists
# (identical computation - SIMPROF is the slow step, no need to repeat it).
full_simprof_file <- file.path(CONFIG$paths$outputs_dir, "clustering_full.rds")
if (file.exists(full_simprof_file)) {
  message(paste("Reusing cached full-dataset clustering from:", full_simprof_file))
  full_results <- readRDS(full_simprof_file)
} else {
  message("No cached full-dataset clustering found - running clustering on full dataset...")
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

# Run hierarchical clustering on the selected proportional subset
message("Running clustering on proportional subset dataset...")
subset_data <- indicator_data[, best_strategy$vars, drop = FALSE]
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

subset_simprof_file <- file.path(CONFIG$paths$outputs_dir, "clustering_subset_proportional.rds")
saveRDS(subset_results, subset_simprof_file)
message(paste("Proportional subset clustering results saved to:", subset_simprof_file))

# Calculate Rand Index
message("Calculating Rand Index...")
rand_idx <- rand_ind(full_results, subset_results)
message(paste("Rand Index:", round(rand_idx, 4)))

# Save Rand Index
rand_index_file <- file.path(CONFIG$paths$outputs_dir, "clustering_rand_index_proportional.csv")
rand_index_df <- data.frame(
  comparison = "Full_vs_ProportionalSubset",
  rand_index = rand_idx,
  n_indicators_full = n_indicators,
  n_indicators_subset = length(best_strategy$vars),
  strategy = best_strategy$name,
  stringsAsFactors = FALSE
)
write.csv(rand_index_df, rand_index_file, row.names = FALSE)
message(paste("Rand Index saved to:", rand_index_file))

message("\n=== Proportional country perspective analysis completed. ===")
