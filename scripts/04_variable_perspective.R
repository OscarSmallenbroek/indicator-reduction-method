###############################################################################
# Variable Perspective Analysis for GII Data
# Produces:
#   - Strategy 2: reflective indicator per sub-pillar (best 1 var/sub-pillar by
#     r_m) - the headline 21-variable subset per todo-list.md
#   - Strategy 1: best 3 variables per pillar (flat), for comparison
#   - Strategy 1a: proportional allocation across pillars (largest-remainder
#     apportionment, budget=42), mirroring the country perspective's 1a
#   - Strategy 2a: proportional allocation across sub-pillars (largest-remainder
#     apportionment, budget=42), mirroring the country perspective's 2a
#   - Flat-allocation sweep (1-4 vars/pillar) -> results_by_dimension.csv
#   - Random-subset r_m curve, cross-referenced against the number of PCs
#     needed for the same cumulative variance
#
# Every candidate subset is chosen by exhaustive search within each pillar/
# sub-pillar using r_m against that group's own sub-PCA (the LOCAL proxy), then
# re-scored with the TRUE r_m against the full 78-indicator PCA - the local
# proxy is a tractable stand-in for a combinatorially infeasible full search;
# the true r_m is what's actually reported and compared.
#
# best_k_within_group / best_allocation_within_group are shared with 03/03b,
# which run the SAME allocation rules under metric = "distance". This script
# relies on their metric = "r_m" default; do not change that default without
# passing the metric explicitly here.
#
# Based on: scripts/legacy/04_variable_perspective_{random,bypillar,exhausitive_search}.R
###############################################################################

# ============================================================================
# 0. SETUP
# ============================================================================
library(dplyr)
source("scripts/R/config.R")
source("scripts/R/functions.R")
source("scripts/R/data_utils.R")

set.seed(CONFIG$seed)

# ============================================================================
# 1. LOAD DATA + FULL-DATASET PCA
# ============================================================================
message("Loading GII data...")
gii_data <- load_gii_data()
imeta <- gii_data$imeta
indicator_data <- gii_data$indicator_data_only
n_indicators <- ncol(indicator_data)

message(paste("Loaded data with", nrow(indicator_data), "countries and", n_indicators, "indicators"))

message("Performing PCA on full dataset...")
pca_info <- get_pca_components(indicator_data, threshold = CONFIG$pca$variance_threshold)
eig_values <- pca_info$eigenvalues

dm <- as.matrix(indicator_data)
pcm <- as.matrix(pca_info$pca$rotation)
data_pc <- dm %*% pcm

message(paste("Components needed for", CONFIG$pca$variance_threshold * 100, "% variance:", pca_info$n_components))

# True r_m for a named subset of indicators against the full-dataset PCA
r_m_named <- function(var_names) {
  idx <- match(var_names, colnames(dm))
  r_m(idx, data_pc, dm, eig_values)
}

# ============================================================================
# 2. STRATEGY 2: REFLECTIVE INDICATOR PER SUB-PILLAR (HEADLINE, 21 VARIABLES)
# ============================================================================
message("\n=== STRATEGY 2: Reflective indicator per sub-pillar ===")
strategy2_results <- best_k_within_group(
  k = CONFIG$exhaustive$strategy2_k_per_subpillar,
  group_codes = CONFIG$subpillars_level2,
  imeta = imeta,
  data = indicator_data,
  level = 2
)
strategy2_vars <- unlist(strsplit(strategy2_results$variables, ","))
strategy2_rm <- r_m_named(strategy2_vars)

message(paste("Strategy 2:", length(strategy2_vars), "variables, r_m =", round(strategy2_rm, 4)))

strategy2_output <- strategy2_results %>% mutate(r_m = strategy2_rm)
strategy2_file <- file.path(CONFIG$paths$outputs_dir, "variable_exhaustive_strategy2.csv")
write.csv(strategy2_output, strategy2_file, row.names = FALSE)
message(paste("Strategy 2 results saved to:", strategy2_file))

best_per_subpillar_file <- file.path(CONFIG$paths$outputs_dir, "best_per_subpillar.csv")
write.csv(strategy2_output, best_per_subpillar_file, row.names = FALSE)
message(paste("Reflective indicators saved to:", best_per_subpillar_file))

# ============================================================================
# 3. STRATEGY 1: BEST 3 VARIABLES PER PILLAR (FLAT), FOR COMPARISON
# ============================================================================
message("\n=== STRATEGY 1: Best 3 variables per pillar ===")
strategy1_results <- best_k_within_group(
  k = CONFIG$exhaustive$strategy1_k_per_pillar,
  group_codes = CONFIG$pillars_level3,
  imeta = imeta,
  data = indicator_data,
  level = 3
)
strategy1_vars <- unlist(strsplit(strategy1_results$variables, ","))
strategy1_rm <- r_m_named(strategy1_vars)

message(paste("Strategy 1:", length(strategy1_vars), "variables, r_m =", round(strategy1_rm, 4)))

strategy1_output <- strategy1_results %>% mutate(r_m = strategy1_rm)
strategy1_file <- file.path(CONFIG$paths$outputs_dir, "variable_exhaustive_strategy1.csv")
write.csv(strategy1_output, strategy1_file, row.names = FALSE)
message(paste("Strategy 1 results saved to:", strategy1_file))

# ============================================================================
# 4. STRATEGY 1a: PROPORTIONAL ALLOCATION ACROSS PILLARS
# ============================================================================
message("\n=== STRATEGY 1a: Proportional allocation across pillars ===")
message(paste("Budget:", CONFIG$exhaustive$proportional_target, "indicators (largest-remainder allocation)"))

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
  level = 3
)
strategy1a_vars <- unlist(strsplit(strategy1a_results$variables, ","))
strategy1a_rm <- r_m_named(strategy1a_vars)

message(paste("Strategy 1a:", length(strategy1a_vars), "variables, r_m =", round(strategy1a_rm, 4)))

strategy1a_output <- strategy1a_results %>% mutate(r_m = strategy1a_rm)
strategy1a_file <- file.path(CONFIG$paths$outputs_dir, "variable_exhaustive_strategy1a.csv")
write.csv(strategy1a_output, strategy1a_file, row.names = FALSE)
message(paste("Strategy 1a results saved to:", strategy1a_file))

# ============================================================================
# 5. STRATEGY 2a: PROPORTIONAL ALLOCATION ACROSS SUB-PILLARS
# ============================================================================
message("\n=== STRATEGY 2a: Proportional allocation across sub-pillars ===")
message(paste("Budget:", CONFIG$exhaustive$proportional_target, "indicators (largest-remainder allocation)"))

allocation_2a <- largest_remainder_allocation(
  target = CONFIG$exhaustive$proportional_target,
  imeta = imeta,
  group_codes = CONFIG$subpillars_level2,
  level = 2
)
message("Sub-pillar allocation (2a):")
print(allocation_2a)

strategy2a_results <- best_allocation_within_group(
  allocation = allocation_2a,
  imeta = imeta,
  data = indicator_data,
  level = 2
)
strategy2a_vars <- unlist(strsplit(strategy2a_results$variables, ","))
strategy2a_rm <- r_m_named(strategy2a_vars)

message(paste("Strategy 2a:", length(strategy2a_vars), "variables, r_m =", round(strategy2a_rm, 4)))

strategy2a_output <- strategy2a_results %>% mutate(r_m = strategy2a_rm)
strategy2a_file <- file.path(CONFIG$paths$outputs_dir, "variable_exhaustive_strategy2a.csv")
write.csv(strategy2a_output, strategy2a_file, row.names = FALSE)
message(paste("Strategy 2a results saved to:", strategy2a_file))

# ============================================================================
# 5b. STRATEGY 1b: PROPORTIONAL ACROSS PILLARS AT THE FLAT BUDGET
# ============================================================================
# Matched-budget control for the flat-vs-proportional comparison: same rule as
# 1a, same 21-indicator budget as Strategy 1. No sub-pillar counterpart - at 21
# that apportionment is exactly 1 per sub-pillar, i.e. Strategy 2 itself.
message("
=== STRATEGY 1b: Proportional across pillars at the flat budget ===")
message(paste("Budget:", CONFIG$exhaustive$proportional_target_matched, "indicators (largest-remainder allocation)"))

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
  level = 3
)
strategy1b_vars <- unlist(strsplit(strategy1b_results$variables, ","))
strategy1b_rm <- r_m_named(strategy1b_vars)

message(paste("Strategy 1b:", length(strategy1b_vars), "variables, r_m =", round(strategy1b_rm, 4)))

strategy1b_output <- strategy1b_results %>% mutate(r_m = strategy1b_rm)
strategy1b_file <- file.path(CONFIG$paths$outputs_dir, "variable_exhaustive_strategy1b.csv")
write.csv(strategy1b_output, strategy1b_file, row.names = FALSE)
message(paste("Strategy 1b results saved to:", strategy1b_file))

# ============================================================================
# 6. FLAT-ALLOCATION SWEEP (1-4 VARIABLES PER PILLAR)
# ============================================================================
message("\n=== FLAT-ALLOCATION SWEEP (1-4 variables per pillar) ===")
dimension_results <- list()
for (k in 1:4) {
  if (k == CONFIG$exhaustive$strategy1_k_per_pillar) {
    # Section 3 already searched this k - identical arguments, identical result.
    message("Reusing Strategy 1 search for ", k, " indicator(s) per pillar...")
    k_results <- strategy1_results
    k_rm <- strategy1_rm
  } else {
    message("Searching for ", k, " indicator(s) per pillar...")
    k_results <- best_k_within_group(
      k = k, group_codes = CONFIG$pillars_level3, imeta = imeta,
      data = indicator_data, level = 3
    )
    k_rm <- r_m_named(unlist(strsplit(k_results$variables, ",")))
  }
  dimension_results[[as.character(k)]] <- k_results %>%
    transmute(best_rm, indicators = variables, component = group,
              n_indicators = k, total_rm = k_rm)
}
results_by_dimension <- bind_rows(dimension_results)
dimension_file <- file.path(CONFIG$paths$outputs_dir, "results_by_dimension.csv")
write.csv(results_by_dimension, dimension_file, row.names = FALSE)
message(paste("Flat-allocation sweep saved to:", dimension_file))

# ============================================================================
# 7. RANDOM-SUBSET r_m CURVE (+ COMPARISON AGAINST #PCS NEEDED)
# ============================================================================
message("\n=== RANDOM-SUBSET r_m CURVE ===")
n_draws <- CONFIG$sampling$n_random_subsets
k_values <- CONFIG$fig2$k_values
message(paste("Testing", length(k_values), "subset sizes x", n_draws, "random draws each"))

curve_results <- list()
for (k in k_values) {
  message(paste("k =", k, "(", round(k / n_indicators * 100), "% of indicators)"))
  scores <- numeric(n_draws)
  for (i in 1:n_draws) {
    idx <- sample(1:n_indicators, k, replace = FALSE)
    scores[i] <- r_m(idx, data_pc, dm, eig_values)
  }
  curve_results[[as.character(k)]] <- data.frame(
    k = k,
    mean_r_m = mean(scores),
    median_r_m = median(scores),
    min_r_m = min(scores),
    max_r_m = max(scores)
  )
  message(paste("  Mean r_m:", round(mean(scores), 4)))
}
r_m_subset_size_comparison <- bind_rows(curve_results)
curve_file <- file.path(CONFIG$paths$outputs_dir, "r_m_subset_size_comparison.csv")
write.csv(r_m_subset_size_comparison, curve_file, row.names = FALSE)
message(paste("Random-subset r_m curve saved to:", curve_file))

# Cross-reference against the number of principal components needed for the
# same cumulative variance (pca_info$cumulative_variance[k] = variance
# explained by the first k PCs - unconstrained linear combinations, not raw
# indicators, so this is an upper bound rather than a like-for-like target)
r_m_vs_pca <- r_m_subset_size_comparison %>%
  mutate(
    subset_variance_captured = mean_r_m^2,
    pca_variance_at_k_components = pca_info$cumulative_variance[k]
  ) %>%
  select(k, mean_r_m, subset_variance_captured, pca_variance_at_k_components)
pca_comparison_file <- file.path(CONFIG$paths$outputs_dir, "r_m_vs_pca_components.csv")
write.csv(r_m_vs_pca, pca_comparison_file, row.names = FALSE)
message(paste("r_m vs. #PCs comparison saved to:", pca_comparison_file))

# ============================================================================
# 8. HEADLINE SUBSET (STRATEGY 2, 21 REFLECTIVE INDICATORS)
# ============================================================================
message("\n=== HEADLINE SUBSET ===")
selected_vars_df <- data.frame(variable = strategy2_vars, stringsAsFactors = FALSE)
selected_vars_file <- file.path(CONFIG$paths$outputs_dir, "selected_variables.csv")
write.csv(selected_vars_df, selected_vars_file, row.names = FALSE)
message(paste("Selected variables saved to:", selected_vars_file))

optimal_subset_info <- strategy2_results %>%
  transmute(subpillar = group, variable = variables, subpillar_r_m = best_rm) %>%
  mutate(overall_r_m = strategy2_rm)
optimal_subset_file <- file.path(CONFIG$paths$outputs_dir, "optimal_21_variable_subset.csv")
write.csv(optimal_subset_info, optimal_subset_file, row.names = FALSE)
message(paste("Optimal 21-variable subset saved to:", optimal_subset_file))

# ============================================================================
# 9. PILLAR ALLOCATION FOR BUDGET = 20 (feeds the QMD combinatorics example)
# ============================================================================
message("\n=== PILLAR ALLOCATION (budget = 20) ===")
pillar_counts_20 <- group_indicator_counts(CONFIG$pillars_level3, imeta, level = 3)
allocation_20 <- largest_remainder_allocation(20, imeta, CONFIG$pillars_level3, level = 3)
allocation20_df <- data.frame(
  pillar = names(pillar_counts_20),
  n_indicators = as.integer(pillar_counts_20),
  allocation = as.integer(allocation_20[names(pillar_counts_20)]),
  stringsAsFactors = FALSE
)
allocation20_file <- file.path(CONFIG$paths$outputs_dir, "allocation20.rds")
saveRDS(allocation20_df, allocation20_file)
print(allocation20_df)
message(paste("Pillar allocation (budget=20) saved to:", allocation20_file))

message("\n=== Variable perspective analysis completed. ===")
