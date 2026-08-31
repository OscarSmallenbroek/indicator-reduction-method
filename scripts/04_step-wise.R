###############################################################################
# Step-wise Optimization for GII Data
#
# An alternative to the exhaustive search in 03/03b/04: instead of enumerating
# every combination within each group (best_allocation_within_group()), start
# from a random per-group selection and greedily swap indicators in and out,
# scoring each candidate 21/42-indicator subset against the WHOLE dataset -
# not locally per group - keeping any swap that improves the score, until a
# full pass makes no improvement. See stepwise_search() in
# scripts/R/functions.R, which reproduces best_n_step() from
# replication/Oscar_stepwise.R (verified to match it bit-for-bit under a fixed
# seed) generalized to any allocation, level and metric. Exhaustive search
# stays authoritative wherever it is feasible (flat 21 at pillar/sub-pillar
# level); step-wise is what makes the LARGER budgets (42) and finer
# allocations tractable without a combinatorial blow-up.
#
# Produces every combination of:
#   - budget:     21 or 42 indicators total (CONFIG$stepwise$budgets)
#   - allocation: flat (equal per group) or proportional (largest-remainder,
#                 see largest_remainder_allocation())
#   - level:      pillar (Level 3, 7 groups) or sub-pillar (Level 2, 21 groups)
#   - metric:     "distance" (country perspective - rank correlation of
#                 Euclidean distance matrices, Example/Matrix_Work.2.R) or
#                 "r_m" (variable perspective - PCA reconstruction,
#                 Example/PCA_Work.2.R)
# = 2 x 2 x 2 x 2 = 16 combinations, minus one: proportional allocation across
# sub-pillars at budget 21 coincides EXACTLY with the flat sub-pillar
# allocation (1 indicator/sub-pillar either way - see CLAUDE.md's note on
# Strategy 1b/proportional_target_matched), so re-running it would just
# re-solve the identical search problem under a different label. That leaves
# 7 allocation combinations x 2 metrics = 14 searches, each run from
# CONFIG$stepwise$n_restarts random restarts (best-scoring restart kept), as
# in Oscar_stepwise.R's `for (i in 1:5) best_n_step(3, 100)`.
#
# Unlike 03/03b's exhaustive strategies, results here are NOT re-clustered
# with SIMPROF: that comparison already exists for the flat/proportional
# budget-21/42 allocations in outputs/clustering_rand_index*.csv, and SIMPROF
# is by far the slowest step in the pipeline (see scripts/R/simprof_fast.R) -
# repeating it for 14 more subsets would dominate this script's run time for
# no new information about the search method itself.
#
# Outputs (outputs/):
#   - stepwise_search_results.csv    one row per (combination, metric): the
#                                     best restart's local score, the
#                                     assembled subset's score against the
#                                     FULL dataset (independent check, as in
#                                     03/03b/04), and its variable list.
#   - stepwise_selected_variables.csv one row per selected indicator, tagged
#                                     with its combination, metric and group -
#                                     the long-format companion to the summary
#                                     above, for report tables/plots.
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
# 1. LOAD DATA + FULL-DATASET REFERENCES
# ============================================================================
message("Loading GII data...")
gii_data <- load_gii_data()
imeta <- gii_data$imeta
indicator_data <- gii_data$indicator_data_only
n_indicators <- ncol(indicator_data)
indicator_matrix <- as.matrix(indicator_data)

message(paste("Loaded data with", nrow(indicator_data), "countries and", n_indicators, "indicators"))

full_distance_matrix <- sim_matrix(indicator_matrix)
pca_info <- get_pca_components(indicator_data, threshold = CONFIG$pca$variance_threshold)
eig_values <- pca_info$eigenvalues
data_pc <- indicator_matrix %*% as.matrix(pca_info$pca$rotation)

# Independent check of an assembled subset against the full dataset, in the
# metric it was selected under - mirrors strategyX_correlation /
# strategyX_rm in 03/03b/04, computed once per assembled subset rather than
# during the search itself.
score_against_full <- function(vars, metric) {
  if (metric == "distance") {
    spearman_dist_cor(full_distance_matrix, sim_matrix(indicator_matrix[, vars, drop = FALSE]))
  } else {
    idx <- match(vars, colnames(indicator_matrix))
    r_m(idx, data_pc, indicator_matrix, eig_values)
  }
}

# ============================================================================
# 2. BUILD THE ALLOCATION COMBINATIONS
# ============================================================================
message("\n=== BUILDING ALLOCATIONS ===")

build_combinations <- function() {
  combos <- list()
  for (budget in CONFIG$stepwise$budgets) {

    # --- Pillar level (Level 3, 7 groups) ---
    stopifnot(budget %% length(CONFIG$pillars_level3) == 0)
    flat_pillar_k <- budget %/% length(CONFIG$pillars_level3)
    combos[[paste0("flat_pillar_", budget)]] <- list(
      label = sprintf("Flat allocation by pillar (%d/pillar, %d total)", flat_pillar_k, budget),
      level = 3,
      allocation = setNames(rep(flat_pillar_k, length(CONFIG$pillars_level3)), CONFIG$pillars_level3)
    )
    combos[[paste0("prop_pillar_", budget)]] <- list(
      label = sprintf("Proportional allocation across pillars (%d total)", budget),
      level = 3,
      allocation = largest_remainder_allocation(budget, imeta, CONFIG$pillars_level3, level = 3)
    )

    # --- Sub-pillar level (Level 2, 21 groups) ---
    stopifnot(budget %% length(CONFIG$subpillars_level2) == 0)
    flat_subpillar_k <- budget %/% length(CONFIG$subpillars_level2)
    combos[[paste0("flat_subpillar_", budget)]] <- list(
      label = sprintf("Flat allocation by sub-pillar (%d/sub-pillar, %d total)", flat_subpillar_k, budget),
      level = 2,
      allocation = setNames(rep(flat_subpillar_k, length(CONFIG$subpillars_level2)), CONFIG$subpillars_level2)
    )
    if (budget == 21) {
      message("Skipping proportional allocation across sub-pillars at budget 21: ",
              "coincides exactly with the flat allocation above (1/sub-pillar either way).")
    } else {
      combos[[paste0("prop_subpillar_", budget)]] <- list(
        label = sprintf("Proportional allocation across sub-pillars (%d total)", budget),
        level = 2,
        allocation = largest_remainder_allocation(budget, imeta, CONFIG$subpillars_level2, level = 2)
      )
    }
  }
  combos
}

combinations <- build_combinations()
message(paste("Built", length(combinations), "allocation combinations:", paste(names(combinations), collapse = ", ")))

# ============================================================================
# 3. RUN STEP-WISE SEARCH FOR EVERY COMBINATION x METRIC
# ============================================================================
message("\n=== STEP-WISE SEARCH ===")
metrics <- c("distance", "r_m")
n_restarts <- CONFIG$stepwise$n_restarts
n_rounds <- CONFIG$stepwise$n_rounds

results_rows <- list()
selected_vars_rows <- list()

for (combo_name in names(combinations)) {
  combo <- combinations[[combo_name]]
  # Group label for each slot in the vars vector stepwise_search() returns -
  # its order follows combo$allocation, so this lines up positionally.
  slot_group <- rep(names(combo$allocation), combo$allocation)

  for (metric in metrics) {
    message(sprintf("\n--- %s | metric = %s ---", combo$label, metric))

    best <- NULL
    for (restart in seq_len(n_restarts)) {
      run <- stepwise_search(
        allocation = combo$allocation, imeta = imeta, data = indicator_data,
        level = combo$level, metric = metric, n_rounds = n_rounds
      )
      message(sprintf("  restart %d/%d: score = %.6f (%d round%s to converge)",
                       restart, n_restarts, run$score, run$rounds,
                       if (run$rounds == 1) "" else "s"))
      if (is.null(best) || run$score > best$score) best <- run
    }

    full_score <- score_against_full(best$vars, metric)
    budget <- sum(combo$allocation)
    level_label <- if (combo$level == 3) "pillar" else "sub-pillar"
    allocation_type <- if (startsWith(combo_name, "flat_")) "flat" else "proportional"

    message(sprintf("  BEST of %d restarts: local score = %.6f | full-dataset score = %.6f",
                     n_restarts, best$score, full_score))

    row_key <- paste(combo_name, metric, sep = "_")
    results_rows[[row_key]] <- data.frame(
      combination = combo_name,
      description = combo$label,
      level = level_label,
      allocation_type = allocation_type,
      metric = metric,
      budget = budget,
      n_restarts = n_restarts,
      best_restart_score = best$score,
      full_dataset_score = full_score,
      rounds_to_converge = best$rounds,
      variables = paste(sort(best$vars), collapse = ","),
      stringsAsFactors = FALSE
    )

    selected_vars_rows[[row_key]] <- data.frame(
      combination = combo_name,
      metric = metric,
      group = slot_group,
      variable = best$vars,
      stringsAsFactors = FALSE
    )
  }
}

results_df <- bind_rows(results_rows)
selected_vars_df <- bind_rows(selected_vars_rows)

# ============================================================================
# 4. SAVE OUTPUTS
# ============================================================================
message("\n=== SAVING OUTPUTS ===")

results_file <- file.path(CONFIG$paths$outputs_dir, "stepwise_search_results.csv")
write.csv(results_df, results_file, row.names = FALSE)
message(paste("Step-wise search summary saved to:", results_file))

selected_vars_file <- file.path(CONFIG$paths$outputs_dir, "stepwise_selected_variables.csv")
write.csv(selected_vars_df, selected_vars_file, row.names = FALSE)
message(paste("Step-wise selected variables saved to:", selected_vars_file))

message("\n=== Step-wise optimization completed. ===")
