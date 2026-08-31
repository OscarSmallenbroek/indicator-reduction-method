# Shared Functions for GII Index Reduction Analysis
# Consolidated from multiple scripts to avoid duplication

#' Compute similarity/distance matrix
#' @param df Data frame or matrix
#' @param method Distance method (default: "euclidean")
#' @return Distance matrix
sim_matrix <- function(df, method = "euclidean") {
  dist(df, method = method)
}

#' Compute Spearman correlation of distance vectors
#' @param full_dist Full distance matrix (as dist object)
#' @param sub_dist Subset distance matrix (as dist object)
#' @return Spearman correlation coefficient
spearman_dist_cor <- function(full_dist, sub_dist) {
  cor(as.vector(full_dist), as.vector(sub_dist), method = "spearman")
}

#' Build a scorer for the distance ("country perspective") selection metric.
#'
#' This is the criterion from Example/Matrix_Work.2.R (`subset_test` inside
#' `find_best_c`/`find_best_d`): score a candidate subset by the Spearman
#' correlation between the *group's own* Euclidean distance matrix and the
#' candidate's. The comparison is local to the group, exactly as in the
#' original - the assembled subset is only scored against the full index once,
#' by the calling script.
#'
#' The original recomputed dist(group) inside the combination loop. Spearman is
#' Pearson on ranks and the group's distance ranks never change, so they are
#' ranked once here and reused; the result is identical.
#' @param group_matrix Numeric matrix of one group's indicators
#' @return Function taking column indices and returning the Spearman score
make_dist_scorer <- function(group_matrix) {
  full_ranks <- rank(as.vector(sim_matrix(group_matrix)))
  function(subset_idx) {
    cor(full_ranks,
        rank(as.vector(sim_matrix(group_matrix[, subset_idx, drop = FALSE]))))
  }
}

#' Build a scorer for the r_m ("variable perspective") selection metric.
#'
#' This is the criterion from Example/PCA_Work.2.R (`best_k_comp`): score a
#' candidate subset by how well it reconstructs the group's own sub-PCA.
#' @param group_matrix Numeric matrix of one group's indicators
#' @return Function taking column indices and returning the r_m score
make_rm_scorer <- function(group_matrix) {
  pca_result_sub <- prcomp(group_matrix, center = TRUE, scale. = TRUE)
  eig_values_sub <- pca_result_sub$sdev^2
  # Scores as group_matrix %*% rotation rather than pca_result_sub$x: the two
  # are equivalent here (the data is already z-scored), but exact ties do occur -
  # a 2-indicator sub-pillar scores identically whichever one is picked - so the
  # arithmetic is kept bit-for-bit as it was to keep selections stable.
  data_pc_sub <- group_matrix %*% as.matrix(pca_result_sub$rotation)
  function(subset_idx) {
    r_m(subset_idx, data_pc_sub, group_matrix, eig_values_sub)
  }
}

#' Compute multiple correlation (R-squared) between a PC and a subset of variables
#' @param pc_idx Index of principal component
#' @param subset_idx Indices of variables in subset
#' @param data_pc Matrix of PC scores
#' @param dm Matrix of original/scaled data
#' @return R-squared value
PC_cor <- function(pc_idx, subset_idx, data_pc, dm) {
  subset_r_squared(subset_idx, data_pc[, pc_idx, drop = FALSE], dm)
}

#' R-squared of every column of `targets` regressed on the same subset of `dm`.
#' The design matrix is identical for every target, so one QR factorisation
#' (via .lm.fit) serves them all - far cheaper than a formula-lm per target.
#' @param subset_idx Indices of variables in subset
#' @param targets Matrix whose columns are regressed on the subset
#' @param dm Matrix of original/scaled data
#' @return Numeric vector of R-squared values, one per column of `targets`
subset_r_squared <- function(subset_idx, targets, dm) {
  X <- cbind(1, dm[, subset_idx, drop = FALSE])
  rss <- colSums(.lm.fit(X, targets)$residuals^2)
  tss <- colSums(sweep(targets, 2, colMeans(targets))^2)
  1 - rss / tss
}

#' Compute weighted r_m metric for a proposed subset of variables
#'
#' Used both against the full-dataset PCA and against a single group's sub-PCA
#' (the "local" r_m proxy) - the two are the same calculation, so there is one
#' implementation. The number of components scored is taken from the inputs.
#' @param subset_idx Indices of variables in subset
#' @param data_pc Matrix of PC scores
#' @param dm Matrix of original/scaled data
#' @param eig_values Eigenvalues from PCA
#' @return Weighted r_m value
r_m <- function(subset_idx, data_pc, dm, eig_values) {
  n_components <- min(length(eig_values), ncol(data_pc))
  weights <- eig_values[seq_len(n_components)]
  r_squared <- subset_r_squared(
    subset_idx, data_pc[, seq_len(n_components), drop = FALSE], dm
  )
  sqrt(sum(weights * r_squared) / sum(weights))
}

#' Get PCA components and number needed to reach variance threshold
#' @param data Data matrix (scaled)
#' @param threshold Variance threshold (default: 0.95)
#' @return List with PCA results and number of components
get_pca_components <- function(data, threshold = 0.95) {
  pca_result <- prcomp(data, center = TRUE, scale. = TRUE)
  eig_values <- pca_result$sdev^2
  cumulative_variance <- cumsum(eig_values) / sum(eig_values)
  n_components <- which(cumulative_variance >= threshold)[1]
  
  list(
    pca = pca_result,
    eigenvalues = eig_values,
    n_components = n_components,
    cumulative_variance = cumulative_variance
  )
}

#' Exhaustive search for the best k variables within every group.
#' A flat k is just an allocation that gives every group the same budget, so
#' this delegates to best_allocation_within_group() rather than repeating the
#' search. The `allocated` column is dropped to keep the flat-strategy CSV
#' schema (group, <score>, variables) that the report expects, where <score> is
#' `best_rm` or `best_spearman` depending on `metric`.
#' @param k Number of variables to select per group
#' @param group_codes Vector of group codes
#' @param imeta Metadata data frame
#' @param data Standardized indicator data
#' @param level Level to group by (3 for pillar, 2 for sub-pillar)
#' @param metric Selection criterion: "r_m" (default) or "distance"
#' @return Data frame with best variables and their score
best_k_within_group <- function(k, group_codes, imeta, data, level = 3,
                                metric = c("r_m", "distance")) {
  allocation <- setNames(rep(k, length(group_codes)), group_codes)
  result <- best_allocation_within_group(allocation, imeta, data, level,
                                         metric = match.arg(metric))
  result[, setdiff(names(result), "allocated"), drop = FALSE]
}

#' Get the Level-1 indicator codes belonging to a pillar or sub-pillar group
#' @param group_code Group code (e.g. "IN.1" for level 3, "SP1.1" for level 2)
#' @param imeta Metadata data frame
#' @param level Level to group by (3 for pillar, 2 for sub-pillar)
#' @return Vector of indicator codes
group_variables <- function(group_code, imeta, level = 3) {
  if (level == 3) {
    # NUM (e.g. "IN.1.1.1") embeds the pillar code; Parent only holds the
    # direct sub-pillar parent, so it can't be used to match pillar level.
    imeta %>%
      filter(Level == 1 & grepl(group_code, NUM, fixed = TRUE)) %>%
      pull(iCode)
  } else if (level == 2) {
    imeta %>%
      filter(Level == 1 & Parent == group_code) %>%
      pull(iCode)
  } else {
    stop("Level must be 2 (sub-pillar) or 3 (pillar)")
  }
}

#' Count Level-1 indicators per group
#' @param group_codes Vector of group codes
#' @param imeta Metadata data frame
#' @param level Level to group by (3 for pillar, 2 for sub-pillar)
#' @return Named integer vector of counts, one per group_code
group_indicator_counts <- function(group_codes, imeta, level = 3) {
  counts <- sapply(group_codes, function(g) length(group_variables(g, imeta, level)))
  names(counts) <- group_codes
  counts
}

#' Largest-remainder (Hamilton) apportionment of a target total across groups,
#' proportional to each group's number of Level-1 indicators.
#' No group is guaranteed a minimum of 1: a group with a very small share of
#' the total can receive zero, and groups with the largest fractional remainder
#' get first claim on the leftover budget.
#' @param target Target total number of indicators to allocate
#' @param imeta Metadata data frame
#' @param group_codes Vector of group codes (pillar or sub-pillar codes)
#' @param level Level to group by (3 for pillar, 2 for sub-pillar)
#' @return Named integer vector of allocations, one per group_code (may include 0s)
largest_remainder_allocation <- function(target, imeta, group_codes, level = 3) {
  counts <- group_indicator_counts(group_codes, imeta, level)
  total <- sum(counts)

  if (target > total) {
    stop("Target (", target, ") exceeds total available indicators (", total, ")")
  }

  share <- counts / total * target
  allocation <- floor(share)
  remainder <- target - sum(allocation)

  if (remainder > 0) {
    frac <- share - allocation
    # Groups with the largest fractional remainder get the leftover budget first.
    # Since remainder < length(group_codes) and floor(share) < counts whenever
    # frac > 0, this never allocates more than a group's indicator count.
    top_groups <- names(counts)[order(frac, decreasing = TRUE)][seq_len(remainder)]
    allocation[top_groups] <- allocation[top_groups] + 1
  }

  allocation
}

#' Exhaustive search for the best subset within each group, using a
#' per-group allocation (e.g. from largest_remainder_allocation()) instead
#' of a single k applied uniformly to every group. Groups allocated 0 are skipped.
#'
#' The allocation rule (flat vs proportional) and the selection metric are
#' independent choices: any allocation can be searched under either metric.
#'   metric = "r_m"      - variable perspective, after Example/PCA_Work.2.R
#'                         (`best_k_comp`). Score column is `best_rm`.
#'   metric = "distance" - country perspective, after Example/Matrix_Work.2.R
#'                         (`find_best_c`/`find_best_d`). Score column is
#'                         `best_spearman`.
#' Both score locally, against the group's own sub-PCA / own distance matrix.
#' A group with no more indicators than its budget takes all of them and scores
#' 1.0 under either metric, as in both originals.
#' @param allocation Named vector of group_code -> k (number of variables to select)
#' @param imeta Metadata data frame
#' @param data Standardized indicator data
#' @param level Level to group by (3 for pillar, 2 for sub-pillar)
#' @param metric Selection criterion: "r_m" (default) or "distance"
#' @return Data frame with best variables and their score per group
best_allocation_within_group <- function(allocation, imeta, data, level = 3,
                                         metric = c("r_m", "distance")) {
  metric <- match.arg(metric)
  score_col <- if (metric == "r_m") "best_rm" else "best_spearman"
  make_scorer <- if (metric == "r_m") make_rm_scorer else make_dist_scorer

  best_results <- list()

  for (group_code in names(allocation)) {
    k <- allocation[[group_code]]
    if (k == 0) next

    variables <- group_variables(group_code, imeta, level)
    nc_var <- length(variables)

    if (nc_var > k) {
      # Need to search
      dm_sub <- as.matrix(data[, variables, drop = FALSE])
      score <- make_scorer(dm_sub)

      combos <- combn(1:nc_var, k, simplify = FALSE)
      best_score <- 0
      best_vars <- NULL

      for (i in seq_along(combos)) {
        combo <- combos[[i]]
        s <- score(combo)
        if (s > best_score) {
          best_score <- s
          best_vars <- variables[combo]
        }
      }

      allocated <- k
    } else {
      # Take all variables if the group has <= k indicators
      best_score <- 1.0  # Perfect score when taking all
      best_vars <- variables
      allocated <- nc_var
    }

    row <- data.frame(
      group = group_code,
      allocated = allocated,
      score = best_score,
      variables = paste(sort(best_vars), collapse = ","),
      stringsAsFactors = FALSE
    )
    names(row)[names(row) == "score"] <- score_col
    best_results[[group_code]] <- row
  }

  bind_rows(best_results)
}

#' Build a scorer that evaluates a candidate subset of indicator codes
#' against the WHOLE dataset, as opposed to make_dist_scorer()/
#' make_rm_scorer(), which each score against a single group's own
#' sub-matrix/sub-PCA. This is the criterion stepwise_search() optimizes: a
#' candidate combination is always compared to the full index, matching
#' subset_test() in replication/Oscar_stepwise.R (which always compares to
#' sm_broad_ranks) and the r_m_named() pattern in 04_variable_perspective.R.
#' @param data Full standardized indicator data (data frame or matrix,
#'   columns named by indicator code)
#' @param metric Selection criterion: "r_m" (default) or "distance"
#' @return Function taking a character vector of indicator codes and
#'   returning the score of that subset against the full dataset
make_full_scorer <- function(data, metric = c("r_m", "distance")) {
  metric <- match.arg(metric)
  dm <- as.matrix(data)
  if (metric == "distance") {
    full_ranks <- rank(as.vector(sim_matrix(dm)))
    function(subset_codes) {
      cor(full_ranks,
          rank(as.vector(sim_matrix(dm[, subset_codes, drop = FALSE]))))
    }
  } else {
    pca_result <- prcomp(dm, center = TRUE, scale. = TRUE)
    eig_values <- pca_result$sdev^2
    data_pc <- dm %*% as.matrix(pca_result$rotation)
    function(subset_codes) {
      r_m(subset_codes, data_pc, dm, eig_values)
    }
  }
}

#' Step-wise local-search optimization for the best subset of indicators
#' within each group, scored against the WHOLE dataset rather than locally
#' per group - contrast with best_allocation_within_group(), which is
#' exhaustive but local (each group is scored only against its own
#' sub-matrix/sub-PCA, per the CLAUDE.md "how the search works" note).
#'
#' Replicates best_n_step() in replication/Oscar_stepwise.R: start from a
#' random per-group selection, then repeatedly try swapping each selected
#' indicator for every unselected indicator from its own group, keeping any
#' swap that improves the FULL-dataset score, until a full pass over every
#' group makes no improving swap or `n_rounds` is reached. A swap is never
#' allowed to duplicate an indicator already in the selection - a duplicated
#' column would be double-weighted in the Euclidean distance / PCA and let
#' the climber buy score with fewer than the intended number of distinct
#' indicators.
#'
#' Unlike the exhaustive search, this scales to allocations too large to
#' enumerate by combination (e.g. proportional budgets), at the cost of only
#' finding a local optimum - the original was run from several random starts
#' for this reason (`for (i in 1:5) best_n_step(3, 100)`).
#'
#' @param allocation Named vector of group_code -> number of indicators to
#'   select from that group (e.g. from largest_remainder_allocation(), or
#'   `setNames(rep(n, length(group_codes)), group_codes)` for a flat n - see
#'   stepwise_k_within_group())
#' @param imeta Metadata data frame
#' @param data Standardized indicator data (data frame or matrix, columns
#'   named by indicator code)
#' @param level Level to group by (3 = pillar, 2 = sub-pillar)
#' @param metric Selection criterion: "r_m" (default) or "distance"
#' @param n_rounds Maximum number of passes over all groups (default 100,
#'   matching Oscar_stepwise.R)
#' @param init_vars Optional starting selection: a named list/vector of
#'   group_code -> character vector of indicator codes (one per allocated
#'   slot). If NULL (default), the starting selection is drawn at random from
#'   each group with sample() - set a seed before calling for reproducibility.
#' @param verbose Print the score after each round (default FALSE)
#' @return List with `vars` (character vector of selected indicator codes,
#'   ordered by group as in `allocation`), `score` (final full-dataset score)
#'   and `rounds` (number of rounds actually run before convergence)
stepwise_search <- function(allocation, imeta, data, level = 3,
                             metric = c("r_m", "distance"), n_rounds = 100,
                             init_vars = NULL, verbose = FALSE) {
  metric <- match.arg(metric)
  group_codes <- names(allocation)
  score_fn <- make_full_scorer(data, metric)

  group_vars <- setNames(
    lapply(group_codes, group_variables, imeta = imeta, level = level),
    group_codes
  )

  if (is.null(init_vars)) {
    curt_vars <- unlist(
      Map(function(g, k) sample(group_vars[[g]], k), group_codes, allocation),
      use.names = FALSE
    )
  } else {
    curt_vars <- unlist(init_vars[group_codes], use.names = FALSE)
  }
  stopifnot(length(curt_vars) == sum(allocation), !anyDuplicated(curt_vars))

  # Which group each slot in curt_vars belongs to, so a swap only ever draws
  # a replacement from that slot's own group.
  slot_group <- rep(group_codes, allocation)

  to_beat <- score_fn(curt_vars)
  if (verbose) message("Round 0: score = ", round(to_beat, 6))

  rounds_run <- 0
  for (r in seq_len(n_rounds)) {
    rounds_run <- r
    improved <- FALSE
    for (g in group_codes) {
      pot_vars <- group_vars[[g]]
      for (slot in which(slot_group == g)) {
        for (t_v in pot_vars) {
          # Re-read the slot's current occupant every candidate, not once per
          # slot: an accepted swap earlier in this same inner loop changes
          # it, and later candidates must be tested against the new
          # occupant, exactly as in Oscar_stepwise.R's best_n_step().
          c_v <- curt_vars[slot]
          if (t_v == c_v || t_v %in% curt_vars) next
          test_vars <- curt_vars
          test_vars[slot] <- t_v
          s_test <- score_fn(test_vars)
          if (s_test > to_beat) {
            curt_vars <- test_vars
            to_beat <- s_test
            improved <- TRUE
          }
        }
      }
    }
    if (verbose) message("Round ", r, ": score = ", round(to_beat, 6))
    if (!improved) break
  }

  list(vars = curt_vars, score = to_beat, rounds = rounds_run)
}

#' Step-wise search with the same number of indicators n per group - the flat
#' wrapper around stepwise_search(), mirroring best_k_within_group()'s
#' relationship to best_allocation_within_group().
#' @param n Number of indicators to select per group
#' @param group_codes Vector of group codes
#' @inheritParams stepwise_search
#' @return See stepwise_search()
stepwise_k_within_group <- function(n, group_codes, imeta, data, level = 3,
                                     metric = c("r_m", "distance"),
                                     n_rounds = 100, init_vars = NULL,
                                     verbose = FALSE) {
  allocation <- setNames(rep(n, length(group_codes)), group_codes)
  stepwise_search(allocation, imeta, data, level, match.arg(metric),
                   n_rounds, init_vars, verbose)
}

# Clustering comparison. simprof() comes from the clustsig package, which the
# calling scripts attach; partitionComparison must also be attached for
# rand_ind() to dispatch randIndex().

#' Flatten a simprof() result into a partition vector.
#' simprof reports its groups as a list of character vectors of unit indices;
#' the Rand Index needs one cluster label per unit instead.
#' @param cluster_results A simprof() result
#' @return Integer vector of cluster labels, one per unit
simprof_partition <- function(cluster_results) {
  clusters <- cluster_results$significantclusters
  # Size from the clustering itself rather than a hardcoded unit count, so this
  # works for any dataset, and assert that every unit really got a label.
  n_units <- max(as.integer(unlist(clusters)))
  partition <- integer(n_units)
  for (i in seq_len(cluster_results$numgroups)) {
    partition[as.integer(clusters[[i]])] <- i
  }
  if (any(partition == 0)) {
    stop("simprof result leaves ", sum(partition == 0), " unit(s) unassigned")
  }
  partition
}

#' Calculate the Rand Index between two simprof() clusterings
#'
#' @param cluster_results1 First simprof() result
#' @param cluster_results2 Second simprof() result
#' @param exclude_units Indices of units to drop before comparing
#' @return Rand Index
rand_ind <- function(cluster_results1, cluster_results2) {
  partition1 <- simprof_partition(cluster_results1)
  partition2 <- simprof_partition(cluster_results2)
  if (length(partition1) != length(partition2)) {
    stop("Clusterings cover different numbers of units: ",
         length(partition1), " vs ", length(partition2))
  }
  # Register the measures to accept plain integer vectors
  registerPartitionVectorSignatures(environment())
  partitionComparison::randIndex(partition1, partition2)
}
