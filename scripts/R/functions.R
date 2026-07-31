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

#' Compute multiple correlation (R-squared) between a PC and a subset of variables
#' @param pc_idx Index of principal component
#' @param subset_idx Indices of variables in subset
#' @param data_pc Matrix of PC scores
#' @param dm Matrix of original/scaled data
#' @return R-squared value
PC_cor <- function(pc_idx, subset_idx, data_pc, dm) {
  out <- data_pc[, pc_idx]
  inputs <- dm[, subset_idx, drop = FALSE]
  matrix_for_analysis <- cbind(out, inputs)
  data <- as.data.frame(matrix_for_analysis)
  model <- lm(out ~ ., data = data)
  summary(model)$r.squared
}

#' Compute weighted r_m metric for a proposed subset of variables
#' @param subset_idx Indices of variables in subset
#' @param data_pc Matrix of PC scores
#' @param dm Matrix of original/scaled data
#' @param eig_values Eigenvalues from PCA
#' @return Weighted r_m value
r_m <- function(subset_idx, data_pc, dm, eig_values) {
  r_sum <- 0
  n_components <- min(length(eig_values), ncol(data_pc))
  
  for (i in 1:n_components) {
    r_m.i.squared <- PC_cor(i, subset_idx, data_pc, dm)
    r_sum <- r_sum + eig_values[i] * r_m.i.squared
  }
  
  eig_sum <- sum(eig_values[1:n_components])
  sqrt(r_sum / eig_sum)
}

#' Sub-PCA versions of PC_cor and r_m (for within-group analysis)
#' @param pc_idx Index of principal component
#' @param subset_idx Indices of variables in subset
#' @param data_pc_sub Matrix of PC scores from sub-dataset
#' @param dm_sub Matrix of original/scaled data from sub-dataset
#' @return R-squared value
PC_cor_sub <- function(pc_idx, subset_idx, data_pc_sub, dm_sub) {
  out <- data_pc_sub[, pc_idx]
  inputs <- dm_sub[, subset_idx, drop = FALSE]
  matrix_for_analysis <- cbind(out, inputs)
  data <- as.data.frame(matrix_for_analysis)
  model <- lm(out ~ ., data = data)
  summary(model)$r.squared
}

#' Sub-PCA version of r_m
#' @param subset_idx Indices of variables in subset
#' @param n_var_sub Number of variables in sub-dataset
#' @param data_pc_sub Matrix of PC scores from sub-dataset
#' @param dm_sub Matrix of original/scaled data from sub-dataset
#' @param eig_values_sub Eigenvalues from sub-dataset PCA
#' @return Weighted r_m value
r_m_sub <- function(subset_idx, n_var_sub, data_pc_sub, dm_sub, eig_values_sub) {
  r_sum <- 0
  
  for (i in 1:n_var_sub) {
    r_m.i.squared <- PC_cor_sub(i, subset_idx, data_pc_sub, dm_sub)
    r_sum <- r_sum + eig_values_sub[i] * r_m.i.squared
  }
  
  eig_sum <- sum(eig_values_sub)
  sqrt(r_sum / eig_sum)
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

#' Get indicators under a specific pillar
#' @param pillar_code Pillar code (e.g., "IN.1")
#' @param imeta Metadata data frame
#' @return Vector of indicator codes
pillar_indicators <- function(pillar_code, imeta) {
  # NUM (e.g. "IN.1.1.1") embeds the pillar code; Parent only holds the
  # direct sub-pillar parent, so it can't be used to match pillar level.
  imeta %>%
    filter(Level == 1 & grepl(pillar_code, NUM, fixed = TRUE)) %>%
    pull(iCode)
}

#' Get indicators under a specific sub-pillar
#' @param subpillar_code Sub-pillar code (e.g., "SP1.1")
#' @param imeta Metadata data frame
#' @return Vector of indicator codes
subpillar_indicators <- function(subpillar_code, imeta) {
  imeta %>%
    filter(Level == 1 & grepl(subpillar_code, Parent)) %>%
    pull(iCode)
}

#' Exhaustive search for best k variables within a group
#' @param k Number of variables to select
#' @param group_codes Vector of group codes ( grid column values)
#' @param imeta Metadata data frame
#' @param data Standardized indicator data
#' @param level Level to group by (3 for pillar, 2 for sub-pillar)
#' @return Data frame with best variables and their r_m score
best_k_within_group <- function(k, group_codes, imeta, data, level = 3) {
  # Determine the column to filter based on level
  if (level == 3) {
    filter_col <- "NUM"  # Pillar level (e.g., "IN.1")
    group_prefix_length <- 3  # First 3 chars like "IN."
  } else if (level == 2) {
    filter_col <- "Parent"  # Sub-pillar level
    group_prefix_length <- nchar(group_codes[1])  # Full sub-pillar code length
  } else {
    stop("Level must be 2 (sub-pillar) or 3 (pillar)")
  }
  
  best_results <- list()
  
  for (group_code in group_codes) {
    # Select variables under this group
    if (level == 3) {
      variables <- imeta %>%
        filter(grepl(group_code, !!sym(filter_col)) & Level == 1) %>%
        pull(iCode)
    } else {  # level == 2
      variables <- imeta %>%
        filter(!!sym(filter_col) == group_code & Level == 1) %>%
        pull(iCode)
    }
    
    nc_var <- length(variables)
    
    if (nc_var > k) {
      # Need to search
      dm_sub <- as.matrix(data[, variables, drop = FALSE])
      n_var_sub <- min(nc_var, 27)  # Limit as in original code
      pca_result_sub <- prcomp(dm_sub, center = TRUE, scale. = TRUE)
      eig_values_sub <- pca_result_sub$sdev^2
      pcm_sub <- as.matrix(pca_result_sub$rotation)
      data_pc_sub <- dm_sub %*% pcm_sub
      
      combos <- combn(1:nc_var, k, simplify = FALSE)
      best_rm <- 0
      best_vars <- NULL
      
      for (i in seq_along(combos)) {
        combo <- combos[[i]]
        s <- r_m_sub(combo, n_var_sub, data_pc_sub, dm_sub, eig_values_sub)
        if (s > best_rm) {
          best_rm <- s
          best_vars <- variables[combo]
        }
      }
      
      best_results[[group_code]] <- data.frame(
        group = group_code,
        best_rm = best_rm,
        variables = paste(sort(best_vars), collapse = ","),
        stringsAsFactors = FALSE
      )
    } else {
      # Take all variables if not enough
      best_results[[group_code]] <- data.frame(
        group = group_code,
        best_rm = 1.0,  # Perfect score when taking all
        variables = paste(sort(variables), collapse = ","),
        stringsAsFactors = FALSE
      )
    }
  }
  
  # Combine results
  result_df <- bind_rows(best_results)
  return(result_df)
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
#' Unlike proportional_allocation(), no group is guaranteed a minimum of 1 -
#' a group with a very small share of the total can receive zero, and groups
#' with the largest share get first claim on leftover budget.
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
#' @param allocation Named vector of group_code -> k (number of variables to select)
#' @param imeta Metadata data frame
#' @param data Standardized indicator data
#' @param level Level to group by (3 for pillar, 2 for sub-pillar)
#' @return Data frame with best variables and their r_m score per group
best_allocation_within_group <- function(allocation, imeta, data, level = 3) {
  best_results <- list()

  for (group_code in names(allocation)) {
    k <- allocation[[group_code]]
    if (k == 0) next

    variables <- group_variables(group_code, imeta, level)
    nc_var <- length(variables)

    if (nc_var > k) {
      # Need to search
      dm_sub <- as.matrix(data[, variables, drop = FALSE])
      n_var_sub <- min(nc_var, 27)  # Limit as in original code
      pca_result_sub <- prcomp(dm_sub, center = TRUE, scale. = TRUE)
      eig_values_sub <- pca_result_sub$sdev^2
      pcm_sub <- as.matrix(pca_result_sub$rotation)
      data_pc_sub <- dm_sub %*% pcm_sub

      combos <- combn(1:nc_var, k, simplify = FALSE)
      best_rm <- 0
      best_vars <- NULL

      for (i in seq_along(combos)) {
        combo <- combos[[i]]
        s <- r_m_sub(combo, n_var_sub, data_pc_sub, dm_sub, eig_values_sub)
        if (s > best_rm) {
          best_rm <- s
          best_vars <- variables[combo]
        }
      }

      best_results[[group_code]] <- data.frame(
        group = group_code,
        allocated = k,
        best_rm = best_rm,
        variables = paste(sort(best_vars), collapse = ","),
        stringsAsFactors = FALSE
      )
    } else {
      # Take all variables if the group has <= k indicators
      best_results[[group_code]] <- data.frame(
        group = group_code,
        allocated = nc_var,
        best_rm = 1.0,  # Perfect score when taking all
        variables = paste(sort(variables), collapse = ","),
        stringsAsFactors = FALSE
      )
    }
  }

  bind_rows(best_results)
}

# use simprof from the clustig package


#' Calculate Rand Index between two cluster assignments
## Define a function to calculate the Rand Index for two clusters
## Uses the results from simprof()
## Cluster format needs adjusting

rand_ind <- function(cluster_results1,cluster_results2) {
  nclust1 <- cluster_results1$numgroups
  cluster1 <- cluster_results1$significantclusters
  nclust2 <- cluster_results2$numgroups
  cluster2 <- cluster_results2$significantclusters
  cluster_array1 <- rep(0,27)
  cluster_array2 <- rep(0,27)
  for (i in 1:nclust1) {
    for (j in 1:length(cluster1[[i]])){
      sample <- as.integer(cluster1[[i]][j])
      cluster_array1[sample] <- i
    }
  }
  for (i in 1:nclust2) {
    for (j in 1:length(cluster2[[i]])){
      sample <- as.integer(cluster2[[i]][j])
      cluster_array2[sample] <- i
    }
  }
  # Register the measures to take ANY input (no clue)
  registerPartitionVectorSignatures(environment())
  # Compare the clusters without EU27 (11th item)
  return(partitionComparison::randIndex(cluster_array1[-11], cluster_array2[-11]))
}


#' Proportional allocation of samples across groups
#' @param target Target total sample size
#' @param imeta Metadata data frame
#' @param groups Vector of group codes (pillar or sub-pillar codes)
#' @param min_per Minimum number to allocate per group
#' @return Named vector of allocations
proportional_allocation <- function(target, imeta, groups, min_per = 1) {
  # Count indicators per group
  group_counts <- sapply(groups, function(group) {
    sum(imeta$Level == 1 & grepl(group, imeta$Parent))
  })
  
  # Start with proportional allocation
  allocation <- pmax(min_per, round(group_counts / sum(group_counts) * target))
  
  # Adjust to exactly match target
  current_sum <- sum(allocation)
  
  if (current_sum < target) {
    # Need to add more - add to groups with most indicators first
    while (sum(allocation) < target) {
      # Order groups by remaining capacity (descending)
      remaining_capacity <- group_counts - allocation
      ordered_groups <- groups[order(remaining_capacity, decreasing = TRUE)]
      
      for (group in ordered_groups) {
        if (allocation[group] < group_counts[group]) {
          allocation[group] <- allocation[group] + 1
          if (sum(allocation) >= target) break
        }
      }
    }
  } else if (current_sum > target) {
    # Need to remove some - remove from groups with most overallocation first
    while (sum(allocation) > target) {
      overallocation <- allocation - group_counts
      overallocation[overallocation < 0] <- 0  # Don't go below actual count
      
      # Only consider groups that can still give up items (above min_per)
      reducible_groups <- names(allocation)[allocation > min_per & overallocation > 0]
      
      if (length(reducible_groups) == 0) {
        # If we can't reduce without going below min_per, break
        break
      }
      
      # Order by amount we can reduce (descending)
      reducible_groups <- reducible_groups[order(
        allocation[reducible_groups] - min_per, 
        decreasing = TRUE
      )]
      
      for (group in reducible_groups) {
        if (allocation[group] > min_per) {
          allocation[group] <- allocation[group] - 1
          if (sum(allocation) <= target) break
        }
      }
    }
  }
  
  allocation
}