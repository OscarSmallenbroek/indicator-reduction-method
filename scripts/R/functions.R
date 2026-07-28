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
  imeta %>%
    filter(Level == 1 & grepl(pillar_code, Parent)) %>%
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

#' Run hierarchical clustering with silhouette-based optimal k
#' Uses base R + cluster package (no clustsig required)
#' @param data Data matrix (rows = observations, cols = variables)
#' @param method.cluster Clustering method for hclust (default: "average")
#' @param method.distance Distance method for dist (default: "euclidean")
#' @param max_k Maximum number of clusters to consider (default: 15)
#' @return List with hclust result, cluster assignments, and number of groups
run_simprof <- function(data, 
                        method.cluster = "average",
                        method.distance = "euclidean",
                        max_k = 15, ...) {
  # Compute distance matrix
  d <- dist(data, method = method.distance)
  
  # Perform hierarchical clustering
  hc <- hclust(d, method = method.cluster)
  
  # Determine optimal number of clusters using silhouette width
  n_obs <- nrow(data)
  k_range <- 2:min(max_k, n_obs - 1)
  
  if (length(k_range) < 1) {
    n_clusters <- 1
  } else {
    sil_widths <- sapply(k_range, function(k) {
      cl <- cutree(hc, k = k)
      if (length(unique(cl)) < 2) return(0)
      sil <- cluster::silhouette(cl, d)
      mean(sil[, 3])
    })
    n_clusters <- k_range[which.max(sil_widths)]
  }
  
  # Cut tree at optimal number of clusters
  clusters <- cutree(hc, k = n_clusters)
  
  # Build cluster list in simprof-like format
  cluster_list <- list()
  for (i in 1:n_clusters) {
    cluster_list[[i]] <- as.character(which(clusters == i))
  }
  
  # Return object with simprof-compatible structure
  result <- list(
    hclust = hc,
    numgroups = n_clusters,
    significantclusters = cluster_list,
    cluster_assignments = clusters,
    nclusters = n_clusters,
    silhouette_widths = if (length(k_range) > 0) {
      setNames(sil_widths, k_range)
    } else {
      NA
    }
  )
  class(result) <- "simprof_compat"
  
  return(result)
}

#' Calculate Rand Index between two cluster assignments
#' Uses base R only (no partitionComparison required)
#' @param res1 First clustering result (must have cluster_assignments)
#' @param res2 Second clustering result
#' @return Rand Index value (numeric between 0 and 1)
rand_index_from_simprof <- function(res1, res2) {
  # Extract cluster assignments
  cl1 <- res1$cluster_assignments
  cl2 <- res2$cluster_assignments
  
  n <- length(cl1)
  total_pairs <- n * (n - 1) / 2
  
  # Count agreements: pairs in same cluster in both OR different in both
  agreements <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      if ((cl1[i] == cl1[j] && cl2[i] == cl2[j]) || 
          (cl1[i] != cl1[j] && cl2[i] != cl2[j])) {
        agreements <- agreements + 1
      }
    }
  }
  
  agreements / total_pairs
}

#' Calculate within-cluster equivalence from two clustering results
#' Uses base R only (no partitionComparison required)
#' @param res1 First clustering result (must have cluster_assignments)
#' @param res2 Second clustering result
#' @return Data frame with cluster equivalence measures
cluster_equivalence <- function(res1, res2) {
  cl1 <- res1$cluster_assignments
  cl2 <- res2$cluster_assignments
  
  n_clusters <- max(cl1)
  equivalence_results <- list()
  
  for (i in 1:n_clusters) {
    items_in_cluster_i <- which(cl1 == i)
    
    if (length(items_in_cluster_i) > 0) {
      # What clusters are these items in res2?
      clusters_in_res2 <- unique(cl2[items_in_cluster_i])
      
      # Calculate proportion in each res2 cluster
      prop_in_each <- sapply(clusters_in_res2, function(clust) {
        sum(cl2[items_in_cluster_i] == clust) / length(items_in_cluster_i)
      })
      
      # Maximum proportion (purity)
      max_prop <- max(prop_in_each)
      dominant_cluster <- clusters_in_res2[which.max(prop_in_each)]
      
      equivalence_results[[i]] <- data.frame(
        cluster_res1 = i,
        size = length(items_in_cluster_i),
        purity = max_prop,
        dominant_cluster_res2 = dominant_cluster,
        stringsAsFactors = FALSE
      )
    }
  }
  
  do.call(rbind, equivalence_results)
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