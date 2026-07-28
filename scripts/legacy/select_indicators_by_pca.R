#' Select Most Representative Indicators Using PCA
#' 
#' This function selects the most representative subset of indicators that create an index
#' using PCA-based selection. It selects indicators which have the strongest loading on one 
#' PC and lowest loading on all others, aiming for K=1 indicator per retained PC.
#'
#' @param data Data frame or matrix of indicators (rows = observations, columns = indicators)
#' @param target_k Number of indicators to select per PC (default = 1)
#' @param eigenvalue_threshold Threshold for Kaiser criterion (default = 1)
#' @param epsilon Small value to avoid division by zero
#' @return List containing:
#'   - selected_indicators: Indices/names of selected indicators
#'   - selected_loadings: Loadings of selected indicators
#'   - r_squared: R-squared values for each retained PC
#'   - retained_pcs: Number of retained PCs
#'   - explained_variance: Explained variance for each retained PC
#'   - loadings_matrix: Full loadings matrix
#'   - eigenvalues: Eigenvalues from PCA
#' @importFrom stats prcomp sd
#' @importFrom stats lm summary
#' @examples
#' # Generate synthetic data
#' set.seed(42)
#' n_obs <- 200
#' n_factors <- 5
#' n_indicators <- 20
#' 
#' # Create loading matrix
#' loadings_matrix <- matrix(0, nrow = n_indicators, ncol = n_factors)
#' for (i in 1:n_factors) {
#'   start_idx <- (i - 1) * 4 + 1
#'   end_idx <- i * 4
#'   loadings_matrix[start_idx:end_idx, i] <- runif(4, 0.6, 0.9)
#' }
#' 
#' # Generate latent factors
#' latent_factors <- matrix(rnorm(n_obs * n_factors), nrow = n_obs)
#' 
#' # Generate indicators
#' indicators <- latent_factors %*% t(loadings_matrix) + 
#'               matrix(rnorm(n_obs * n_indicators, sd = 0.3), nrow = n_obs)
#' 
#' # Run indicator selection
#' result <- select_indicators_by_pca(indicators)
#' print(result$selected_indicators)
#' print(result$r_squared)

select_indicators_by_pca <- function(data, target_k = 1, eigenvalue_threshold = 1) {
  # Convert to matrix if data frame
  if (is.data.frame(data)) {
    data <- as.matrix(data)
  }
  
  # Step 1: Standardize the data (mean-center and scale to unit variance)
  data_std <- scale(data)
  
  # Step 2: Perform PCA on standardized data
  pca_result <- prcomp(data_std, center = TRUE, scale. = FALSE)
  eigenvalues <- pca_result$sdev^2  # eigenvalues from PCA
  
  # Step 3: Apply Kaiser criterion to retain PCs with eigenvalues > threshold
  retained_pcs <- which(eigenvalues > eigenvalue_threshold)
  n_retained <- length(retained_pcs)
  
  if (n_retained == 0) {
    warning("No PCs retained with current eigenvalue threshold. Try lowering the threshold.")
    return(list(error = "No PCs retained"))
  }
  
  # Extract loadings matrix for retained PCs
  loadings_matrix <- pca_result$rotation[, retained_pcs, drop = FALSE]
  explained_variance <- eigenvalues[retained_pcs] / sum(eigenvalues)
  
  # Step 4: Implement indicator selection algorithm
  # For each retained PC, select indicator with highest score
  # Score = |loading_on_target_PC| / max(|loadings_on_other_PCs|)
  
  n_indicators <- nrow(loadings_matrix)
  selected_indices <- integer(n_retained)
  selected_scores <- numeric(n_retained)
  
  for (pc_idx in seq_len(n_retained)) {
    # Calculate score for each indicator
    target_loadings <- abs(loadings_matrix[, pc_idx])
    
    # For each indicator, find max loading on other PCs
    max_other_loading <- numeric(n_indicators)
    for (ind_idx in seq_len(n_indicators)) {
      other_loadings <- abs(loadings_matrix[ind_idx, -pc_idx])
      max_other_loading[ind_idx] <- max(other_loadings, na.rm = TRUE)
    }
    
    # Calculate score = |loading_on_target| / max(|loadings_on_other|)
    scores <- target_loadings / (max_other_loading)
    
    # Select indicator with highest score
    best_idx <- which.max(scores)
    selected_indices[pc_idx] <- best_idx
    selected_scores[pc_idx] <- scores[best_idx]
  }
  
  # Step 5: Calculate R-squared values
  # Regress full dataset PC scores on selected indicators
  
  # Get PC scores for retained PCs
  pc_scores <- pca_result$x[, retained_pcs, drop = FALSE]
  
  # Prepare data for regression (selected indicators from original standardized data)
  selected_data <- data_std[, selected_indices, drop = FALSE]
  
  # Calculate R² for each retained PC
  r_squared <- numeric(n_retained)
  for (pc_idx in seq_len(n_retained)) {
    # Regress PC scores on selected indicators
    model <- lm(pc_scores[, pc_idx] ~ selected_data)
    r_squared[pc_idx] <- summary(model)$r.squared
  }
  
  # Step 6: Prepare output
  # Get indicator names if available
  indicator_names <- colnames(data)
  if (is.null(indicator_names)) {
    indicator_names <- as.character(selected_indices)
  }
  
  # Build result list
  result <- list(
    selected_indicators = selected_indices,
    selected_indicator_names = indicator_names[selected_indices],
    selected_loadings = loadings_matrix[selected_indices, , drop = FALSE],
    r_squared = r_squared,
    retained_pcs = n_retained,
    retained_pc_indices = retained_pcs,
    explained_variance = explained_variance,
    loadings_matrix = loadings_matrix,
    eigenvalues = eigenvalues[retained_pcs],
    total_eigenvalues = eigenvalues,
    r_squared_by_indicator = list()
  )
  
  # Add detailed information about each selected indicator
  for (i in seq_len(n_retained)) {
    result$r_squared_by_indicator[[i]] <- list(
      indicator_index = selected_indices[i],
      indicator_name = indicator_names[selected_indices[i]],
      r_squared = r_squared[i],
      pc_index = retained_pcs[i]
    )
  }
  
  return(result)
}


#' Generate Synthetic Data for Testing
#'
#' Creates synthetic indicator data with known underlying factor structure
#' and realistic cross-loadings
#'
#' @param n_obs Number of observations
#' @param n_factors Number of underlying factors
#' @param n_indicators Total number of indicators
#' @param loading_range Range for primary loading values (min, max)
#' @param cross_loading_range Range for cross-loading values (min, max)
#' @param noise_sd Standard deviation of noise
#' @param seed Random seed for reproducibility
#' @return List containing data matrix and true loading matrix

generate_synthetic_data <- function(n_obs = 200, n_factors = 5, n_indicators = 20,
                                     loading_range = c(0.6, 0.9), cross_loading_range = c(-0.3, 0.3),
                                     noise_sd = 0.3, seed = 42) {
  set.seed(seed)

  # Create loading matrix with known structure
  # Each factor loads on approximately n_indicators/n_factors indicators
  loadings_matrix <- matrix(0, nrow = n_indicators, ncol = n_factors)

  indicators_per_factor <- floor(n_indicators / n_factors)

  # Assign primary loadings
  for (i in 1:n_factors) {
    start_idx <- (i - 1) * indicators_per_factor + 1
    end_idx <- min(i * indicators_per_factor, n_indicators)

    # Assign primary loadings in the range
    loadings_matrix[start_idx:end_idx, i] <- runif(end_idx - start_idx + 1,
                                                      loading_range[1],
                                                      loading_range[2])
  }

  # Add random cross-loadings (small loadings on other factors)
  for (i in 1:n_indicators) {
    for (j in 1:n_factors) {
      if (loadings_matrix[i, j] == 0) {
        # Add small random cross-loading
        loadings_matrix[i, j] <- runif(1, cross_loading_range[1], cross_loading_range[2])
      } else {
        # With some probability, reduce the primary loading slightly
        # to make structure less rigid
        if (runif(1) < 0.3) {
          loadings_matrix[i, j] <- loadings_matrix[i, j] * runif(1, 0.7, 0.95)
        }
      }
    }
  }

  # Ensure diagonal dominance (primary loadings should be stronger than cross-loadings)
  for (i in 1:n_indicators) {
    factor_idx <- ceiling(i / indicators_per_factor)
    if (factor_idx > n_factors) factor_idx <- n_factors
    primary_loading <- loadings_matrix[i, factor_idx]

    # Scale cross-loadings to be smaller than primary loading
    for (j in 1:n_factors) {
      if (j != factor_idx && abs(loadings_matrix[i, j]) > primary_loading * 0.5) {
        loadings_matrix[i, j] <- sign(loadings_matrix[i, j]) * primary_loading * 0.5 * runif(1, 0.3, 0.7)
      }
    }
  }

  # Generate latent factors
  latent_factors <- matrix(rnorm(n_obs * n_factors), nrow = n_obs)

  # Generate indicators as linear combination of factors + noise
  indicators <- latent_factors %*% t(loadings_matrix) +
                matrix(rnorm(n_obs * n_indicators, sd = noise_sd), nrow = n_obs)

  # Add column names
  colnames(indicators) <- paste0("Indicator_", 1:n_indicators)
  rownames(loadings_matrix) <- paste0("Indicator_", 1:n_indicators)
  colnames(loadings_matrix) <- paste0("Factor_", 1:n_factors)

  return(list(
    data = indicators,
    true_loadings = loadings_matrix,
    true_factors = latent_factors
  ))
}


#' Helper function to print summary of indicator selection results
#' @param results Output from select_indicators_by_pca()
#' @importFrom utils capture.output

print_summary <- function(results) {
  cat("=== PCA-based Indicator Selection Summary ===\n\n")
  
  cat("Retained PCs (Kaiser criterion, eigenvalue > 1):", results$retained_pcs, "\n\n")
  
  cat("Selected Indicators:\n")
  for (i in seq_along(results$selected_indicators)) {
    cat(sprintf("  PC %d: Indicator %s (index %d), R² = %.4f\n",
                results$retained_pc_indices[i],
                results$selected_indicator_names[i],
                results$selected_indicators[i],
                results$r_squared[i]))
  }
  
  cat("\nExplained Variance:\n")
  for (i in seq_along(results$explained_variance)) {
    cat(sprintf("  PC %d: %.2f%%\n", results$retained_pc_indices[i], 
                results$explained_variance[i] * 100))
  }
  
  cat("\nOverall R² (regression of PC scores on selected indicators):\n")
  cat(sprintf("  Mean R²: %.4f\n", mean(results$r_squared)))
  cat(sprintf("  Min R²: %.4f\n", min(results$r_squared)))
  cat(sprintf("  Max R²: %.4f\n", max(results$r_squared)))
  
  cat("\nSelected Indicator Loadings on Retained PCs:\n")
  print(round(results$selected_loadings, 3))
}


# Example usage (commented out, uncomment to test)
# data <- generate_synthetic_data(n_obs = 200, n_factors = 5, n_indicators = 20)
# result <- select_indicators_by_pca(data$data)
# print_summary(result)

