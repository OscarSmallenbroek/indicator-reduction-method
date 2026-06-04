library(dplyr)
library(stats)

print("Loading prepared GII data for variable perspective analysis...")

# Load standardized data
data_standardized <- read.csv("data/gii_data_standardized.csv", stringsAsFactors = FALSE)
print(paste("Loaded data with", nrow(data_standardized), "countries and", ncol(data_standardized)-1, "indicators"))

# Prepare data for analysis (exclude Country column)
country_names <- data_standardized$Country
indicator_data <- data_standardized[, -1]  # Remove Country column
variable_names <- colnames(indicator_data)

# Load the metadata to identify pillars
imeta <- read.csv("gii-data/imeta.csv", stringsAsFactors = FALSE)

# Extract indicators (rows where Level == 1)
indicators <- imeta[imeta$Level == 1, ]

# Create a mapping from iCode to NUM
indicator_num_map <- indicators$NUM
names(indicator_num_map) <- indicators$iCode

# Get the indices of indicators in the data that correspond to each pillar (IN.1 to IN.5 and OUT.6 to OUT.7)
pillar_indicators <- list()
for (pillar in c("IN.1", "IN.2", "IN.3", "IN.4", "IN.5", "OUT.6", "OUT.7")) {
  # Find all indicators that belong to this pillar (NUM starts with the pillar code)
  pillar_indicator_codes <- indicators[grep(paste0("^", pillar), indicators$NUM), "iCode"]
  
  # Get the corresponding column indices in our data
  pillar_indices <- which(variable_names %in% pillar_indicator_codes)
  pillar_indicators[[pillar]] <- pillar_indices
}

print("Performing PCA analysis...")

# Perform PCA on the full dataset
pca_result <- prcomp(indicator_data, center = TRUE, scale. = TRUE)
summary_pca <- summary(pca_result)

# Eigenvalues in decreasing order
eig_values <- pca_result$sdev^2

print("PCA Results:")
print(paste("Number of principal components:", length(eig_values)))
print("Proportion of variance explained by first 10 components:")
print(round(summary_pca$importance[2, 1:min(10, length(eig_values))], 4))

# Calculate cumulative variance explained
cumulative_variance <- cumsum(summary_pca$importance[2, ])
components_80_pct <- which(cumulative_variance >= 0.8)[1]
components_90_pct <- which(cumulative_variance >= 0.9)[1]
components_95_pct <- which(cumulative_variance >= 0.95)[1]

print(paste("Components needed for 80% variance:", components_80_pct))
print(paste("Components needed for 90% variance:", components_90_pct))
print(paste("Components needed for 95% variance:", components_95_pct))

# Function to calculate weighted r_m for a proposed subset of variables
# This implements Equation 1 from the manuscript:
# r_m = sqrt((sum(i=1 to p)(λ_i * (r_m)_i^2)) / (sum(j=1 to p)λ_j))

# First, we need matrices for multiplication
dm <- as.matrix(indicator_data)
pcm <- as.matrix(pca_result$rotation)
data_pc <- dm %*% pcm

# Function to find the multiple correlation between a PC and a subset of variables
PC_cor <- function(PC, subset_indices, data_pc_matrix, dm_matrix) {
  out <- data_pc_matrix[, PC]
  inputs <- dm_matrix[, subset_indices]
  matrix_for_analysis <- cbind(out, inputs)
  data <- as.data.frame(matrix_for_analysis)
  model <- lm(out ~ inputs, data = data)
  model_summary <- summary(model)
  mult_r.squared <- model_summary$r.squared
  return(mult_r.squared)
}


# Create a data frame with pillar information
pillar_df <- data.frame(
  pillar_name = names(pillar_indicators),
  n_indicators = sapply(pillar_indicators, length)
)
# Calculate percentage of total indicators for each pillar
pillar_df$percent_of_total <- pillar_df$n_indicators / sum(pillar_df$n_indicators) * 100

print("Pillar information:")
print(pillar_df)

# Function to calculate allocation for a given target size
calculate_allocation <- function(target_size) {

  if (target_size < 7){
    stop('Target size less than 7. Not able to allocate at least one indicator for each pillar')
  }
  # Start with proportional allocation, ensuring at least 1 from each pillar
  allocation <- pmax(1, round(pillar_df$percent_of_total/100 * target_size))

  # Adjust to exactly match target size
  while (sum(allocation) < target_size) {
    # Add one to the pillar with largest number of indicators first
    ordered_pillars <- pillar_df$pillar_name[order(pillar_df$n_indicators, decreasing = TRUE)]
    for (p in ordered_pillars) {
      p_idx <- which(pillar_df$pillar_name == p)
      if (allocation[p_idx] < pillar_df$n_indicators[p_idx]) {
        allocation[p_idx] <- allocation[p_idx] + 1
        break
      }
    }
    if (sum(allocation) >= target_size) break
  }
  names(allocation)<-pillar_df$pillar_name
  return(allocation)
}
calculate_allocation(20)
    
# Add the allocation to the pillar data frame
# This part will be used in the sampling loop
pillar_df$allocation <- 0  # Initialize

# Corrected, fully implemented r_m function
r_m <- function(subset_indices, data_pc_matrix, dm_matrix, eig_values_vector) {
  r_sum <- 0
  n_components <- min(length(eig_values_vector), ncol(data_pc_matrix))

  for (i in 1:n_components) {
    r_m.i.squared <- PC_cor(i, subset_indices, data_pc_matrix, dm_matrix)
    r_sum <- r_sum + eig_values_vector[i] * r_m.i.squared 
  }
  
  eig_sum <- sum(eig_values_vector[1:n_components])
  return(sqrt(r_sum / eig_sum))
}

print("Testing with random subsets...")

set.seed(456)  # For reproducibility
subset_sizes <- c(10, 15, 20, 25, 30,35,40)
r_m_results <- data.frame(
  size = subset_sizes,
  mean_r_m = numeric(length(subset_sizes)),
  max_r_m = numeric(length(subset_sizes))
)

for (i in 1:length(subset_sizes)) {
  size <- subset_sizes[i]
  print(paste("Testing subset size:", size))
  
  # Test 20 random subsets of this size (reduced from 50 for faster execution)
  r_m_values <- numeric(20)
  
  for (j in 1:20) {
    selected_indices <- c()

    # 1. Compute proportional allocation: sample from each pillar by size
    allocation<-calculate_allocation(size)
    # 2. Sample from each pillar according to allocation
    for (pillar in names(allocation)) {
      n_to_sample <- allocation[pillar]
      pillar_idx <- pillar_indicators[[pillar]]

      # Handle case where pillar has only one indicator or we need only one
      if (n_to_sample == 1 || length(pillar_idx) == 1) {
        # Just sample one if only needed
        selected_pillar_ind <- sample(pillar_idx, 1)
        selected_indices <- c(selected_indices, selected_pillar_ind)
      } else {
        # Use correlation-based clustering within the pillar
        sub_data <- indicator_data[, pillar_idx, drop = FALSE]
        if (nrow(sub_data) < 2 || ncol(sub_data) < 2) {
          # If not enough data, sample randomly
          selected_indices <- c(selected_indices, sample(pillar_idx, n_to_sample))
        } else {
          # Compute correlation matrix and cluster
          cor_mat <- cor(sub_data)
          # Replace NAs (if any) with 0 for clustering
          cor_mat[is.na(cor_mat)] <- 0
          d <- as.dist((1 - cor_mat)/2)  # Convert to dissimilarity
          hc <- hclust(d, method = "ward.D2")
          clusters <- cutree(hc, k = min(n_to_sample, length(pillar_idx)))

          # Sample one variable from each cluster (as much as possible)
          selected_from_pillar <- c()
          for (cluster_id in sample(unique(clusters))) {  # shuffle cluster order
            cluster_vars <- pillar_idx[clusters == cluster_id]
            if (length(cluster_vars) > 0) {
              selected_from_pillar <- c(selected_from_pillar, sample(cluster_vars, 1))
              if (length(selected_from_pillar) >= n_to_sample) break
            }
          }

          # If not enough selected, add randomly
          if (length(selected_from_pillar) < n_to_sample) {
            remaining <- setdiff(pillar_idx, selected_from_pillar)
            additional <- sample(remaining, n_to_sample - length(selected_from_pillar), replace = FALSE)
            selected_from_pillar <- c(selected_from_pillar, additional)
          }

          selected_indices <- c(selected_indices, selected_from_pillar)
        }
      }
    }

    r_m_values[j] <- r_m(selected_indices, data_pc, dm, eig_values)
  }
  
  r_m_results$mean_r_m[i] <- mean(r_m_values)
  r_m_results$max_r_m[i] <- max(r_m_values)
  
  print(paste("  Mean r_m:", round(mean(r_m_values), 4)))
  print(paste("  Max r_m:", round(max(r_m_values), 4)))
}

print("Variable subset size comparison results (r_m):")
print(r_m_results)

# Save results
write.csv(r_m_results, "outputs/r_m_subset_size_comparison.csv", row.names = FALSE)

# Identify a good subset using a simpler approach (rather than full greedy)
print("Identifying a good subset of 15 variables...")

# For computational efficiency, we'll use a correlation-based approach
# Select variables that are most representative of the dataset, ensuring one from each pillar first

# Calculate the average absolute correlation of each variable with all others
cor_matrix <- cor(indicator_data)
mean_abs_cor <- apply(cor_matrix, 2, function(x) mean(abs(x)))

# Start with one indicator from each pillar that has the highest mean correlation
selected_indices <- c()
for (pillar in c("IN.1", "IN.2", "IN.3", "IN.4", "IN.5", "OUT.6", "OUT.7")) {
  pillar_indicator_indices <- pillar_indicators[[pillar]]
  pillar_var_names <- variable_names[pillar_indicator_indices]
  # Get mean correlation for these variables
  pillar_mean_cor <- mean_abs_cor[pillar_var_names]
  # Select the variable with highest mean correlation
  best_var <- names(which.max(pillar_mean_cor))
  selected_indices <- c(selected_indices, which(variable_names == best_var))
}

# Now select the remaining 8 variables (to reach 15) with highest mean correlation 
# among the remaining variables
remaining_vars <- setdiff(variable_names, variable_names[selected_indices])
remaining_mean_cor <- mean_abs_cor[remaining_vars]
top_remaining <- names(sort(remaining_mean_cor, decreasing = TRUE)[1:8])
selected_indices <- c(selected_indices, which(variable_names %in% top_remaining))

top_15_r_m <- r_m(selected_indices, data_pc, dm, eig_values)

print(paste("Top 15 variables by correlation r_m:", round(top_15_r_m, 4)))
print("Selected variables (by index):")
print(selected_indices)

# Get variable names for the selected indices
print("Selected variables (by name):")
selected_variable_names <- variable_names[selected_indices]
print(selected_variable_names)

# Save the subset information
optimal_subset_info <- data.frame(
  index = selected_indices,
  variable_name = selected_variable_names
)
write.csv(optimal_subset_info, "outputs/optimal_15_variable_subset_bypillar.csv", row.names = FALSE)

# Also save just the variable names for easy reference
write.csv(selected_variable_names, "outputs/selected_variables_bypillar.csv",
 row.names = FALSE)

print("Variable perspective analysis completed.")
