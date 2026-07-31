library(dplyr)
library(stats)

print("Loading prepared GII data for variable perspective analysis...")

# Load standardized data
data_standardized <- read.csv("data/gii_data_standardized.csv", stringsAsFactors = FALSE)
print(paste("Loaded data with", nrow(data_standardized), "countries and", ncol(data_standardized)-1, "indicators"))

# Prepare data for analysis (exclude Country column)
country_names <- data_standardized$Country
indicator_data <- data_standardized[, -1]  # Remove Country column
# Load the metadata to identify pillars
imeta <- read.csv("gii-data/imeta.csv", stringsAsFactors = FALSE)
# Extract indicators (rows where Level == 1)
indicators <- imeta[imeta$Level == 1, ]

# Ensure that alll indicators have a positive direction
rev_indicators<-imeta$iCode[imeta$Direction== -1]
indicator_data<-indicator_data |> 
  mutate(across(all_of(rev_indicators), ~ max(.) - .))

# remove all aggregates from the dataset - pillars, sub-index etc.
indicator_data<-indicator_data |> 
  select(all_of(indicators$iCode))
variable_names <- colnames(indicator_data)


print("Performing PCA analysis...")

# Perform PCA on the full dataset -----
pca_result <- prcomp(indicator_data, center = TRUE, scale. = TRUE)
summary_pca <- summary(pca_result)

saveRDS(pca_result, 
      file.path( 'outputs', 'pca_fullGII.rds'))
saveRDS(summary_pca, 
      file.path( 'outputs', 'summary_pca_fullGII.rds'))

# Eigenvalues in decreasing order
eig_values <- pca_result$sdev^2

print("PCA Results:")
print(paste("Number of principal components:", length(eig_values)))
print("Proportion of variance explained by first 10 components:")
print(round(summary_pca$importance[2, 1:min(10, length(eig_values))], 4))


# testing subsets  ---------------


# Function to calculate weighted r_m for a proposed subset of variables
# This implements Equation 1 from the manuscript:
# r_m = sqrt((sum(i=1 to p)(λ_i * (r_m)_i^2)) / (sum(j=1 to p)λ_j))

# First, we need matrices for multiplication
dm <- as.matrix(indicator_data)
pcm <- as.matrix(pca_result$rotation)
data_pc <- dm %*% pcm  # these are scores on all PCs 

# Function to find the multiple correlation between a PC and a subset of variables
PC_cor <- function(PC, subset_indices, data_pc_matrix, dm_matrix) {
  out <- data_pc_matrix[, PC]
  inputs <- dm_matrix[, subset_indices]
  matrix_for_analysis <- cbind(out, inputs)
  data <- as.data.frame(matrix_for_analysis)

  model <- lm(out ~ . , data = data)
  model_summary <- summary(model)
  mult_r.squared <- model_summary$r.squared
  return(mult_r.squared)
}


# Create a data frame with pillar information
pillar_df <- imeta |> 
  filter(Level == 1) |>  # only indicators
  mutate(pillar_name = substr(Parent, 2,3)) |> 
  group_by(pillar_name) |> 
  summarise(n_indicators  = n())

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

# Create an exmple for the manuscript. 
allocate20<-cbind(pillar_df, calculate_allocation(20))
names(allocate20)[4]<-'allocation'
saveRDS(allocate20, file.path('outputs', 'allocation20.rds'))


# Add the allocation to the pillar data frame
# This part will be used in the sampling loop
pillar_df$allocation <- 0  # Initialize

# r_m function ------
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
s
# a function to get the pillar' indicators indexes in the data set
pillar_indicators<-function(pillar){

  imeta |> 
    filter(Level == 1) |> 
    filter(grepl(pillar, Parent)) |> 
    pull(iCode)

}

#test function 
indicator_data[,pillar_indicators("P1")]

# Have to make sure not sample with replacement! 
# IN correlation and cluster based sampling it can select same indicator twice. 
for (i in 1:length(subset_sizes)) {
  size <- subset_sizes[i]
  print(paste("Testing subset size:", size))
  
  # Test 20 random subsets of this size (reduced from 50 for faster execution)
  r_m_values <- numeric(20)
  
  for (j in 1:20) {
    selected_indices <- c()

    # 1. Compute proportional allocation: sample from each pillar by size
    allocation <- calculate_allocation(size)

    # 2. Sample from each pillar according to allocation, ensuring no replacement within pillar
    for (pillar in names(allocation)) {
      n_to_sample <- allocation[pillar]
      pillar_idx <- pillar_indicators(pillar)

      # Handle case where pillar has only one indicator or we need only one
      if (n_to_sample == 1 || length(pillar_idx) == 1) {
        # Just sample one if only needed
        selected_pillar_ind <- sample(pillar_idx, 1, replace = FALSE)
        selected_indices <- c(selected_indices, selected_pillar_ind)
      } else {
        # Use correlation-based clustering within the pillar
        sub_data <- indicator_data[, pillar_idx, drop = FALSE]
        if (nrow(sub_data) < 2 || ncol(sub_data) < 2) {
          # If not enough data, sample randomly without replacement
          selected_indices <- c(selected_indices, sample(pillar_idx, n_to_sample, replace = FALSE))
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
              selected_from_pillar <- c(selected_from_pillar, sample(cluster_vars, 1, replace = FALSE))
              if (length(selected_from_pillar) >= n_to_sample) break
            }
          }

          # If not enough selected, add remaining ones without replacement
          if (length(selected_from_pillar) < n_to_sample) {
            remaining <- setdiff(pillar_idx, selected_from_pillar)
            additional <- sample(remaining, n_to_sample - length(selected_from_pillar), replace = FALSE)
            selected_from_pillar <- c(selected_from_pillar, additional)
          }

          selected_indices <- c(selected_indices, selected_from_pillar)
        }
      }
    }

    # Compute r_m for this subset and store
    r_m_values[j]$r_m <- r_m(selected_indices, data_pc, dm, eig_values)
    r_m_values[j]$indicators<-selected_indices
  }

  r_m_results$mean_r_m[i] <- mean(r_m_values$r_m )
  r_m_results$max_r_m[i] <- max(r_m_values$r_m )    
  r_m_results_full[i]<=r_m_values
  
  print(paste("  Mean r_m:", round(mean(r_m_values), 4)))
  print(paste("  Max r_m:", round(max(r_m_values), 4)))
}
print("Variable subset size comparison results (r_m):")
print(r_m_results)

# Save results
write.csv(r_m_results, "outputs/r_m_subset_size_comparison.csv", row.names = FALSE)
