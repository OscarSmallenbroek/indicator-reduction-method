###############################################
#### Preliminaries ############################
###############################################

## Load data
library(clustsig)
library(readxl)
library(dplyr)
library(partitionComparison)

Data_P <- read_excel("G:/My Drive/DeskTop on H/Research/JRC Workshop/Project/Analysis/Peter_Test_data.xlsx")
Var_Sel <- read_excel("G:/My Drive/DeskTop on H/Research/JRC Workshop/Project/Analysis/Variable_Selection.xlsx")

## Create a data frame with one row per country and average values
## over all years for each variables

Countries <- unique(Data_P$Country)
Variables <- c("Country",unique(Data_P$Indicatorname))
n_C = length(Countries)
n_V = length(Variables)

data <- matrix(rep("na",(n_C * n_V)),
                   nrow = n_C, ncol = n_V)
data[1:n_C,1] <- Countries
for (i in 1:n_C){
  for (j in 2:n_V){
    sub1 <- Data_P[Data_P$Country  == Countries[i],]
    sub <- sub1[sub1$Indicatorname == Variables[j],]
    avg = mean(sub$Value)
    data[i,j] <- avg
  }
}

data <- as.data.frame(data)
colnames(data) <- Variables
data[,2:ncol(data)] <- lapply(data[,2:ncol(data)], as.numeric)

## Remove EU27
data <- data[data$Country != "EU27",]

## Standardize data
data[, 2:n_V] <- scale(data[, 2:n_V])

## Reduce variables to remove missing values

Broad_vars = pull(Var_Sel[Var_Sel[,2]==1,1],"Variable")
Narrow_vars = pull(Var_Sel[Var_Sel[,3]==1,1],"Variable")

data_broad <- data[,Broad_vars]
data_narrow <- data[,Narrow_vars]


###############################################
#### Compare subset of 45 #####################
###############################################

## Test #1 - correlation between sim matrices
## Need to drop Country variable

# Function for finding the similarity/distance matrix (fast) 

sim_matrix <- function(df) {
  # Note that other methods are possible
  sim_mat <- dist(df, method = "euclidean")
  return(sim_mat)
}

sm_broad <- sim_matrix(data_broad[,-1])
sm_narrow <- sim_matrix(data_narrow[,-1])
s <- cor(c(sm_broad), c(sm_narrow), method = "spearman")
s

## Compare s to a random selection of 45 variables

scores = rep(0,1000)
for (i in 1:1000) {
  col_subset = sample(2:129,44,replace = FALSE)
  ## 2:129 avoids the country column
  reduced <- data_broad[ , col_subset]
  reduced_mat <- sim_matrix(reduced)
  s <- cor(c(sm_broad), c(reduced_mat), method = "spearman")
  scores[i] <- s
}

hist(scores)
narrow_s <- cor(c(sm_broad), c(sm_narrow), method = "spearman")

sum(scores < narrow_s)/1000

## What about significantly fewer

size_test <- function(n_var){
  scores = rep(0,10000)
  for (i in 1:10000) {
    col_subset = sample(2:129,n_var,replace = FALSE)
    ## 2:129 avoids the country column
    reduced <- data_broad[ , col_subset]
    reduced_mat <- sim_matrix(reduced)
    s <- cor(c(sm_broad), c(reduced_mat), method = "spearman")
    scores[i] <- s
  }
  return(max(scores))
}

results = rep(0, 46)
for (i in 5:50){
  results[i-4] = size_test(i)
}
 
plot(5:50,results, xlab = "Number of Variables", ylab = "Max Correlation")

## Compare clusters

broad_results <- simprof(data_broad[,-1], num.expected=1000, num.simulated=999,
                         method.cluster="average", method.distance="euclidean", 
                         method.transform="identity", alpha=0.05)

narrow_results <- simprof(data_narrow[,-1], num.expected=1000, num.simulated=999,
                          method.cluster="average", method.distance="euclidean", 
                          method.transform="identity", alpha=0.05)

## Pull out list of clusters
n_broad <- broad_results$numgroups
broad_clusters <- broad_results$significantclusters
n_narrow <- narrow_results$numgroups
narrow_clusters <- narrow_results$significantclusters

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

rand_ind(broad_results,narrow_results)

###############################################
#### Find best set within Components ##########
###############################################

## Drop first row (Country) and second row (GDP)

Components <- unique(Var_Sel[-(1:2),]$Component)
Dimensions <- unique(Var_Sel[-(1:2),]$Dimension)
Subdimensions <- unique(Var_Sel[-(1:2),]$Subdimension)

## Drop Spillovers
Components <- Components[-6]

## Drop Spillovers and Leisure and Social Interactions
Dimensions <- Dimensions[-c(4,20)]

### Find best set of three for each component

# Function to test a subset of columns for fit with
# the full dataset with regard to the sim matrices
subset_test <- function(col_subset,data) {
  reduced <- data[ , col_subset]
  reduced_mat <- sim_matrix(reduced)
  full_matrix <- sim_matrix(data)
  r <- cor(c(full_matrix), c(reduced_mat), method = "pearson")
  s <- cor(c(full_matrix), c(reduced_mat), method = "spearman")
  #  cat("Pearson =",r,", Spearman =",s,"    ")
  return(c(r,s))
}

# Function to find the best subset of variables within
# a given component of the given size.
find_best_c <- function(size,component) {
  test1 = Var_Sel$Component == component
  test2 = Var_Sel$Broad != 0
  ## Selects all variables under the given component
  variables <- Var_Sel[test1 & test2,]$Variable
  variables <- variables[!is.na(variables)]
  if (length(variables) > size){ # No search needed if there are not enough variables
    data_sub <- data_broad[,variables]
    n_vars <- length(variables)
    combos <- combn(1:n_vars, size, simplify = FALSE)
    best_cor <- 0
    for (i in combos) {
      s <- subset_test(i,data_sub)[2]
      if (s > best_cor) {
        best_cor <- s
        vars = i
      }
    }
  return(c(best_cor,variables[vars]))
  } else {
    return(c(1.0,variables))
  }
}

# Use the above function to find the best 3 variables
# within each component. Save variables and their scores.
# Note that scores only measure relative to the reduced
# matrix for that component
scores = c()
best_vars = c()
for (i in Components[1:6]){
  out <- find_best_c(3, i)
  scores <- c(scores,out[1])
  best_vars <- c(best_vars,out[2:length(out)])
}

# Determine how well this set of variables matches the full set

best_comp_data <- data_broad[ , best_vars]
reduced_mat <- sim_matrix(best_comp_data)
full_matrix <- sim_matrix(data_broad)
s <- cor(c(full_matrix), c(reduced_mat), method = "spearman")
s

# Compare to random set of 18
size_scores <- function(n_var){
  scores = rep(0,1000)
  for (i in 1:1000) {
    col_subset = sample(2:129,n_var,replace = FALSE)
    ## 2:129 avoids the country column
    reduced <- data_broad[ , col_subset]
    reduced_mat <- sim_matrix(reduced)
    s <- cor(c(sm_broad), c(reduced_mat), method = "spearman")
    scores[i] <- s
  }
  return(scores)
}

sum(s < size_scores(18))/1000
size_test(18)

# Determine variation in information for this set

bcd_results <- simprof(best_comp_data, num.expected=1000, num.simulated=999,
                          method.cluster="average", method.distance="euclidean", 
                          method.transform="identity", alpha=0.05)

## Pull out list of clusters
n_broad <- broad_results$numgroups
broad_clusters <- broad_results$significantclusters
n_bcd <- bcd_results$numgroups
bcd_clusters <- bcd_results$significantclusters

rand_ind(broad_results,bcd_results)


###############################################
#### Find best within dimension ###############
###############################################

# Function to find the best subset of variables within
# a given dimension with size ceiling1 variable per dimension
find_best_d <- function(dimension,size) {
  test1 = Var_Sel$Dimension == dimension
  test2 = Var_Sel$Broad != 0
  variables <- Var_Sel[test1 & test2,]$Variable
  variables <- variables[!is.na(variables)]
  n_vars <- length(variables)
  if (n_vars > 1) {
    data_sub <- data_broad[,variables]
    combos <- combn(1:n_vars, size, simplify = FALSE)
    best_cor <- 0
    for (i in combos) {
      s <- subset_test(i,data_sub)[2]
      if (s > best_cor) {
        best_cor <- s
        vars = i
      }
    }
    return(c(best_cor,variables[vars]))
  } else {
    if (n_vars == 1) {return(c(1.0,variables))}
  }
}

# Use the above function to find the best variables
# within each dimension. Save variables and their scores.
# Note that scores only measure relative to the reduced
# matrix for that component
scores = c()
best_vars = c()
for (i in Dimensions){
  out <- find_best_d(i,1)
  if (length(out) > 0){
    scores <- c(scores,out[1])
    best_vars <- c(best_vars,out[2:length(out)])
  }
}

# Determine how well this set of variables matches the full set

best_dim_data <- data_broad[ , best_vars]
reduced_mat <- sim_matrix(best_dim_data)
full_matrix <- sm_broad
s <- cor(c(full_matrix), c(reduced_mat), method = "spearman")
s

# Determine Rand Index for this set

bdd_results <- simprof(best_dim_data, num.expected=1000, num.simulated=999,
                       method.cluster="average", method.distance="euclidean", 
                       method.transform="identity", alpha=0.05)

## Pull out list of clusters
n_bdd <- bdd_results$numgroups
bdd_clusters <- bdd_results$significantclusters

# Compare the clusters
rand_ind(broad_results,bdd_results)

write.csv(best_dim_data, 
          file = "G:/My Drive/DeskTop on H/Research/JRC Workshop/Project/Analysis/Data_Out/bdd_data.csv", 
          row.names = TRUE)

###############################################################
# What is a good Rand Index for matching the broad data set? ##
###############################################################

# Make a clustering with a set of random variables and compare

rand_RI <- function(n_var){
  col_subset = sample(2:129,n_var,replace = FALSE)
  ## 2:129 avoids the country column
  reduced <- data_broad[ , col_subset]
  rand_results <- simprof(reduced, num.expected=1000, num.simulated=999,
                         method.cluster="average", method.distance="euclidean", 
                         method.transform="identity", alpha=0.05)
  RI <- rand_ind(broad_results, rand_results)
  return(RI)
}

rand_RI(19)

ri <- c()
for (i in 1:100){
  out <- rand_RI(19)
  ri <- c(ri,out)
}

hist(ri)
sum(ri>0.88)
###########################################################
## Repeat much of this, but now treat each year and #######
###### country as one sample. #############################
###########################################################

Countries <- unique(Data_P$Country)
Variables <- c("Country","Year",unique(Data_P$Indicatorname))
Years <- unique(Data_P$Year)
n_C = length(Countries)
n_V = length(Variables)
n_Y = length(Years)

data <- matrix(rep("na",(n_C * n_V * n_Y)),
               nrow = (n_C * n_Y), ncol = n_V)

## Populate first two columns with combos of year and country
data[,1] <- rep(Countries, n_Y)
col_2 = c()
for (i in 1:n_Y){
  col_2 = c(col_2,rep((i+2010),n_C))
}
data[,2] <- col_2

for (i in 1:n_C){
  for (j in 1:n_Y){
    for (k in 3:n_V){
      sub1 <- Data_P[Data_P$Country  == Countries[i],]
      sub2 <- sub1[sub1$Year  == Years[j],]
      sub <- sub2[sub2$Indicatorname == Variables[k],]
      loc = (j-1)*n_C + i
      data[loc,k] <- sub$Value[1]
    }
  }
}

data <- as.data.frame(data)
colnames(data) <- Variables
data[,3:ncol(data)] <- lapply(data[,3:ncol(data)], as.numeric)

## Remove EU27
data <- data[data$Country != "EU27",]

## Standardize data
data[, 3:n_V] <- scale(data[, 3:n_V])

## Reduce variables to remove missing values

data_broad <- data[,c("Country","Year",Broad_vars[-1])]
data_narrow <- data[,c("Country","Year",Narrow_vars[-1])]

###############################################
#### Compare subset of 45 #####################
###############################################

## Test #1 - correlation between sim matrices
## Need to drop Country and Year variables

sm_broad <- sim_matrix(data_broad[-(1:2)])
sm_narrow <- sim_matrix(data_narrow[-(1:2)])
s <- cor(c(sm_broad), c(sm_narrow), method = "spearman")
s


## Compare to a random selection of 45 variables

scores = rep(0,1000)
for (i in 1:1000) {
  col_subset = sample(2:129,44,replace = FALSE)
  ## 2:129 avoids the country column
  reduced <- data_broad[ , col_subset]
  reduced_mat <- sim_matrix(reduced)
  s <- cor(c(sm_broad), c(reduced_mat), method = "spearman")
  scores[i] <- s
}

hist(scores)
narrow_s <- cor(c(sm_broad), c(sm_narrow), method = "spearman")

sum(scores < narrow_s)/1000

## What about significantly fewer
results = rep(0, 46)
for (i in 5:50){
  results[i-4] = size_test(i)
}

plot(5:50,results)

###############################################
#### Find best set within Components ##########
###############################################

### Find best set of three for each component
# Note that scores only measure relative to the reduced
# matrix for that component

scores = c()
best_vars = c()
for (i in Components[1:6]){
  out <- find_best_c(3, i)
  scores <- c(scores,out[1])
  best_vars <- c(best_vars,out[2:length(out)])
}

# Determine how well this set of variables matches the full set

best_comp_data <- data_broad[ , best_vars]
reduced_mat <- sim_matrix(best_comp_data)
full_matrix <- sm_broad
s <- cor(c(full_matrix), c(reduced_mat), method = "spearman")
s


# Compare to random set of 19
size_scores <- function(n_var){
  scores = rep(0,1000)
  for (i in 1:1000) {
    col_subset = sample(2:129,n_var,replace = FALSE)
    ## 2:129 avoids the country column
    reduced <- data_broad[ , col_subset]
    reduced_mat <- sim_matrix(reduced)
    s <- cor(c(sm_broad), c(reduced_mat), method = "spearman")
    scores[i] <- s
  }
  return(scores)
}

sum(s < size_scores(19))/1000

###############################################
#### Find best within dimension ###############
###############################################

find_best_d <- function(dimension,size) {
  test1 = Var_Sel$Dimension == dimension
  test2 = Var_Sel$Broad != 0
  variables <- Var_Sel[test1 & test2,]$Variable
  variables <- variables[!is.na(variables)]
  n_vars <- length(variables)
  if (n_vars > 1) {
    data_sub <- data_broad[,variables]
    combos <- combn(1:n_vars, size, simplify = FALSE)
    best_cor <- 0
    for (i in combos) {
      s <- subset_test(i,data_sub)[2]
      if (s > best_cor) {
        best_cor <- s
        vars = i
      }
    }
    return(c(best_cor,variables[vars]))
  } else {
    if (n_vars == 1) {return(c(1.0,variables))}
  }
}

# Use the above function to find the best variable
# within each dimension. Save variables and their scores.
# Note that scores only measure relative to the reduced
# matrix for that subdimension
scores = c()
best_vars = c()
for (i in Dimensions){
  out <- find_best_d(i,1)
  if (length(out) > 0){
    scores <- c(scores,out[1])
    best_vars <- c(best_vars,out[2:length(out)])
  }
}

# Determine how well this set of variables matches the full set

best_dim_data <- data_broad[ , best_vars]
reduced_mat <- sim_matrix(best_dim_data)
full_matrix <- sm_broad
s <- cor(c(full_matrix), c(reduced_mat), method = "spearman")
s

write.csv(best_dim_data, 
          file = "G:/My Drive/DeskTop on H/Research/JRC Workshop/Project/Analysis/Data_Out/bdd_data.CY.csv", 
          row.names = TRUE)

