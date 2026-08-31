## Preliminaries 

library(dplyr)
library(StatMatch)
library(MASS)

######################################################################################
########### Necessary functions. #####################################################
######################################################################################

sim_matrix <- function(df) {
  # Note that other methods are possible
  sim_mat <- dist(df, method = "euclidean")
  return(sim_mat)
}

## For all methods, variables are assumed broken into groups
## where the size of each group is given by the array n_var

### Method 1: Average values within dimension for each sample and rescale. 
### Find distance matrix using rescaled, average values. 

DM_1 <- function(data, n_var) {
  n_samp <- nrow(data)
  comp_avg_data <- matrix(nrow = n_samp, ncol = length(n_var))
  dividers <- cumsum(n_var)
  breaks <- c(0,dividers)
  ## For each component, average over all variables
  for (i in 1:length(n_var)) {
    data_sub <- data[,(breaks[i]+1):breaks[i+1]]
    comp_avg_data[,i] <- rowMeans(data_sub)
  }

  ## rescale
  comp_avg_data <- scale(comp_avg_data)

  ## make sim matrix and find correlation
  DM_comp_avg <- sim_matrix(as.data.frame(comp_avg_data))
  DM_data <- sim_matrix(data)  
  s <- cor(c(DM_data), c(DM_comp_avg), method = "spearman")
  return(list(
    corr = s,
    DM = DM_comp_avg
    ))
}

### Method 2: Calculate distance matrix within each component.
### Take the squared distance matrix.  Divide by number of variables with the component.
### Sum the results across all components and take the square root. 

DM_2 <- function(data, n_var) {
  ## For each component, find the distance matrix using only those variables
  ## Sum the squared distances divided by the number of variables

  n_samp <- nrow(data)
  total_matrix <- matrix(0,nrow = n_samp, ncol = n_samp)
  dividers <- cumsum(n_var)
  breaks <- c(0,dividers)

  for (i in 1:length(n_var)) {
    data_sub <- data[,(breaks[i]+1):breaks[i+1]]
    dist_mat <- as.matrix(sim_matrix(data_sub))
    total_matrix <- total_matrix + dist_mat^2/n_var[i]
  }

  ## Find distance matrix as the root of the average squared distances
  comp_DM <- as.dist(sqrt(total_matrix))

  ## Find correlation
  DM_data <- sim_matrix(data)
  s <- cor(c(DM_data), c(comp_DM), method = "spearman") 
  return(list(
    corr = s,
    DM = comp_DM
    ))
}

### Method 3: Mahalanobis Distance (MD) 

DM_3 <- function(data) {
  MD_DM <- as.dist(mahalanobis.dist(data))
  DM_data <- sim_matrix(data)
  s <- cor(c(DM_data), c(MD_DM), method = "spearman") 
  return(list(
    corr = s, 
    DM = MD_DM
    ))
}

######################################################################################
####### Test 1. Multivariate normally distributed data.###############################
######################################################################################

## Create simulated data

## Dataset properties

# array showing the number of variables in each category
n_var = c(8,10,5,3,5,5,11,2,15,7,7,5,10,11,7,4,5,2,7)
w_in_corr = 0.7 ## Correlation coefficient within a category
btw_corr = 0.3 ## Correlation coefficient between categories
n_samp = 200 # number of samples

run_model_1() # Run function below first

run_model_1 <- function() {
  ## Set mean for all variables to be 0
  mu <- rep(0,sum(n_var))
  
  ## Set covariance matrix.  All variables to have variance = 1
  Sigma <- matrix(rep(0,sum(n_var)^2), nrow = sum(n_var), ncol = sum(n_var))
  dividers <- cumsum(n_var)
  for (i in 1:sum(n_var)){
    for (j in 1:sum(n_var)){
      if ((which.max(dividers >= i)) == (which.max(dividers >= j))) {
        if (i == j) { Sigma[i,j] <- 1}
        else {Sigma[i,j] <- w_in_corr}
      }
      else {Sigma[i,j] <- btw_corr}
    }
  }
  
  ## Random draws
  data <- as.data.frame(mvrnorm(n_samp, mu, Sigma))
  
  ## standardize data
  data <- scale(data)
  
  print("Correlations with the original matrix for all three methods:")
  print(c(
    DM_1(data,n_var)$corr,
    DM_2(data,n_var)$corr,
    DM_3(data)$corr
  ))
}  

######################################################################################
####### Test 2. Simulated set of latent variables each with a set of indicators ######
######################################################################################

## Create simulated data

n_lat <- 5  ## Number of latent variables for each sample
n_var <- c(10,10,5,5,2) ## Arraying showing number of indicators for each variable
ind_cor <- 0.8  ## Desired correlation between indicators and the latent variable
lat_cor <- 0.5  ## correlation between latent variables
n_samp <- 200 ## Number of samples in data set

run_model_2() # Run function below first

run_model_2 <- function() {
  ## create truth  
  ## Set mean for all variables to be 0
  mu <- rep(0,n_lat)
  
  ## Set covariance matrix.  All variables to have variance = 1
  Sigma <- matrix(rep(0,n_lat^2), nrow = n_lat, ncol = n_lat)
  for (i in 1:n_lat){
    for (j in 1:n_lat){
      if (i == j) { Sigma[i,j] <- 1}
      else {Sigma[i,j] <- lat_cor}
    }
  }
  
  ## Random draws
  data_lat <- as.data.frame(mvrnorm(n_samp, mu, Sigma))
  
  ## scale and make distance matrix
  data_lat <- scale(data_lat)
  DM_lat <- sim_matrix(data_lat)
  
  
  ## create the observed indicators
  data_ind <- matrix(rep(0, n_samp * sum(n_var)), nrow = n_samp, ncol = sum(n_var))
  
  dividers <- cumsum(n_var)
  breaks <- c(0,dividers)
  s_ind <- sqrt(1/ind_cor^2 - 1)
  ## For each component, populate the necessary columns
  ## with draws that are correlated with the latent variable
  for (i in 1:length(n_var)) {
    for (j in (breaks[i]+1):breaks[i+1]) {
      data_ind[,j] <- data_lat[,i] + rnorm(n_samp,0,s_ind)
    }
  }
  
  ## scale and make distance matrix
  data_ind <- scale(data_ind)
  DM_ind <- sim_matrix(data_ind)
  
  print("Matches between adjusted matrices and indicator matrix:")
  
  print(c(
    DM_1(data_ind,n_var)$corr,
    DM_2(data_ind,n_var)$corr,
    DM_3(data_ind)$corr
  ))
  
  print("Matches with latent variable matrix:")
  
  print(c(
    cor(c(DM_1(data_ind,n_var)$DM),c(DM_lat)),
    cor(c(DM_2(data_ind,n_var)$DM),c(DM_lat)),
    cor(c(DM_3(data_ind)$DM),c(DM_lat))
  ))
  cat(
    "Match between indicator matrix and the latent variable matrix:",
    cor(c(DM_lat),c(DM_ind))
  )
}
######################################################################################
##### Test 3. Simulated set of latent variables assumed correlated with development ##
##### Development is a latent variable with latent correlates. #######################
##### Each correlate has a given set of correlated indicators. #######################
######################################################################################

## Create simulated data

n_lat <- 5  ## Number of latent variables for each sample
n_var <- c(30,5,5,5,5) ## Arraying showing number of indicators for each variable
ind_cor <- 0.8  ## Desired correlation between indicators and the latent variable
lat_cor <- 0.8  ## Correlation with development variable
n_samp <- 200 ## Number of samples in data set

run_model_3() # Need to run the function below first

run_model_3 <- function() {
  ####################################
  ## Create Truth
  ####################################
  
  ## Determine development and create distance matrix
  
  development <- rnorm(n_samp,0,1)
  DM_develop <- sim_matrix(development)
  
  ## Create latent component variables
  
  data_lat <- matrix(rep(0, n_samp * n_lat), nrow = n_samp, ncol = n_lat)
  s_lat <- sqrt(1/lat_cor^2 - 1)
  for (i in 1:n_lat){
    data_lat[,i] <- development + rnorm(n_samp,0,s_lat)
  }
  
  ## scale and make distance matrix
  data_lat <- scale(data_lat)
  DM_lat <- sim_matrix(data_lat)
  
  ## create the observed indicators
  data_ind <- matrix(rep(0, n_samp * sum(n_var)), nrow = n_samp, ncol = sum(n_var))
  
  dividers <- cumsum(n_var)
  breaks <- c(0,dividers)
  s_ind <- sqrt(1/ind_cor^2 - 1)
  ## For each component, populate the necessary columns
  ## with draws that are correlated with the latent variable
  for (i in 1:length(n_var)) {
    for (j in (breaks[i]+1):breaks[i+1]) {
      data_ind[,j] <- data_lat[,i] + rnorm(n_samp,0,s_ind)
    }
  }
  
  ## scale and make distance matrix
  data_ind <- scale(data_ind)
  DM_ind <- sim_matrix(data_ind)
  
  print("Matches between adjusted matrices and indicator matrix:")
  print(c(
    DM_1(data_ind,n_var)$corr,
    DM_2(data_ind,n_var)$corr,
    DM_3(data_ind)$corr
  ))
  
  print("Matches with latent variable matrix:")
  
  print(c(
    cor(c(DM_1(data_ind,n_var)$DM),c(DM_lat)),
    cor(c(DM_2(data_ind,n_var)$DM),c(DM_lat)),
    cor(c(DM_3(data_ind)$DM),c(DM_lat))
  ))
  
  print("Matches with development matrix:")
  
  print(c( 
    cor(c(DM_1(data_ind,n_var)$DM),c(DM_develop)),
    cor(c(DM_2(data_ind,n_var)$DM),c(DM_develop)),
    cor(c(DM_3(data_ind)$DM),c(DM_develop))
  ))
  
  cat(
    "Match between latent variable matrix and development matrix:",
    cor(c(DM_lat),c(DM_develop),)
  )
}
