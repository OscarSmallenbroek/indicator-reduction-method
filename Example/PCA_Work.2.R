###############################################
#### Preliminaries ############################
###############################################

## Load data
library(readxl)
library(dplyr)
library(stats)

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
## Average across years for each country
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
#### PCA Analysis #############################
###############################################

pca_result <- prcomp(data_broad[,-1], center = TRUE, scale. = TRUE)
summary(pca_result)
## Eigenvalues in decreasing order
eig_values <- pca_result$sdev^2

## Multiply the data set by the PC results as matrices to get 
## coordinates for each country in the PCs
## Also take the first 19 as a comparison space
## The first 19 PC's account for 95.3% of the variability


dm <- as.matrix(data_broad[,-1])
n_var <- ncol(dm)
pcm <- as.matrix(pca_result$rotation)
data_pc <- dm %*% pcm
data_pc_19 <- data_pc[,1:19]
eig_val_19 <- eig_values[1:19]


## function to find the multiple correlation between a PC and 
## a subset of variables.  This is (r_m)^2_i for the ith PC.
## Inputs are the index of the PC and the names or indices of 
## variables in the subset
PC_cor <- function(PC,subset){
  out <- data_pc[,PC]
  inputs <- dm[,subset]
  matrix_for_analysis <- cbind(out,inputs)
  data <- as.data.frame(matrix_for_analysis)
  model <- lm(out ~ ., data = data)
  model_summary <- summary(model)
  mult_r.squared <- model_summary$r.squared
  return(mult_r.squared)
}

## function to calculate weighted r_m for a proposed subset of variables
r_m <- function(subset){
  r_sum <- 0
  for (i in 1:27) {
    r_m.i.squared <- PC_cor(i,subset)
    r_sum <- r_sum + eig_values[i] * r_m.i.squared 
  }
  eig_sum = sum(eig_values)
  return(sqrt(r_sum / eig_sum))
}

####################################################
## SEARCH WITHIN COMPONENTS TO FIND THE BEST########
####################################################

Components <- unique(Var_Sel[-(1:2),]$Component)
Components <- Components[-6]
## Write each of the key steps in terms of 
## an initial subset of the original data

## Sub versions of the main functions
## Given a component, must subset the data and find the new PCs
#' Compute Multiple Squared Correlation for a Principal Component and Variable Subset
#'
#' This function computes the multiple squared correlation (R-squared) between
#' a specified principal component (PC) and a subset of explanatory variables
#' within a sub-dataset. The R-squared value reflects how much variance in the
#' given PC is explained by the selected variables using a linear model.
#'
#' @param PC Integer index specifying which principal component to analyze.
#' @param subset Character vector or integer indices indicating the subset of variables
#'   to include in the regression model.
#' @param data_pc_sub Matrix or data frame containing principal component scores
#'   (for all components) derived from a sub-dataset.
#' @param dm_sub Matrix or data frame containing the original scaled variables
#'   (excluding response) corresponding to the sub-dataset. Used as predictors.
#'
#' @return Numeric value representing the multiple R-squared from regressing
#'   the specified PC on the given subset of variables.
#'
#' @details This function is intended for internal use during variable selection strategies,
#'   particularly when evaluating how well subsets of variables reconstruct principal
#'   components. It mirrors `PC_cor` but operates on subsetted data structures.
#'
#' @example
#' # Suppose data_pc_sub and dm_sub are already computed for a subset
#' PC_cor_sub(PC = 1, subset = c("Var1", "Var2"), data_pc_sub, dm_sub)
PC_cor_sub <- function(PC,subset,data_pc_sub, dm_sub){
  out <- data_pc_sub[,PC]
  inputs <- dm_sub[,subset]
  matrix_for_analysis <- cbind(out,inputs)
  data <- as.data.frame(matrix_for_analysis)
  model <- lm(out ~ ., data = data)
  model_summary <- summary(model)
  mult_r.squared <- model_summary$r.squared
  return(mult_r.squared)
}

r_m_sub <- function(subset,n_var_sub,data_pc_sub,dm_sub,eig_values_sub){
  r_sum <- 0
  for (i in 1:n_var_sub) {
    r_m.i.squared <- PC_cor_sub(i,subset,data_pc_sub,dm_sub)
    r_sum <- r_sum + eig_values_sub[i] * r_m.i.squared 
  }
  eig_sum = sum(eig_values_sub)
  return(sqrt(r_sum / eig_sum))
}

## Function to find the best k variables for a given component
## Based on correlation with the PCs for that component
best_k_comp <- function(k,component) {
  test1 = Var_Sel$Component == component
  test2 = Var_Sel$Broad != 0
  ## Selects all variables under the given component
  variables <- Var_Sel[test1 & test2,]$Variable
  variables <- variables[!is.na(variables)]
  nc_var <- length(variables)
  if (nc_var > k){ # No search needed if there are not enough variables
    dm_sub <- as.matrix(data_broad[,variables])
    n_var_sub <- min(c(nc_var,27))
    pca_result_sub <- prcomp(dm_sub, center = TRUE, scale. = TRUE)
    ## Eigenvalues in decreasing order
    eig_values_sub <- pca_result_sub$sdev^2
    pcm_sub <- as.matrix(pca_result_sub$rotation)
    data_pc_sub <- dm_sub %*% pcm_sub ## Matrix of PCs
    combos <- combn(1:nc_var, k, simplify = FALSE)
    best_rm <- 0
    test <- 0
    for (i in combos) {
      s <- r_m_sub(i,n_var_sub,data_pc_sub,dm_sub,eig_values_sub)
      test <- test + 1
      if (s > best_rm) {
        best_rm <- s
        vars = i
      }
    }
    return(c(best_rm,colnames(dm_sub)[sort(vars)]))
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
for (i in Components[c(1:6)]){
  out <- best_k_comp(3, i)
  scores <- c(scores,out[1])
  best_vars <- c(best_vars,out[2:length(out)])
}

best_vars
round(as.numeric(scores),3)
r_m(best_vars) ## 0.931
sum(eig_values[1:18])/sum(eig_values)

## Test against a random subset of 18
scores = c()
for (i in 1:1000) {
  rand = sample(1:ncol(dm),18,replace = FALSE)
  scores = c(scores,r_m(rand)) 
}
sum(scores > 0.9)/1000
max(scores)
################################################################
# Function to find the best variable within a given dimension ##
################################################################

Dimensions <- unique(Var_Sel[-(1:2),]$Dimension)
Dimensions <- Dimensions[-c(4,20)]

best_k_dim <- function(k,dimension) {
  test1 = Var_Sel$Dimension == dimension
  test2 = Var_Sel$Broad != 0
  variables <- Var_Sel[test1 & test2,]$Variable
  variables <- variables[!is.na(variables)]
  nd_var <- length(variables)
  if (nd_var > k){ # No search needed if there are not enough variables
    dm_sub <- as.matrix(data_broad[,variables])
    n_var_sub <- min(c(nd_var,27))
    pca_result_sub <- prcomp(dm_sub, center = TRUE, scale. = TRUE)
    ## Eigenvalues in decreasing order
    eig_values_sub <- pca_result_sub$sdev^2
    pcm_sub <- as.matrix(pca_result_sub$rotation)
    data_pc_sub <- dm_sub %*% pcm_sub ## Matrix of PCs
    combos <- combn(1:nd_var, k, simplify = FALSE)
    best_rm <- 0
    test <- 0
    for (i in combos) {
      s <- r_m_sub(i,n_var_sub,data_pc_sub,dm_sub,eig_values_sub)
      test <- test + 1
      if (s > best_rm) {
        best_rm <- s
        vars = i
      }
    }
    return(c(best_rm,colnames(dm_sub)[sort(vars)]))
  } else {
    return(c(1.0,"only one"))
  }
}

# Use the above function to find the best variable
# within each dimension. Save variables and their scores.
# Note that scores only measure relative to the reduced
# matrix for that dimension


scores = c()
best_vars = c()
for (i in Dimensions){
  out <- best_k_dim(1, i)
  scores <- c(scores,out[1])
  best_vars <- c(best_vars,out[2:length(out)])
}

best_vars[-6]
round(as.numeric(scores),3)
r_m(best_vars[-6]) 
sum(eig_values[1:18])/sum(eig_values)

##################################################################
##################################################################
## Repeat, but now treat each year and country as one sample. ####
##################################################################

Countries <- unique(Data_P$Country)
Variables <- c("Country","Year",unique(Data_P$Indicatorname))
Years <- unique(Data_P$Year)
n_C = length(Countries)
n_V = length(Variables)
n_Y = length(Years)

data_cy <- matrix(rep("na",(n_C * n_V * n_Y)),
               nrow = (n_C * n_Y), ncol = n_V)

## Populate first two columns with combos of year and country
data_cy[,1] <- rep(Countries, n_Y)
col_2 = c()
for (i in 1:n_Y){
  col_2 = c(col_2,rep((i+2010),n_C))
}
data_cy[,2] <- col_2

for (i in 1:n_C){
  for (j in 1:n_Y){
    for (k in 3:n_V){
      sub1 <- Data_P[Data_P$Country  == Countries[i],]
      sub2 <- sub1[sub1$Year  == Years[j],]
      sub <- sub2[sub2$Indicatorname == Variables[k],]
      loc = (j-1)*n_C + i
      data_cy[loc,k] <- sub$Value[1]
    }
  }
}

data_cy <- as.data.frame(data_cy)
colnames(data_cy) <- Variables
data_cy[,3:ncol(data_cy)] <- lapply(data_cy[,3:ncol(data_cy)], as.numeric)

## Remove EU27
data_cy <- data_cy[data_cy$Country != "EU27",]

## Standardize data
data_cy[, 3:n_V] <- scale(data_cy[, 3:n_V])

data_broad_cy <- data_cy[,c("Country","Year",Broad_vars[-1])]

###############################################
#### PCA Analysis #############################
###############################################

pca_result <- prcomp(data_broad_cy[,-(1:2)], center = TRUE, scale. = TRUE)
summary(pca_result)
## Eigenvalues in decreasing order
eig_values <- pca_result$sdev^2

## Multiply the data set by the PC results as matrices to get 
## coordinates in the PCs
## First 19 PC's have a cumulative variance of 0.871
## These are the variables to compare with

dm_cy <- as.matrix(data_broad_cy[,-(1:2)])
pcm <- as.matrix(pca_result$rotation)
data_pc_cy <- dm_cy %*% pcm

## function to find the multiple correlation between a PC and 
## a subset of variables.  This is (r_m)^2_i for the ith PC.
PC_cor_cy <- function(PC,subset){
  out <- data_pc_cy[,PC]
  inputs <- dm_cy[,subset]
  matrix_for_analysis <- cbind(out,inputs)
  data <- as.data.frame(matrix_for_analysis)
  model <- lm(out ~ ., data = data)
  model_summary <- summary(model)
  mult_r.squared <- model_summary$r.squared
  return(mult_r.squared)
}

## function to calculate weighted r_m for a proposed subset of variables
r_m_cy <- function(subset){
  r_sum <- 0
  for (i in 1:128) {
    r_m.i.squared <- PC_cor_cy(i,subset)
    r_sum <- r_sum + eig_values[i] * r_m.i.squared 
  }
  eig_sum = sum(eig_values)
  return(sqrt(r_sum / eig_sum))
}

##########################################################
## Test the 44 variable subset ###########################
##########################################################

## Get column indices for 44 variable subset
dm_var <- Var_Sel[Var_Sel$Broad==1,] # Drop missing data
# Select Narrow variables, dropping row 8 (missing data) and country
sub_narrow <- (1:129)[dm_var$Narrow==1][-1] - 1
colnames(dm)[sub_narrow]

r_m_cy(sub_narrow)

## proportion captured by first 44 PC
sum(eig_values[1:44])/sum(eig_values)

# How does a random group of 44 variables work?
r_m_test <- function(k){
  rand_sub <- sample(1:128,k)
  return(r_m_cy(rand_sub))
}

scores <- c()
for (i in 1:100) {
  scores <- c(scores, r_m_test(44))
}

sum(scores > 0.958)

####################################################
## SEARCH WITHIN COMPONENTS TO FIND THE BEST########
####################################################

Components <- unique(Var_Sel[-(1:2),]$Component)
Components <- Components[-6]
## Write each of the key steps in terms of 
## an initial subset of the original data

## Modified versions of the main functions as needed

## Function to find the best k variables for a given component
## Based on correlation with the PCs for that component

best_k_comp_cy <- function(k,component) {
  test1 = Var_Sel$Component == component
  test2 = Var_Sel$Broad != 0
  ## Selects all variables under the given component
  variables <- Var_Sel[test1 & test2,]$Variable
  variables <- variables[!is.na(variables)]
  nc_var <- length(variables)
  if (nc_var > k){ # No search needed if there are not enough variables
    dm_sub <- as.matrix(data_broad_cy[,variables])
    n_var_sub <- nc_var  #No limit on number of variables
    pca_result_sub <- prcomp(dm_sub, center = TRUE, scale. = TRUE)
    ## Eigenvalues in decreasing order
    eig_values_sub <- pca_result_sub$sdev^2
    pcm_sub <- as.matrix(pca_result_sub$rotation)
    data_pc_sub <- dm_sub %*% pcm_sub ## Matrix of PCs
    combos <- combn(1:nc_var, k, simplify = FALSE)
    best_rm <- 0
    test <- 0
    for (i in combos) {
      s <- r_m_sub(i,n_var_sub,data_pc_sub,dm_sub,eig_values_sub)
      test <- test + 1
      if (s > best_rm) {
        best_rm <- s
        vars = i
      }
    }
    return(c(best_rm,colnames(dm_sub)[sort(vars)]))
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
for (i in Components[c(1:6)]){
  out <- best_k_comp_cy(3, i)
  scores <- c(scores,out[1])
  best_vars <- c(best_vars,out[2:length(out)])
}

best_vars
round(as.numeric(scores),3)
r_m_cy(best_vars) ## 0.931
sum(eig_values[1:18])/sum(eig_values)

scores = c()
for (i in 1:1000) {
  rand = sample(1:ncol(dm),18,replace = FALSE)
  scores = c(scores,r_m_cy(rand)) 
}
sum(scores > 0.0.850)/1000
################################################################
# Function to find the best variable within a given dimension ##
################################################################

Dimensions <- unique(Var_Sel[-(1:2),]$Dimension)
Dimensions <- Dimensions[-c(4,20)]

best_k_dim <- function(k,dimension) {
  test1 = Var_Sel$Dimension == dimension
  test2 = Var_Sel$Broad != 0
  variables <- Var_Sel[test1 & test2,]$Variable
  variables <- variables[!is.na(variables)]
  nd_var <- length(variables)
  if (nd_var > k){ # No search needed if there are not enough variables
    dm_sub <- as.matrix(data_broad[,variables])
    n_var_sub <- min(c(nd_var,27))
    pca_result_sub <- prcomp(dm_sub, center = TRUE, scale. = TRUE)
    ## Eigenvalues in decreasing order
    eig_values_sub <- pca_result_sub$sdev^2
    pcm_sub <- as.matrix(pca_result_sub$rotation)
    data_pc_sub <- dm_sub %*% pcm_sub ## Matrix of PCs
    combos <- combn(1:nd_var, k, simplify = FALSE)
    best_rm <- 0
    test <- 0
    for (i in combos) {
      s <- r_m_sub(i,n_var_sub,data_pc_sub,dm_sub,eig_values_sub)
      test <- test + 1
      if (s > best_rm) {
        best_rm <- s
        vars = i
      }
    }
    return(c(best_rm,colnames(dm_sub)[sort(vars)]))
  } else {
    return(c(1.0,"only one"))
  }
}

# Use the above function to find the best variable
# within each dimension. Save variables and their scores.
# Note that scores only measure relative to the reduced
# matrix for that dimension


scores = c()
best_vars = c()
for (i in Dimensions){
  out <- best_k_dim(1, i)
  scores <- c(scores,out[1])
  best_vars <- c(best_vars,out[2:length(out)])
}

best_vars[-6]
round(as.numeric(scores),3)
r_m(best_vars[-6]) 
sum(eig_values[1:18])/sum(eig_values)


###########################################################
## Break out by year ######################################
###########################################################

B18_scores <- c()
p_array <- c()
var_by_year <- array(rep("x",12 * 18), dim = c(18,12))
for (yr in Years){
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
      sub2 <- sub1[sub1$Indicatorname == Variables[j],]
      sub3 <- sub2[sub2$Year == yr,]
      if (length(sub3$Value) > 0) {data[i,j] <- sub3$Value}
      else {data[i,j] <- "NA"}
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
  #### PCA Analysis #############################
  ###############################################
  
  pca_result <- prcomp(data_broad[,-1], center = TRUE, scale. = TRUE)
  summary(pca_result)
  ## Eigenvalues in decreasing order
  eig_values <- pca_result$sdev^2
  
  ## Multiply the data set by the PC results as matrices to get 
  ## coordinates for each country in the PCs

  dm <- as.matrix(data_broad[,-1])
  n_var <- ncol(dm)
  pcm <- as.matrix(pca_result$rotation)
  data_pc <- dm %*% pcm

  
  ####################################################
  ## SEARCH WITHIN COMPONENTS TO FIND THE BEST########
  ####################################################

  scores = c()
  best_vars = c()
  for (i in Components[c(1:6)]){
    out <- best_k_comp(3, i)
    scores <- c(scores,out[1])
    best_vars <- c(best_vars,out[2:length(out)])
  }
  
  var_by_year[,(as.numeric(yr) - 2010)] <- best_vars
  B18_scores <- c(B18_scores,r_m(best_vars)) 
#  sum(eig_values[1:18])/sum(eig_values)
  
  scores <- c()
  for (i in 1:100) {
    rand_sub <- sample(1:128,18)
    scores <- c(scores, r_m(rand_sub))
  }
  
  p <- sum(scores > r_m(best_vars))/100
  p_array <- c(p_array, p)
  
  
}

write.csv(var_by_year, 
          file = "G:/My Drive/DeskTop on H/Research/JRC Workshop/Project/Analysis/Data_Out/Year.PC.Vars.csv", 
          row.names = FALSE)
p_array
round(B18_scores,3)
mean(B18_scores)
