library(dplyr)
library(stats)

print("Loading prepared GII data for variable perspective analysis...")

################################
# Prepare GII data ------
#################################
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


# standardize the data
indicator_data<-indicator_data |> 
  mutate(across(everything(), ~scale(.)))


###############################################
#### PCA Analysis #############################
###############################################

print("Performing PCA analysis...")

pca_result <- prcomp(indicator_data, center = TRUE, scale. = TRUE)
summary(pca_result)
## Eigenvalues in decreasing order
eig_values <- pca_result$sdev^2

## Multiply the data set by the PC results as matrices to get 
## coordinates for each country in the PCs
## Also take the first 19 as a comparison space
## The first N PC's account for 95% of the variability

pca_summary<-summary(pca_result)
pca_summary$importance 
pca95_N<-which(pca_summary$importance[3, ] > 0.95)[[1]]

dm <- as.matrix(indicator_data)
n_var <- ncol(dm)
pcm <- as.matrix(pca_result$rotation)
data_pc <- dm %*% pcm
data_pc_19 <- data_pc[,1:pca95_N]
eig_val_19 <- eig_values[1:pca95_N]


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
  for (i in 1:length(eig_values)) { # not clear if we should use only the PCs that get 95% of variation of all of the PCs
    r_m.i.squared <- PC_cor(i,subset)
    r_sum <- r_sum + eig_values[i] * r_m.i.squared 
  }
  eig_sum = sum(eig_values)
  return(sqrt(r_sum / eig_sum))
}

####################################################
## SEARCH WITHIN COMPONENTS TO FIND THE BEST########
####################################################

Components <- unique(imeta$NUM[imeta$Level ==3])

## Write each of the key steps in terms of 
## an initial subset of the original data

## Sub versions of the main functions
## Given a component, must subset the data and find the new PCs
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

  ## Selects all variables under the given component
  variables<-imeta |> 
    filter(grepl(pattern = component, NUM) & Level ==1 ) |> 
    pull(iCode)
  nc_var <- length(variables)

  #message('variables for component ' , component, ' are ', paste(variables))
  if (nc_var > k){ # No search needed if there are not enough variables
    dm_sub <- as.matrix(indicator_data[,variables])
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
    return(data.frame(
              best_rm = best_rm, 
              variables = colnames(dm_sub)[sort(vars)]
            ))
  } else {
    return(c(1.0,variables))
  }
}

# Use the above function to find the best 3 variables
# within each component. Save variables and their scores.
# Note that scores only measure relative to the reduced
# matrix for that component

# search for taking 1-5 indicators per component. 

results<-data.frame()

for (ind_n in 1:4) {
message( 'searching for ', ind_n, " indicators per component")
scores = c()
best_vars = c()
for (i in Components){
  out <- best_k_comp(ind_n, i)
  out<-cbind(out,i,ind_n)
  names(out)<-c('best_rm', 'indicators', 'component', 'n_indicators')
  best_vars <- rbind(best_vars,out)
}

best_vars$total_rm<-r_m(best_vars$indicators)
results<-rbind(results, best_vars)

}


readr::write_csv(results, 
          file.path('outputs', 'results_by_dimension.csv'))

# calculate what % of indicators within component was selected
# Create a data frame with pillar information
pillar_df <- imeta |> 
  filter(Level == 1) |>  # only indicators
  mutate(pillar_name = substr(Parent, 2,3),
        NUM = paste0(substr(NUM, 1,2),".", substr(pillar_name,2,2)) ) |> 
  group_by(pillar_name, NUM) |> 
  summarise(n_indicators  = n())

# Calculate percentage of total indicators for each pillar
pillar_df$percent_of_total <- pillar_df$n_indicators / sum(pillar_df$n_indicators) * 100

pillar_df <-pillar_df |> 
  mutate(choose_1 = 1 / n_indicators, 
        choose_2 = 2 / n_indicators,
        choose_3 = 3 / n_indicators,
        choose_4 = 4 / n_indicators, 
  )


## Test against a random subset of 18
scores = c()
for (i in 1:1000) {
  rand = sample(1:ncol(dm),21,replace = FALSE)
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
    dm_sub <- as.matrix(indicator_data[,variables])
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
  