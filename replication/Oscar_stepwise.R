###############################################
#### Preliminaries ############################
###############################################

## Load data
library(dplyr)

path<-file.path(here::here(), "replication")
data_broad <- read.csv(file.path(path, "gii_data_standardized.csv"))
pillar_data <- read.csv(file.path(path, "Pillars.csv"))

Countries <- data_broad$Country
Variables <- colnames(data_broad)   # Country + the 78 indicators
## Pillars.csv has one row per indicator and no Country row, so the pillar
## flags line up with Variables[-1], not with Variables. Index the indicators
## by name off pillar_data instead of positionally off Variables.
Indicators <- pillar_data$iCode
stopifnot(identical(Indicators, setdiff(Variables, "Country")))

## Numeric-only copy of the data. dist() coerces the character Country column
## to NA, and the as.matrix() round-trip that follows costs ~5e-9 of precision,
## which is enough to flip exact ties between candidate subsets.
indicator_data <- data_broad[, Indicators]

n_C = length(Countries)
n_V = length(Indicators)

## Pillars.csv numbers Subpillar 1-3 *within* each pillar, so Subpillar alone
## groups across pillars (27/26/25 indicators). Only Pillar is globally unique;
## build an explicit id for sub-pillar level work.
pillar_data$SubpillarID <- paste(pillar_data$Pillar, pillar_data$Subpillar, sep = ".")


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

## Spearman is Pearson on ranks, and the ranks of a reference distance vector do
## not change as candidate subsets are swapped in and out. Rank each reference
## once, here, so that the search loops below never re-rank it.
rank_dist <- function(df) {
  return(rank(c(sim_matrix(df))))
}


sm_broad <- sim_matrix(indicator_data)
sm_broad_ranks <- rank(c(sm_broad))

# Test a random smaller subset and report correlation
# Sample from the indicator names so that every indicator - including the last
# column, AppCrea - can be drawn.

sub_test <- function(n) {
  cols <- sample(Indicators, n)
  data_narrow <- data_broad[,cols]
  s <- cor(sm_broad_ranks, rank_dist(data_narrow))
  return(s)
}

## Distribution of s for a random selection of 42 variables

scores = rep(0,1000)
for (i in 1:1000) {
  s <- sub_test(21)
  scores[i] <- s
}

hist(scores)

#############################################################
#### Find best set with three variables per pillar ##########
#############################################################

pillars <- 1:7

### Find best set of three for each component

# Function to test a subset of columns for fit with
# the full dataset with regard to the sim matrices
# full_ranks is the ranked distance vector of the reference dataset. It is
# constant for a given `data`, so a caller scoring many subsets against the same
# reference passes it in once rather than paying for it on every call. The
# default reproduces the standalone behaviour.
subset_test <- function(col_subset, data, full_ranks = rank_dist(data)) {
  # Guard against a character column (e.g. Country) reaching dist().
  stopifnot(all(vapply(data, is.numeric, logical(1))))
  reduced <- data[ , col_subset]
  s <- cor(full_ranks, rank_dist(reduced))
  return(s)
}

# Function to find the best subset of variables within
# a given pillar (given by number) of the given size.
# Returns a list so that the score stays numeric.
find_best_p <- function(size,pil_num) {
  ## Selects all variables under the given component
  variables <- Indicators[pillar_data$Pillar==pil_num]
  data_sub <- indicator_data[,variables]
  full_ranks <- rank_dist(data_sub)   ## constant across every combination
  n_vars <- length(variables)
  combos <- combn(seq_len(n_vars), size, simplify = FALSE)
  best_cor <- -Inf
  vars <- integer(0)
  for (i in combos) {
    s <- subset_test(i, data_sub, full_ranks)
    if (s > best_cor) {
      best_cor <- s
      vars = i
    }
  }
  return(list(score = best_cor, vars = variables[vars]))
}

find_best_p(3,1)

# Use the above function to find the best 3 variables
# within each component. Save variables and their scores.
# Note that scores only measure relative to the reduced
# matrix for that component
scores = numeric(0)
best_vars = c()
for (i in 1:7){
  out <- find_best_p(3, i)
  scores <- c(scores,out$score)
  best_vars <- c(best_vars,out$vars)
}
print(best_vars)

# Determine how well this set of variables matches the full set
s_best <- subset_test(best_vars, indicator_data, sm_broad_ranks)
s_best

# save the output
imeta <- readxl::read_excel("gii-data/imeta.xlsx")
  # strip stray annotation marks (e.g. "...business†", "...tariff rate*")
  imeta$iName <- gsub("[*†]", "", imeta$iName)
  imeta$iName <- trimws(imeta$iName)
  imeta$iName <- trimws(gsub("\\s+", " ", imeta$iName)) # collapse repeated spaces, trim ends

flat_three_pillar<-
  imeta |> 
  filter(iCode %in% best_vars) |> 
  select(iCode, iName, Parent) |> 
  mutate(Spearman_cor = s_best)

write.csv(flat_three_pillar , 
file = file.path(here::here(), 'replication', 'best3pillar.csv'))

# Compare to random set of 21
scores = rep(0,1000)
for (i in 1:1000) {
  s <- sub_test(21)
  scores[i] <- s
}

hist(scores)
sum(scores>s_best)/1000

#############################################################################
######## Step-wise optimization for best three by pillar. ###################
######## Comparison will always be to the whole dataset. ####################
#############################################################################

best_n_step <- function(n, n_rounds) {
  ## start with n random variables from each pillar
  curt_vars <- c()
  for (i in 1:7){
    pot_vars <- Indicators[pillar_data$Pillar==i]
    curt_vars <- c(curt_vars,sample(pot_vars,n))
  }
  curt_s <- subset_test(curt_vars, indicator_data, sm_broad_ranks)
  print(curt_s)
  ## Swap one at a time and test for improvement
  to_beat <- curt_s
  for (r in 1:n_rounds){
    print(r)
    stop_test <- 0
    for (i in 1:7){
      pot_vars <- Indicators[pillar_data$Pillar==i]
      n_poss <- length(pot_vars)
      for (j in 1:n){
        for (k in 1:n_poss) {
          c_v <- curt_vars[(n*(i-1) + j)] ## current variable under consideration
          t_v <- pot_vars[k] ## variable to consider
          ## Skip a no-op swap, and refuse one that would put the same
          ## indicator in two slots: a duplicated column is double-weighted in
          ## the Euclidean distance, so the climber will otherwise buy score
          ## with a set that has 21 slots but fewer than 21 indicators.
          if (c_v != t_v && !(t_v %in% curt_vars))  {
            test_vars <- curt_vars
            test_vars[(n*(i-1) + j)] <- t_v
            s_test <- subset_test(test_vars, indicator_data, sm_broad_ranks)
            if (s_test > to_beat) {
              curt_vars <- test_vars
              to_beat <- s_test
              stop_test <- 1
            }
          } 
        }
      }
    }
    if (stop_test == 0) {break}
  }
  print(curt_vars)
  print(to_beat)
  return(list(curt_vars, to_beat))
}

results= list()
for (i in 1:5){
  best_n_step(3,100)
}


