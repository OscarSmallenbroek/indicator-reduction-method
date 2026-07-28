
# to do list 

- standardize methodology in line with Kenneth, use as templates:
    - Example/PCA_Work2.R
    - Example/Matrix_work.2.R

- differences with Kenneth
    - we do not have time-series, only do one year, so there is no comparison across time. 

- outputs required are two main categories

1. Variable perspective (choosing a subset)
    - choose one variable for each level-2 aggregate (justify with reflective indicator). Choose the variables in the level-2 aggregate set that has the best R_m. 
    - for each subset: test with R_m how much of variation is captured by compete subset, compare it to the number of PCAs needed. 


2. Country perspective
    1. Plot the rank correlation of 1000 random subsets of  k variables for k 10% of total variables to 90% with the similarity matrix of the full dataset. (fig 2 in methods and results v3). DONE
        - save the list of variables that is the best random selection. DONE
    2. Exhaustive search: find best variable subset using dissimilarity matrix rank correlation as metric. DONE
        - Searched for an optimal subset using two strategies. 
            1. included 3 variables per domain (Level 3) or DONE
            2. 1Include 1 variable per subdomain (Level 2). DONE
    3. SIMPROF + agglomerative clustering: shows which countries are the same of different based on a latent state. [check if we are using the same packages and commands as Ken]
        - clustering countries using full dataset and subset. 
            - derive an equiviliance measure within clusters and compare between datasets. Rand Index. 


# functions required
We can use the COINr imeta file which contains meta data on the indicator names and structur of the index. 

- choose X random units from a list of indicators in an imeta file. 
  -  function(x, level, iMeta)
  - function should select x number of indicators from iMeta$iCode that are grouped at Level = level in the iMeta file. 
  - then return the list of indicators. 

- PC_cor
- r_m function 

