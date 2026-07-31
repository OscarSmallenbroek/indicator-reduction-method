# Data Utilities for GII Index Reduction Analysis
# Contains data loading and preprocessing functions

#' Load and clean GII data (single entry point)
#' @return List containing standardized data, raw data, and metadata
load_gii_data <- function() {
  # Load metadata
  imeta <- readxl::read_excel("gii-data/imeta.xlsx")
  
  # Load data
  idata <- readxl::read_excel("gii-data/idata.xlsx")
  
  # Identify Level 1 indicators
  level1_indicators <- imeta$iCode[imeta$Level == 1 & imeta$Type == "Indicator"]
  
  # Subset data to Level 1 indicators only
  idata_level1 <- idata[, c("uCode", level1_indicators)]
  
  # Check missing data pattern for Level 1 indicators
  missing_counts <- colSums(is.na(idata_level1))
  total_missing_before <- sum(missing_counts)
  num_vars_missing <- sum(missing_counts > 0)
  num_countries_missing <- sum(rowSums(is.na(idata_level1[, -1])) > 0)
  
  # Perform kNN imputation instead of case-wise deletion
  if (!require("VIM", quietly = TRUE)) {
    install.packages("VIM", repos = "https://cloud.r-project.org/")
    library(VIM)
  } else {
    library(VIM)
  }
  
  # Perform kNN imputation for missing data
  imputed_data <- kNN(idata_level1, 
                      variable = names(idata_level1)[-1],
                      k = 5, imp_var = FALSE)
  
  # Check missing data after imputation
  missing_counts_after <- colSums(is.na(imputed_data))
  total_missing_after <- sum(missing_counts_after)
  num_countries_after <- sum(rowSums(is.na(imputed_data[, -1])) > 0)
  
  # Remove rows that still have missing values after imputation (if any)
  complete_cases <- complete.cases(imputed_data[, -1])
  complete_data <- imputed_data[complete_cases, ]
  
  # Standardize data (z-scores) for indicators only
  indicator_data <- complete_data[, -1]  # Exclude uCode column

  indicator_data <-standardize_with_direction(indicator_data,imeta)
  standardized_data <- as.data.frame(scale(indicator_data))
  standardized_data_with_country <- cbind(Country = complete_data$uCode, standardized_data)
  
  # Return results
  list(
    raw_data = idata,
    imeta = imeta,
    complete_data = complete_data,
    standardized_data = standardized_data_with_country,
    indicator_data_only = standardized_data,
    country_names = complete_data$uCode
  )
}

#' Standardize data with direction reversal for negative indicators
#' @param data Data frame with countries as rows, indicators as columns
#' @param imeta Metadata data frame containing Direction column
#' @return Standardized data frame with direction applied
standardize_with_direction <- function(data, imeta) {
  # Make a copy to avoid modifying original
  data_std <- data
  
  # Identify indicators with negative direction
  neg_indicators <- imeta$iCode[imeta$Direction == -1]
  
  # Apply direction reversal: max(x) - x for negative indicators
  for (ind in neg_indicators) {
    if (ind %in% colnames(data_std)) {
      max_val <- max(data_std[[ind]], na.rm = TRUE)
      data_std[[ind]] <- max_val - data_std[[ind]]
    }
  }
  
  # Standardize (z-score) the data
  indicator_cols <- setdiff(colnames(data_std), "Country")
  data_std[, indicator_cols] <- as.data.frame(scale(data_std[, indicator_cols]))
  
  return(data_std)
}