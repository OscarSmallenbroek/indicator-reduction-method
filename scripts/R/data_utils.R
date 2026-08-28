# Data Utilities for GII Index Reduction Analysis
# Contains data loading and preprocessing functions

#' Load and clean GII data (single entry point)
#' @return List with the standardized data, raw data, metadata, and the
#'   missing-data counts describing the imputation that was performed
load_gii_data <- function() {
  # Load metadata
  imeta <- readxl::read_excel(CONFIG$paths$imeta_raw)
  # strip stray annotation marks (e.g. "...business†", "...tariff rate*")
  imeta$iName <- gsub("[*†]", "", imeta$iName)
  imeta$iName <- trimws(imeta$iName)
  imeta$iName <- trimws(gsub("\\s+", " ", imeta$iName)) # collapse repeated spaces, trim ends
  # Load data
  idata <- readxl::read_excel(CONFIG$paths$data_raw)

  # Identify Level 1 indicators
  level1_indicators <- imeta$iCode[imeta$Level == 1 & imeta$Type == "Indicator"]

  # Subset data to Level 1 indicators only
  idata_level1 <- idata[, c("uCode", level1_indicators)]

  # Missing data pattern before imputation. These counts are the source of
  # data/imputation_report.csv, which the report's appendix quotes, so they
  # are returned rather than only printed.
  missing_counts <- colSums(is.na(idata_level1[, -1]))
  missingness <- data.frame(
    Metric = c("Total Missing Values", "Variables with Missing Data",
               "Countries with Missing Data"),
    Count = c(sum(missing_counts), sum(missing_counts > 0),
              sum(rowSums(is.na(idata_level1[, -1])) > 0)),
    stringsAsFactors = FALSE
  )

  # kNN imputation rather than case-wise deletion
  imputed_data <- VIM::kNN(idata_level1,
                           variable = names(idata_level1)[-1],
                           k = 5, imp_var = FALSE)

  # Remove rows that still have missing values after imputation (if any)
  complete_cases <- complete.cases(imputed_data[, -1])
  complete_data <- imputed_data[complete_cases, ]

  # Standardize data (z-scores) for indicators only
  indicator_data <- complete_data[, -1]  # Exclude uCode column

  indicator_data <- standardize_with_direction(indicator_data, imeta)
  # standardize_with_direction() already z-scores; this second pass is a no-op
  # mathematically but is kept because it is not bit-identical, and exhaustive
  # search has exact ties that a 1-ulp change would silently flip.
  standardized_data <- as.data.frame(scale(indicator_data))
  standardized_data_with_country <- cbind(Country = complete_data$uCode, standardized_data)

  # Return results
  list(
    raw_data = idata,
    imeta = imeta,
    complete_data = complete_data,
    standardized_data = standardized_data_with_country,
    indicator_data_only = standardized_data,
    country_names = complete_data$uCode,
    missingness = missingness
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