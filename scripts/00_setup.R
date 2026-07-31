# Setup and Library Loading Script

# Clear environment
rm(list = ls())

# Set CRAN mirror
options(repos = c(CRAN = "https://cran.rstudio.com/"))

# Load required libraries
required_packages <- c("readxl", "dplyr", "stats", "VIM", "COINr")
# Check and install missing packages
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}

# Load additional packages if available, but don't install automatically
additional_packages <- c( "ggplot2", "corrplot")
for (pkg in additional_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    library(pkg, character.only = TRUE)
    print(paste("Loaded package:", pkg))
  } else {
    print(paste("Package", pkg, "not available. Some functionality may be limited."))
  }
}

# Source shared R files
source("scripts/R/config.R")
source("scripts/R/functions.R")
source("scripts/R/data_utils.R")

# Set working directory (if needed)
# setwd("path/to/your/project")

print("Setup complete. Required libraries loaded and shared functions sourced.")
