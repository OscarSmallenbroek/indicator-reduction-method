# Package bootstrap for the GII analysis.
#
# The numbered pipeline scripts (01-04) each source config/functions/data_utils
# themselves and attach the libraries they need, so they run standalone under
# Rscript. This script only needs to be run once, to install the packages.

options(repos = c(CRAN = "https://cran.rstudio.com/"))

# Installed on demand; "stats" is a base package and needs no installation.
required_packages <- c("readxl", "dplyr", "VIM", "COINr", "clustsig",
                       "partitionComparison")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing ", pkg, "...")
    install.packages(pkg, dependencies = TRUE)
  }
}

# Optional, only used for ad hoc plotting - not required by the pipeline.
optional_packages <- c("ggplot2", "corrplot")
for (pkg in optional_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Optional package not installed: ", pkg)
  }
}

message("Setup complete. Run scripts/01_data_preparation.R onwards.")
