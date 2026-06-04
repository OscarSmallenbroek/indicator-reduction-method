#!/usr/bin/env Rscript
# Script to knit QMD report to Word document

# Install required packages if needed
required_packages <- c("rmarkdown", "knitr")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}

# Knit the QMD file to Word
rmarkdown::render(
  input = "analysis-report.qmd",
  output_format = "word_document",
  output_file = "GII_Index_Reduction_Report.docx",
  clean = TRUE
)

cat("Report successfully knitted to: GII_Index_Reduction_Report.docx\n")