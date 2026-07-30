# GII Index Reduction Analysis Project

## Project Context

This project aims to develop and apply a method for assessing to what extent the number of indicators in an index can be reduced without significantly compromising the discriminatory power of the index. 

The methodology has been previously tested on the Sustainable and Inclusive Well-being Index by the European Commission/JRC. This project will replicate that workflow on the Global Innovation Index (GII) data.
##
## Methodology Overview

The approach involves two main perspectives:

1. **Country Perspective**: Assessing how similar countries are using all indicators versus a subset, using distance matrices and clustering techniques
2. **Variable Perspective**: Using Principal Component Analysis (PCA) to measure how much information a subset captures compared to the full set

### Country perspective: flat vs. proportional indicator allocation

Within the country perspective, indicator subsets are chosen via exhaustive within-group search, using two allocation strategies:

- **Strategies 1 & 2** (`scripts/03_country_perspective.R`) allocate a flat number of indicators to every pillar (3 each) or sub-pillar (1 each), ignoring how unbalanced the GII's pillar/sub-pillar sizes actually are (pillars range from 6 to 15 indicators; sub-pillars from 2 to 5).
- **Strategies 1a & 2a** (`scripts/03b_country_perspective_proportional.R`) instead allocate a 42-indicator budget across pillars/sub-pillars proportionally to their size, using a largest-remainder (Hamilton) apportionment — groups with a very small share of the index can receive zero indicators, while larger groups get correspondingly more.

Results from both are documented in the GII Analysis report.

## Objectives

1. Apply the established methodology to GII data for the year 2025
2. Identify a reduced set of indicators that maintain the discriminatory power of the full index
3. Evaluate different subset sizes and compositions to find optimal balance between simplicity and information retention
4. Document the findings and methodology for future applications

## Workflow

The project will follow these steps:
1. Data exploration and preparation
2. Time-averaged analysis (single year = 2025)
3. Country perspective analysis (distance matrices and clustering)
4. Variable perspective analysis (PCA-based methods)
5. Generation of optimal indicator subsets
6. Documentation of results and recommendations

Results are documented in the GII Analysis report in the main directory. 
```
