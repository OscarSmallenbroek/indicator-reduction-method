# GII Index Reduction Analysis - Implementation Plan

## Project Overview
- **Goal**: Replicate Kenneth's PCA/r_m methodology and country perspective analysis on the full GII 2025 dataset (77 indicators, 139 countries)
- **Key difference from Kenneth**: We include ALL 7 pillars (IN.1-IN.5, OUT.6, OUT.7) and ALL 21 sub-pillars (SP1.1-SP7.3) - no spillover exclusions
- **Output**: Results for `GII_Analysis_Report.qmd` + new analyses from todo-list

## Data Structure (from imeta.csv)
| Level | Name | Count | Codes |
|-------|------|-------|-------|
| 1 | Indicators | 77 | Level==1 & Type=="Indicator" |
| 2 | Sub-pillars | 21 | SP1.1–SP7.3 |
| 3 | Pillars | 7 | IN.1–IN.5, OUT.6, OUT.7 |
| 4 | Sub-indices | 2 | Inputs, Outputs |
| 5 | Index | 1 | Global Innovation Index |

## Key Parameters (config.R)
```r
n_indicators = 77
pillars_level3 = c("IN.1","IN.2","IN.3","IN.4","IN.5","OUT.6","OUT.7")
subpillars_level2 = c("SP1.1","SP1.2","SP1.3", ..., "SP7.3")  # 21 total
exhaustive:
  strategy1_k_per_pillar = 3      # 7 × 3 = 21 vars
  strategy2_k_per_subpillar = 1   # 21 × 1 = 21 vars
fig2:
  k_percentiles = seq(0.1, 0.9, 0.1)  # k ∈ {8,15,23,31,39,46,54,62,69}
  n_random_subsets = 1000
```

## Analyses Required (from todo-list)

### 1. Variable Perspective
- **Subset selection**: Choose 1 variable per Level-2 aggregate (sub-pillar) with best R_m
- **Test**: For each subset, compute R_m vs number of PCAs needed
- **Pillar-constrained random subsets**: Proportional allocation + cluster sampling within pillar
- **Exhaustive search**:
  - Strategy 1: Best 3 per pillar (Level 3)
  - Strategy 2: Best 1 per sub-pillar (Level 2)
- **Reflective indicator**: Best single variable per sub-pillar by R_m_sub

### 2. Country Perspective
- **Fig 2**: Rank correlation of 1000 random subsets (k=10%-90% of 77) vs full similarity matrix
- **Exhaustive search** (same two strategies as above) using dissimilarity matrix rank correlation
- **SIMPROF + agglomerative clustering**: Cluster countries on full data vs subset, compare with Rand Index

## File Structure After Refactoring

```
scripts/
├── R/
│   ├── config.R              # All parameters, paths, constants
│   ├── functions.R           # Consolidated shared functions
│   └── data_utils.R          # load_gii_data(), standardize_with_direction()
├── 00_setup.R                # Libraries, source config/functions/data_utils
├── 01_data_preparation.R     # Thin wrapper calling load_gii_data()
├── 02_eda.R                  # PCA components for QMD
├── 03_country_perspective.R  # Fig 2, Exhaustive (2 strategies), SIMPROF+Rand
├── 04_variable_perspective.R # Consolidated: pillar-constrained, exhaustive, reflective
├── 05_indicator_selection_gii.R  # Uses select_indicators_by_pca.R
├── 05_results.R              # Compile outputs, generate plots
├── master_script.R           # Updated sources
├── select_indicators_by_pca.R  # Standalone PCA selection (keep)
└── test_indicator_selection.R  # Test script (keep)

# DELETED (merged into 04_variable_perspective.R):
# 04_variable_perspective_bypillar.R
# 04_variable_perspective_exhaustive_search.R  
# 04_variable_perspective_random.R
```

## Output Files Required by QMD

| File | Produced By | Status |
|------|-------------|--------|
| `outputs/pca_components.csv` | 02_eda.R | 🆕 |
| `outputs/subset_size_comparison.csv` | 03_country_perspective.R | ✅ |
| `outputs/r_m_subset_size_comparison.csv` | 04_variable_perspective.R | 🆕 |
| `outputs/selected_variables.csv` | 04_variable_perspective.R | 🆕 |
| `outputs/optimal_15_variable_subset.csv` | 04_variable_perspective.R | 🆕 |
| `outputs/allocation20.rds` | 04_variable_perspective.R (save 21-var allocation) | 🆕 |
| `outputs/results_by_dimension.csv` | 04_variable_perspective.R | 🆕 |

## New Output Files (from todo-list)

| File | Produced By | Purpose |
|------|-------------|---------|
| `outputs/fig2_rank_correlations.csv` | 03_country_perspective.R | Fig 2 data |
| `outputs/fig2_best_subsets/k_XX.csv` | 03_country_perspective.R | Best subsets per k |
| `outputs/exhaustive_strategy1_pillar.csv` | 03_country_perspective.R | 3 per pillar |
| `outputs/exhaustive_strategy2_subpillar.csv` | 03_country_perspective.R | 1 per sub-pillar |
| `outputs/exhaustive_combined_subsets.csv` | 03_country_perspective.R | Combined 21-var subsets |
| `outputs/clustering_rand_index.csv` | 03_country_perspective.R | Rand Index (all strategies) |
| `outputs/clustering_full.rds` | 03_country_perspective.R | Full simprof |
| `outputs/clustering_subset_<Strategy>.rds` | 03_country_perspective.R | Per-strategy subset simprof |
| `outputs/cluster_equivalence.csv` | 03_country_perspective.R | Within-cluster equiv |
| `outputs/variable_exhaustive_strategy1.csv` | 04_variable_perspective.R | Var perspective strat 1 |
| `outputs/variable_exhaustive_strategy2.csv` | 04_variable_perspective.R | Var perspective strat 2 |
| `outputs/best_per_subpillar.csv` | 04_variable_perspective.R | Reflective indicators |
| `outputs/r_m_vs_pca_components.csv` | 04_variable_perspective.R | R_m vs # PCs |

## Consolidated Functions (functions.R)

| Function | Purpose |
|----------|---------|
| `sim_matrix(df, method)` | Distance matrix |
| `spearman_dist_cor(full_dist, sub_dist)` | Spearman corr of distance vectors |
| `PC_cor(pc_idx, subset_idx, data_pc, dm)` | R² of PC vs subset |
| `r_m(subset_idx, data_pc, dm, eig_values)` | Weighted R_m metric |
| `PC_cor_sub / r_m_sub` | Sub-PCA versions |
| `get_pca_components(data, threshold)` | PCA + n components |
| `pillar_indicators(pillar_code, imeta)` | Level-1 indicators under pillar |
| `subpillar_indicators(subpillar_code, imeta)` | Level-1 indicators under sub-pillar |
| `best_k_within_group(k, group_codes, imeta, data, level)` | Exhaustive within group |
| `run_simprof(data, ...)` | clustsig::simprof wrapper |
| `rand_index_from_simprof(res1, res2)` | Rand Index comparison |
| `cluster_equivalence(res1, res2)` | Within-cluster equivalence |
| `proportional_allocation(target, imeta, groups, min_per)` | Sampling allocation |
| `load_gii_data()` | Single entry point for data |
| `standardize_with_direction(data, imeta)` | Reverse Direction==-1 + scale |

## Computational Notes
- **True exhaustive across all groups is impossible** (combinatorial explosion)
- **Kenneth's approach**: Exhaustive WITHIN each group (pillar or sub-pillar), then combine best from each group
- **Fig 2**: 9 k-values × 1000 subsets = 9,000 distance matrices (~2-5 min)
- **SIMPROF**: 2 runs × 1000 simulations (~1-2 min)
- **Exhaustive within group**: Max C(15,3)=455 combos per pillar, C(5,1)=5 per sub-pillar - very fast

## Direction Reversal
- Indicators with `Direction == -1` in imeta: `max(.) - .` transformation
- Applied in `standardize_with_direction()` before scaling
- Used by ALL analysis scripts via shared `data_utils.R`

## Next Implementation Steps
1. Create `scripts/R/config.R` ✅
2. Create `scripts/R/functions.R` (consolidate all duplicated functions)
3. Create `scripts/R/data_utils.R` (load_gii_data, standardize_with_direction)
4. Update `00_setup.R` to source new files
5. Refactor `01_data_preparation.R` → thin wrapper
6. Refactor `02_eda.R` → thin wrapper + pca_components.csv
7. Refactor `03_country_perspective.R` → Fig 2, Exhaustive, SIMPROF
8. Consolidate `04_variable_perspective.R` → merge 3 scripts
9. Update `05_indicator_selection_gii.R` and `05_results.R`
10. Update `master_script.R`
11. Delete 3 old 04_* scripts
12. Run master script, verify all outputs exist for QMD