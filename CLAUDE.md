# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

This is an R-based statistical analysis project (not software with a build/test pipeline). It develops and applies a method for determining how far the number of indicators in a composite index can be reduced without significantly compromising the index's discriminatory power. The method was previously applied to the EU JRC's Sustainable and Inclusive Well-being Index; this project replicates it on the **Global Innovation Index (GII)** 2025 data (77 Level-1 indicators, ~139 countries).

There is no package manager / lockfile (no `renv`, no `.Rproj`). Scripts are run directly with `Rscript` or sourced interactively; missing packages are auto-installed by `scripts/00_setup.R`.

## Running the analysis

Scripts must be run with the repo root as the working directory (they use relative paths like `"scripts/R/config.R"` and `"gii-data/idata.xlsx"`).

```bash
Rscript scripts/00_setup.R
```

Then run the numbered pipeline scripts in order — each sources `00_setup.R`'s config/functions/data_utils indirectly and expects to be run from repo root:

```bash
Rscript scripts/01_data_preparation.R
Rscript scripts/02_eda.R
Rscript scripts/03_country_perspective.R
```

There is no `04_variable_perspective.R` yet (planned per [IMPLEMENTATION_PLAN.md](IMPLEMENTATION_PLAN.md) but not committed) — the variable-perspective analysis currently lives only in the legacy scripts under `scripts/legacy/`.

To render the report (requires Quarto + a LaTeX/PDF engine, since output format is `pdf`):

```bash
quarto render GII_Analysis_Report.qmd
```

There is no automated test suite; `Example/clustsig-optimization/test_clustsig.R`, `test_runner.R`, and `test_windows_fix.R` are ad hoc test scripts for the vendored `clustsig` package fork only, not for this project's own analysis code.

## Architecture

### Two analytical perspectives

The methodology (see [src/method-explanation.md](src/method-explanation.md) and [todo-list.md](todo-list.md)) evaluates candidate indicator subsets from two angles:

1. **Country perspective** (`scripts/03_country_perspective.R`) — do countries end up "the same distance apart" using the subset as with the full indicator set? Implemented via:
   - Spearman rank correlation between full vs. subset Euclidean distance matrices, evaluated over 1000 random subsets at k = 10%–90% of indicators (produces `outputs/fig2_rank_correlations.csv` and `outputs/fig2_best_subsets/k_*.csv`).
   - Exhaustive within-group search for the best-scoring subset (see below) — Strategy 1: best 3 indicators per pillar (Level 3); Strategy 2: best 1 indicator per sub-pillar (Level 2).
   - SIMPROF (`clustsig::simprof`) + agglomerative clustering on the full vs. selected-subset data, compared via Rand Index (`rand_ind()` in `scripts/R/functions.R`).

2. **Variable perspective** — does a subset preserve the PCA structure (i.e., "information content") of the full indicator set? Uses the **r_m criterion**:

   ```
   r_m = sqrt( Σ λ_i · r_mi²  /  Σ λ_j )
   ```

   where `λ_i` is the eigenvalue of principal component *i* from the full dataset and `r_mi²` is the R² from regressing PC *i* on the candidate subset's indicators (`PC_cor` / `r_m` in `scripts/R/functions.R`). This is a *subset*-selection criterion, not a per-variable importance score — it measures how well an entire subset reconstructs the full PCA ordination, and has built-in redundancy protection (near-duplicate indicators don't all get selected). This is currently only implemented in `scripts/legacy/04_variable_perspective_*.R` and needs consolidating into a new `scripts/04_variable_perspective.R` per the implementation plan.

Full mathematical derivation and interpretation caveats (e.g., that this treats "information" as PCA variance, which may not match a substantive policy concept) are in [src/method-explanation.md](src/method-explanation.md).

### Data hierarchy (GII / COINr `imeta` structure)

The index has 5 levels, defined in `gii-data/imeta.csv`/`.xlsx` and mirrored in `scripts/R/config.R`:

| Level | Meaning | Count | Example codes |
|---|---|---|---|
| 1 | Indicators | 77 | individual `iCode`s |
| 2 | Sub-pillars | 21 | `SP1.1`–`SP7.3` |
| 3 | Pillars | 7 | `IN.1`–`IN.5`, `OUT.6`, `OUT.7` |
| 4 | Sub-indices | 2 | Inputs, Outputs |
| 5 | Index | 1 | Global Innovation Index |

Unlike the original (Kenneth's) analysis, this project includes **all** pillars/sub-pillars with no spillover exclusions, and analyzes a single year (2025) rather than a time series.

### Code layout

- `scripts/R/config.R` — all parameters/paths/constants (`CONFIG` list): indicator counts, pillar/sub-pillar codes, exhaustive-search group sizes, Fig-2 k-values, SIMPROF params, PCA variance threshold, random seed. Read this first when changing analysis parameters.
- `scripts/R/functions.R` — shared analysis functions: `sim_matrix`, `spearman_dist_cor`, `PC_cor`/`r_m` (and `_sub` variants for within-group analysis), `get_pca_components`, `pillar_indicators`/`subpillar_indicators`, `best_k_within_group` (exhaustive search within a pillar/sub-pillar), `rand_ind` (Rand Index from two `simprof` results), `proportional_allocation`.
- `scripts/R/data_utils.R` — `load_gii_data()` is the single entry point for loading+cleaning data: reads `gii-data/idata.xlsx` + `imeta.xlsx`, subsets to Level-1 indicators, kNN-imputes missing values (`VIM::kNN`, k=5), applies `standardize_with_direction()` (reverses indicators with `Direction == -1` via `max(x) - x`, since GII imeta encodes polarity), then z-scores. Any new script needing data should call this rather than re-implementing loading logic.
- `scripts/00_setup.R` → `01_data_preparation.R` → `02_eda.R` → `03_country_perspective.R` — numbered pipeline; each script prints progress and writes its outputs to `outputs/` or `data/`.
- `scripts/legacy/` — pre-refactor versions of the numbered scripts, including the not-yet-ported variable-perspective logic (`04_variable_perspective_bypillar.R`, `04_variable_perspective_exhausitive_search.R`, `04_variable_perspective_random.R`, `05_indicator_selection_gii.R`, `05_results.R`, `master_script.R`). Treat these as reference/source-of-truth for logic not yet migrated into `scripts/R/functions.R`, not as code to run directly.
- `Example/Matrix_Work.2.R`, `Example/PCA_Work.2.R` — the original ("Kenneth's") reference implementation this project's methodology was standardized against; consult these when reconciling behavior differences.
- `Example/clustsig/`, `Example/clustsig-optimization/` — vendored copies of the CRAN `clustsig` package (provides `simprof`, used for statistically significant hierarchical clustering per Clarke, Somerfield & Gorley 2008). `clustsig-optimization` is a performance-oriented fork/variant; prefer installing `clustsig` from CRAN or GitHub (`remotes::install_github("douglaswhitaker/clustsig")`) rather than editing these vendored copies unless specifically working on the package itself.
- `GII_Analysis_Report.qmd` — Quarto document (PDF output) that consumes files from `outputs/` and `data/` to produce the final report; must be re-rendered after pipeline scripts change their output files.
- `IMPLEMENTATION_PLAN.md` — current refactor plan, including the target file structure, the full list of expected `outputs/*.csv`/`.rds` files and which script should produce each, and the consolidated function list. Consult this when adding new outputs or scripts to keep naming/location consistent.
- `todo-list.md` — task-level status of the two analysis perspectives (marks what's DONE vs. outstanding); check before assuming a described analysis hasn't been implemented.

### Output conventions

Pipeline scripts write intermediate/cleaned data to `data/` and analysis results to `outputs/` (CSVs for tabular results, `.rds` for R objects like `simprof` clustering results). `outputs/fig2_best_subsets/` holds one CSV per k-value tested in the Fig-2 random-subset analysis. File names are pre-specified in `IMPLEMENTATION_PLAN.md`'s output table — match them when adding new analysis outputs so the QMD report can find them.
