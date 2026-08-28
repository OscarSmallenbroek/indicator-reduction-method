# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

This is an R-based statistical analysis project (not software with a build/test pipeline). It develops and applies a method for determining how far the number of indicators in a composite index can be reduced without significantly compromising the index's discriminatory power. The method was previously applied to the EU JRC's Sustainable and Inclusive Well-being Index; this project replicates it on the **Global Innovation Index (GII)** 2025 data (78 Level-1 indicators, 139 countries).

There is no package manager / lockfile (no `renv`, no `.Rproj`). Scripts are run directly with `Rscript` or sourced interactively; missing packages are auto-installed by `scripts/00_setup.R`.

## Running the analysis

Scripts must be run with the repo root as the working directory (they use relative paths like `"scripts/R/config.R"` and `"gii-data/idata.xlsx"`).

```bash
Rscript scripts/00_setup.R
```

Run `00_setup.R` once to install packages. Each numbered script then sources
`scripts/R/config.R`, `functions.R` and `data_utils.R` itself and attaches its own
libraries, so they run standalone under `Rscript` in any order (01 before the report):

```bash
Rscript scripts/01_data_preparation.R
Rscript scripts/02_eda.R
Rscript scripts/03_country_perspective.R
Rscript scripts/03b_country_perspective_proportional.R
Rscript scripts/04_variable_perspective.R
```

`03` and `03b` both reuse `outputs/clustering_full.rds` when it exists rather than
recomputing the full-dataset SIMPROF clustering, which is the pipeline's slowest step by
a wide margin. Run `03` first so `03b` finds it, and delete that file to force a fresh
clustering after the input data or the clustering parameters change. `01` and `02` write
only report inputs — `03`/`03b`/`04` do not read `data/`, they each call
`load_gii_data()` themselves.

To render the report (requires Quarto + a LaTeX/PDF engine, since output format is `pdf`):

```bash
quarto render GII_Analysis_Report.qmd
```

There is no automated test suite; `Example/clustsig-optimization/test_clustsig.R`, `test_runner.R`, and `test_windows_fix.R` are ad hoc test scripts for the vendored `clustsig` package fork only, not for this project's own analysis code.

## Architecture

### Two selection metrics x two allocation rules

Two things vary independently, and conflating them has already caused one
replication failure — keep them separate when reading or changing the code.

**The selection metric** is what we are replicating from the reference project
in `Example/`. Each perspective has its own criterion, and each was implemented
in its own reference script:

| Metric | Perspective | Question | Reference implementation | Our scorer | Score column |
|---|---|---|---|---|---|
| `"distance"` | Country | Do countries end up "the same distance apart"? | `find_best_c`/`find_best_d` → `subset_test` in [Example/Matrix_Work.2.R](Example/Matrix_Work.2.R) | `make_dist_scorer()` | `best_spearman` |
| `"r_m"` | Variable | Does the subset preserve the PCA structure? | `best_k_comp` in [Example/PCA_Work.2.R](Example/PCA_Work.2.R) | `make_rm_scorer()` | `best_rm` |

**The allocation rule** — how many indicators each group gets — is *our own
extension and adaptation*; it has no counterpart in the reference project. Four
rules are defined: Strategy 1 (flat, 3 per pillar), Strategy 2 (flat, 1 per
sub-pillar), and Strategies 1a/2a (a 42-indicator budget apportioned across
pillars/sub-pillars by `largest_remainder_allocation()`).

Any allocation can be searched under either metric. `best_allocation_within_group()`
takes `metric = c("r_m", "distance")` and dispatches to the matching scorer;
`best_k_within_group()` is its flat-k wrapper and forwards `metric` unchanged.
The default is `"r_m"`, so **a call that omits `metric` gets the variable
perspective** — `03`/`03b` must pass `metric = "distance"` explicitly.

Which script uses which:

| Script | Allocation | Metric |
|---|---|---|
| `03_country_perspective.R` | flat (Strategies 1, 2) | `"distance"` |
| `03b_country_perspective_proportional.R` | proportional (1a, 2a) | `"distance"` |
| `04_variable_perspective.R` | all four | `"r_m"` |

So `03`/`03b` and `04` run the same allocations over the same data and differ
only in the metric — that contrast is the point, and the two sets of outputs
are expected to disagree (currently 8/21 overlap for Strategy 1).

> **History:** until 2026-08-28, `03`/`03b` also selected on `r_m` and merely
> *reported* a Spearman correlation on the result, so `03` and `04` emitted
> identical variable lists and the country perspective was never actually
> optimised for its own criterion. A coauthor's independent port of
> `Matrix_Work.2.R` failed to reproduce our Strategy 1 for this reason. If a
> future selection looks suspiciously identical across `03` and `04`, check
> that `metric = "distance"` is still being passed.

### Naming: internal codes vs. report prose

`strategy1`, `2`, `1a`, `2a`, `1b` are **internal identifiers only** - variable
names, output filenames, log messages. They must never appear in the report's
prose, headings, table captions, or column values. A reader has no way to know
what "Strategy 2a" means, and the labels carry no information about what
distinguishes them.

In `GII_Analysis_Report.qmd`, refer to each by what it does:

| Internal code | Prose name |
|---|---|
| Strategy 1 | flat allocation by pillar (3 indicators per pillar) |
| Strategy 2 | flat allocation by sub-pillar (1 indicator per sub-pillar) |
| Strategy 1a | proportional allocation across pillars (42 indicators) |
| Strategy 2a | proportional allocation across sub-pillars (42 indicators) |
| Strategy 1b | proportional allocation across pillars at the flat budget (21 indicators) |

Shorten to "flat allocation by pillar" / "proportional allocation across
pillars" once the budget is clear from context. The `strategy` column in the
Rand Index CSVs stores internal codes (`Strategy1_3perPillar` etc.); the report
maps them to prose names before display - keep that mapping current when adding
a strategy.

### How the search works, and what "local" means

Both metrics score **locally**: a candidate subset is compared against its own
group's sub-PCA (`r_m`) or its own group's distance matrix (`distance`), never
against the whole index. This is faithful to both reference scripts — Kenneth's
own comment in `Matrix_Work.2.R` flags it ("scores only measure relative to the
reduced matrix for that component"). The assembled 21/42-indicator subset is
scored against the full index exactly once, by the calling script, after the
per-group winners are concatenated.

`make_dist_scorer()` ranks the group's distance vector once and reuses it across
combinations (Spearman is Pearson on ranks, and those ranks never change);
the original recomputed `dist()` inside the loop. The results are identical —
verified against a literal transcription of `subset_test` on all 7 pillars.

A group with no more indicators than its budget takes all of them and scores
`1.0` under either metric, matching both originals.

The `r_m` criterion itself:

```
r_m = sqrt( Σ λ_i · r_mi²  /  Σ λ_j )
```

where `λ_i` is the eigenvalue of principal component *i* and `r_mi²` is the R²
from regressing PC *i* on the candidate subset's indicators (`PC_cor` / `r_m` in
`scripts/R/functions.R`). It is a *subset*-selection criterion, not a
per-variable importance score — it measures how well an entire subset
reconstructs the PCA ordination, and has built-in redundancy protection
(near-duplicate indicators don't all get selected). `r_m()` scores every
principal component from a single QR factorisation (`subset_r_squared()`), since
the design matrix is identical across components; the same function serves both
the full-dataset PCA and a group's own sub-PCA.

Full mathematical derivation and interpretation caveats (e.g., that `r_m` treats
"information" as PCA variance, which may not match a substantive policy concept)
are in [src/method-explanation.md](src/method-explanation.md).

### Other country-perspective outputs

Beyond the exhaustive search, `03` also produces:

- Spearman rank correlation between full vs. subset Euclidean distance matrices
  over 1000 random subsets at k = 10%–90% of indicators, the baseline the
  exhaustive strategies are judged against (`outputs/fig2_rank_correlations.csv`,
  `outputs/fig2_best_subsets/k_*.csv`).
- SIMPROF (`clustsig::simprof`) + agglomerative clustering on the full vs.
  selected-subset data, compared via Rand Index (`rand_ind()` in
  `scripts/R/functions.R`).

### Data hierarchy (GII / COINr `imeta` structure)

The index has 5 levels, defined in `gii-data/imeta.csv`/`.xlsx` and mirrored in `scripts/R/config.R`:

| Level | Meaning | Count | Example codes |
|---|---|---|---|
| 1 | Indicators | 78 | individual `iCode`s |
| 2 | Sub-pillars | 21 | `SP1.1`–`SP7.3` |
| 3 | Pillars | 7 | `IN.1`–`IN.5`, `OUT.6`, `OUT.7` |
| 4 | Sub-indices | 2 | Inputs, Outputs |
| 5 | Index | 1 | Global Innovation Index |

Unlike the original (Kenneth's) analysis, this project includes **all** pillars/sub-pillars with no spillover exclusions, and analyzes a single year (2025) rather than a time series.

### Code layout

- `scripts/R/config.R` — all parameters/paths/constants (`CONFIG` list): pillar/sub-pillar codes, exhaustive-search group sizes, Fig-2 k-values, SIMPROF params, PCA variance threshold, random seed. Read this first when changing analysis parameters.
  - **`CONFIG$n_indicators = 77` is a search parameter, not a dataset count.** The data has 78 Level-1 indicators; 77 is the basis for the Fig-2 k grid (`k_percentiles * n_indicators`, capped at `n_indicators - 1`) so that every sampled k is a *strict subset* of the full list. Do not "correct" it to 78 — that is not what it measures. The actual indicator count is always derived from the data (`ncol(indicator_data)`) or metadata (`group_indicator_counts()`), never from this constant.
- `scripts/R/functions.R` — shared analysis functions: `sim_matrix`, `spearman_dist_cor`, `PC_cor`/`subset_r_squared`/`r_m` (one implementation, used for both the full and per-group PCA), `get_pca_components`, `group_variables`/`group_indicator_counts`, `make_dist_scorer`/`make_rm_scorer` (the two selection metrics), `best_allocation_within_group` (exhaustive search within each pillar/sub-pillar given a per-group budget and a `metric`) with `best_k_within_group` as its flat-k wrapper, `largest_remainder_allocation`, `rand_ind`/`simprof_partition` (Rand Index from two `simprof` results — `simprof_partition()` derives the unit count from the clustering and asserts every unit was assigned; all 139 countries are scored, with no positional exclusions).
- `scripts/R/simprof_fast.R` — `use_fast_simprof()` patches four `clustsig` internals (`genSimilarityProfile`, `genProfile`, `computeAverage`, `tsComparison`) in place with faster rewrites. A similarity profile is just a *sorted* distance vector, so it is built with a Gram-matrix identity and `sort()` instead of `dist()` + `rank()` + `order()`; `genProfile` drops the list of every permuted data copy that the package builds and never reads, keeps profiles in a matrix rather than a list of two-row matrices, and reduces the expected profiles to a running sum since only their mean is ever used. `03`/`03b` call it right after sourcing the shared files. SIMPROF is by far the slowest step in the pipeline, so this is where run time goes; the exhaustive searches take ~2 s and Fig 2 ~20 s. Results are unchanged: nodes of 12 or fewer countries keep using `dist()` (measurably faster at that size, and exact where the permutation test is degenerate - at a 2-country node no column permutation can change the single pairwise distance, so the test statistic is exactly zero and rounding noise would decide the tie), and above that the Gram identity's ~1e-14 departure from `dist()` was verified not to move any p-value. Non-Euclidean distances fall through to the package implementation.
- `scripts/R/data_utils.R` — `load_gii_data()` is the single entry point for loading+cleaning data: reads `gii-data/idata.xlsx` + `imeta.xlsx`, subsets to Level-1 indicators, kNN-imputes missing values (`VIM::kNN`, k=5), applies `standardize_with_direction()` (reverses indicators with `Direction == -1` via `max(x) - x`, since GII imeta encodes polarity), then z-scores. Any new script needing data should call this rather than re-implementing loading logic.
- `scripts/00_setup.R` (one-off package install) then `01_data_preparation.R`, `02_eda.R`, `03_country_perspective.R`, `03b_country_perspective_proportional.R`, `04_variable_perspective.R` — each sources `scripts/R/{config,functions,data_utils}.R` itself and runs standalone under `Rscript`; each reports progress via `message()` and writes to `outputs/` or `data/`.
- `scripts/legacy/` — pre-refactor versions of the numbered scripts. The variable-perspective logic they held (`04_variable_perspective_{bypillar,exhausitive_search,random}.R`) is now consolidated in `scripts/04_variable_perspective.R`; `05_indicator_selection_gii.R`, `05_results.R` and `master_script.R` were never ported. Treat all of these as historical reference, not as code to run.
- `Example/Matrix_Work.2.R`, `Example/PCA_Work.2.R` — the original ("Kenneth's") reference implementation this project's methodology was standardized against; consult these when reconciling behavior differences. `Matrix_Work.2.R` is the source of the **distance** metric (`find_best_c`/`find_best_d`, scoring via `subset_test`), `PCA_Work.2.R` the source of the **r_m** metric (`best_k_comp`). Only the metrics are replicated from these — the allocation rules are ours. Note both originals hardcode details specific to that dataset (27 units, dropping EU27 at position 11, excluding Spillovers) that deliberately do **not** carry over here.
- `Example/clustsig/`, `Example/clustsig-optimization/` — vendored copies of the CRAN `clustsig` package (provides `simprof`, used for statistically significant hierarchical clustering per Clarke, Somerfield & Gorley 2008). `clustsig-optimization` is a performance-oriented fork/variant; prefer installing `clustsig` from CRAN or GitHub (`remotes::install_github("douglaswhitaker/clustsig")`) rather than editing these vendored copies unless specifically working on the package itself.
- `GII_Analysis_Report.qmd` — Quarto document (PDF output) that consumes files from `outputs/` and `data/` to produce the final report; must be re-rendered after pipeline scripts change their output files.
- `todo-list.md` — task-level status of the two analysis perspectives (marks what's DONE vs. outstanding); check before assuming a described analysis hasn't been implemented.

### Output conventions

The per-group score column in an exhaustive-search CSV is named for the metric that produced it: `best_rm` (r_m) or `best_spearman` (distance). `GII_Analysis_Report.qmd`'s `format_group_table()` maps both to display headers, so adding a metric means adding a mapping there too. The separate `correlation` / `r_m` column in those files is the *assembled* subset's global score, repeated on every row — not the per-group score.

Pipeline scripts write intermediate/cleaned data to `data/` and analysis results to `outputs/` (CSVs for tabular results, `.rds` for R objects like `simprof` clustering results). `outputs/fig2_best_subsets/` holds one CSV per k-value tested in the Fig-2 random-subset analysis. There is no separate output manifest: the authoritative list is the `write.csv()`/`saveRDS()` calls in the pipeline scripts and the matching `read.csv()`/`readRDS()` calls in the `load-results` chunk of `GII_Analysis_Report.qmd`. When adding an output, add both ends — an output nothing reads, or a read of a file nothing writes, is the failure mode to avoid (the report's `data/imputation_report.csv` was orphaned this way).
