# Configuration for GII Index Reduction Analysis
# All magic numbers, paths, and parameters centralized here

CONFIG <- list(
  # File paths
  paths = list(
    data_raw      = "gii-data/idata.xlsx",
    imeta_raw     = "gii-data/imeta.xlsx",
    data_std      = "data/gii_data_standardized.csv",
    data_complete = "data/gii_data_complete_cases.csv",
    imputation    = "data/imputation_report.csv",
    stats_original = "data/descriptive_stats_original.csv",
    stats_complete = "data/descriptive_stats_complete.csv",
    outputs_dir   = "outputs"
  ),

  # Indicator counts for Fig-2 k grid random selection. 
  # Is not 78 as we want to search only for subsets of the full idnicator list
  n_indicators = 77,  

  # Hierarchical structure (from imeta$Level)
  # Level 1: Indicators (78)
  # Level 2: Sub-pillars (21) - SP1.1 through SP7.3
  # Level 3: Pillars (7) - IN.1-IN.5, OUT.6, OUT.7
  # Level 4: Sub-indices (2) - Inputs, Outputs
  # Level 5: Index (1) - Global Innovation Index
  pillars_level3 = c("IN.1", "IN.2", "IN.3", "IN.4", "IN.5", "OUT.6", "OUT.7"),
  subpillars_level2 = c(
    "SP1.1", "SP1.2", "SP1.3",
    "SP2.1", "SP2.2", "SP2.3",
    "SP3.1", "SP3.2", "SP3.3",
    "SP4.1", "SP4.2", "SP4.3",
    "SP5.1", "SP5.2", "SP5.3",
    "SP6.1", "SP6.2", "SP6.3",
    "SP7.1", "SP7.2", "SP7.3"
  ),

  # Exhaustive search strategies
  exhaustive = list(
    strategy1_k_per_pillar = 3,   # 3 vars per pillar (Level 3) -> 7 * 3 = 21 vars
    strategy2_k_per_subpillar = 1,  # 1 var per sub-pillar (Level 2) -> 21 * 1 = 21 vars
    # Strategies 1a/2a: total budget allocated proportionally to group size
    # (largest-remainder method) instead of a flat k per group. 42 (~2/sub-pillar
    # on average) gives sub-pillar allocation real room to vary instead of
    # collapsing to 1-per-group.
    proportional_target = 42,
    # Matched-budget control (Strategy 1b): the same proportional apportionment
    # at 21, the flat strategies' budget. Without it, "flat vs proportional" is
    # confounded with "21 vs 42 indicators" and neither can be attributed.
    # PILLAR LEVEL ONLY. At 21 the pillar apportionment is 2,3,2,3,4,4,3, which
    # genuinely differs from a flat 3-per-pillar. The sub-pillar apportionment
    # at 21 is exactly 1 per sub-pillar - identical to flat Strategy 2 - so
    # there is no sub-pillar counterpart to run; the two rules coincide by
    # construction at that budget, which is why proportional_target is 42.
    proportional_target_matched = 21
  ),

  # Figure 2: Random subset correlations
  fig2 = list(
    k_percentiles = seq(0.1, 0.9, 0.1),  # 10% to 90% of 77 indicators
    n_random_subsets = 1000
  ),

  # SIMPROF clustering parameters
  clustering = list(
    simprof_expected = 1000,
    simprof_simulated = 999,
    method_cluster = "average",
    method_dist = "euclidean",
    alpha = 0.05
  ),

  # PCA parameters
  pca = list(
    variance_threshold = 0.95
  ),

  # Pillar allocation for constrained sampling
  sampling = list(
    subset_sizes = c(10, 15, 20, 25, 30),
    n_random_subsets = 20,
    min_per_pillar = 1
  ),

  # Random seed for reproducibility
  seed = 42
)

# Derived values
CONFIG$fig2$k_values <- round(CONFIG$fig2$k_percentiles * CONFIG$n_indicators)
CONFIG$fig2$k_values <- unique(pmax(1, pmin(CONFIG$n_indicators - 1, CONFIG$fig2$k_values)))
# Add k=42 explicitly: matches the proportional strategies' budget (Strategies 1a/2a),
# so the random-subset baseline can be compared directly at that exact size.
CONFIG$fig2$k_values <- sort(unique(c(CONFIG$fig2$k_values, CONFIG$exhaustive$proportional_target)))

# Print config summary on load
message("GII Analysis Config loaded:")
message("  - Indicators: ", CONFIG$n_indicators)
message("  - Pillars (L3): ", length(CONFIG$pillars_level3))
message("  - Sub-pillars (L2): ", length(CONFIG$subpillars_level2))
message("  - Fig 2 k-values: ", paste(CONFIG$fig2$k_values, collapse = ", "))