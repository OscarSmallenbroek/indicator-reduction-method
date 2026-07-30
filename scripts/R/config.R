# Configuration for GII Index Reduction Analysis
# All magic numbers, paths, and parameters centralized here

CONFIG <- list(
  # File paths
  paths = list(
    data_raw      = "gii-data/idata.xlsx",
    imeta_raw     = "gii-data/imeta.csv",
    data_std      = "data/gii_data_standardized.csv",
    imeta_clean   = "data/imeta_clean.csv",
    outputs_dir   = "outputs"
  ),

  # Indicator counts
  n_indicators = 77,  # Level 1 indicators from imeta

  # Hierarchical structure (from imeta$Level)
  # Level 1: Indicators (77)
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
    proportional_target = 42
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

# Print config summary on load
message("GII Analysis Config loaded:")
message("  - Indicators: ", CONFIG$n_indicators)
message("  - Pillars (L3): ", length(CONFIG$pillars_level3))
message("  - Sub-pillars (L2): ", length(CONFIG$subpillars_level2))
message("  - Fig 2 k-values: ", paste(CONFIG$fig2$k_values, collapse = ", "))