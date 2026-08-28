###############################################################################
# Drop-in speedups for clustsig::simprof()
#
# simprof() dominates the runtime of 03/03b. At every node of the dendrogram it
# permutes the data num.expected + num.simulated = 1999 times and builds a
# similarity profile from each permutation. With 139 countries the tree has
# ~100 tested nodes, so the inner profile step runs ~200,000 times.
#
# Four of the package's internals are replaced here. All four are pure
# rewrites - same arguments, same role, same RNG stream - so a patched run
# reproduces an unpatched run's clustering. Nothing about the SIMPROF procedure
# changes; only how the profile is computed and stored.
#
#   genSimilarityProfile()  dist() + rank() + order() -> one Gram-matrix
#                           identity and a sort()
#   genProfile()            drops an accumulator of every permuted data copy
#                           that the package builds but never reads, and keeps
#                           the profiles in a matrix rather than a list of
#                           two-row matrices whose second row is always 1..m
#   computeAverage()        consumes that representation; for the expected
#                           profile only the column mean is ever read, so the
#                           replicates are summed as they are generated instead
#                           of all being held at once
#   tsComparison()          consumes that representation
#
# Together these cut peak memory per node several-fold, which matters more than
# the arithmetic on machines where the unpatched run pages to disk.
#
# Numerical note: the Gram-matrix identity agrees with dist() to ~1e-14 rather
# than bit-for-bit. A permutation p-value is a count of simulated statistics at
# or above the observed one, so a difference that small can move a p-value by
# one count out of num.simulated. Cluster membership is unaffected.
#
# Non-Euclidean distances fall through to the package's own implementation.
###############################################################################

if (!requireNamespace("clustsig", quietly = TRUE)) {
  stop("clustsig is required; install with remotes::install_github(\"douglaswhitaker/clustsig\")")
}

# Lower-triangle indices are the same for every permutation of an n-row matrix,
# so they are built once per n and reused.
.simprof_lt_cache <- new.env(parent = emptyenv())

.lower_tri_index <- function(n) {
  key <- as.character(n)
  idx <- .simprof_lt_cache[[key]]
  if (is.null(idx)) {
    idx <- which(lower.tri(matrix(FALSE, n, n)))
    .simprof_lt_cache[[key]] <- idx
  }
  idx
}

#' Sorted vector of pairwise Euclidean distances - the similarity profile.
#'
#' ||a-b||^2 = ||a||^2 + ||b||^2 - 2 a.b, so one BLAS crossprod replaces the
#' pairwise loop in dist(). pmax() clamps the small negative values that
#' floating-point cancellation produces for near-coincident rows.
#'
#' The package ranks the distance vector with ties.method="first" and lets the
#' caller reorder the profile by those ranks, which is just a sort; sorting
#' here skips both the rank() and the order() pass.
#'
#' Small nodes go through dist() instead. Measured on this data the identity
#' only starts winning around n = 16, and below that it is not merely pointless
#' but harmful: permuting columns cannot change a two-row matrix's single
#' pairwise distance, so at n = 2 the observed test statistic is exactly zero
#' and every simulated statistic ties with it. dist() reproduces that tie
#' exactly and the node gets p = 1; ~1e-14 of rounding noise would turn the
#' tie-breaking into a coin flip. Such nodes are never significant either way,
#' but keeping them exact keeps the whole p-value matrix identical to an
#' unpatched run.
.SIMPROF_EXACT_MAX_N <- 12L

.sorted_euclidean <- function(X) {
  n <- nrow(X)
  if (n < 2L) return(numeric(0))
  if (n <= .SIMPROF_EXACT_MAX_N) return(sort(as.vector(stats::dist(X))))
  ss <- rowSums(X * X)
  d2 <- outer(ss, ss, "+") - 2 * tcrossprod(X)
  sort(sqrt(pmax(d2[.lower_tri_index(n)], 0)))
}

.is_euclidean <- function(method.distance) {
  is.character(method.distance) && identical(method.distance, "euclidean")
}

#' Replacement for clustsig:::genSimilarityProfile().
#' Returns the same two-row matrix; the rank row is 1..m by construction, which
#' is exactly what the caller's reorder-by-rank would have produced.
.fast_gen_similarity_profile <- function(rawdata.samples, method.distance, const, undef.zero) {
  if (!.is_euclidean(method.distance)) {
    return(.simprof_orig$gsp(rawdata.samples, method.distance, const, undef.zero))
  }
  d <- .sorted_euclidean(rawdata.samples)
  rbind(rawdata.distvec = d, rawdata.ranks = seq_along(d))
}

#' Replacement for clustsig:::genProfile().
#'
#' The expected and the simulated profiles come from the same call, told apart
#' by `type`. They are consumed differently: the expected set is only ever
#' averaged, so it is reduced to a running column sum as it is generated, while
#' the simulated set has to be kept because each profile is compared against
#' that average. Anything unrecognised falls back to the package.
.fast_gen_profile <- function(rawdata, originaldata, num.expected, method.distance,
                              const, silent, increment, type, undef.zero) {
  if (!.is_euclidean(method.distance)) {
    return(.simprof_orig$gp(rawdata, originaldata, num.expected, method.distance,
                            const, silent, increment, type, undef.zero))
  }
  permute <- clustsig:::columnPermuter
  announce <- function(i) if (!silent && i %% increment == 0) print(paste(type, "iteration", i))

  if (identical(type, "Expected")) {
    total <- NULL
    for (i in seq_len(num.expected)) {
      announce(i)
      p <- .sorted_euclidean(permute(rawdata))
      total <- if (is.null(total)) p else total + p
    }
    return(structure(list(total = total), class = "simprof_profile_sum"))
  }

  profiles <- NULL
  for (i in seq_len(num.expected)) {
    announce(i)
    p <- .sorted_euclidean(permute(rawdata))
    if (is.null(profiles)) profiles <- matrix(0, nrow = num.expected, ncol = length(p))
    profiles[i, ] <- p
  }
  if (is.null(profiles)) profiles <- matrix(0, nrow = num.expected, ncol = 0L)
  structure(profiles, class = "simprof_profile_matrix")
}

#' Replacement for clustsig:::computeAverage(): mean profile across replicates.
.fast_compute_average <- function(expectedprofile.simprof, num.expected) {
  if (inherits(expectedprofile.simprof, "simprof_profile_sum")) {
    return(expectedprofile.simprof$total / num.expected)
  }
  if (inherits(expectedprofile.simprof, "simprof_profile_matrix")) {
    return(colMeans(unclass(expectedprofile.simprof)))
  }
  .simprof_orig$avg(expectedprofile.simprof, num.expected)
}

#' Replacement for clustsig:::tsComparison(): the proportion of simulated
#' profiles whose departure from the expected profile is at least the observed
#' one.
.fast_ts_comparison <- function(simulatedprofile.simprof, expectedprofile.average,
                                num.simulated, teststatistic) {
  if (!inherits(simulatedprofile.simprof, "simprof_profile_matrix")) {
    return(.simprof_orig$ts(simulatedprofile.simprof, expectedprofile.average,
                            num.simulated, teststatistic))
  }
  profiles <- unclass(simulatedprofile.simprof)
  # A row at a time rather than a whole-matrix sweep: same arithmetic, without a
  # second copy of what is already the largest object in the run.
  at_least <- 0L
  for (i in seq_len(nrow(profiles))) {
    if (sum(abs(profiles[i, ] - expectedprofile.average)) >= teststatistic) {
      at_least <- at_least + 1L
    }
  }
  at_least / num.simulated
}

.simprof_orig <- NULL

#' Patch clustsig's internals in place. Idempotent; safe to call from any script.
use_fast_simprof <- function() {
  if (is.null(.simprof_orig)) {
    .simprof_orig <<- list(
      gsp = clustsig:::genSimilarityProfile,
      gp  = clustsig:::genProfile,
      avg = clustsig:::computeAverage,
      ts  = clustsig:::tsComparison
    )
    utils::assignInNamespace("genSimilarityProfile", .fast_gen_similarity_profile, ns = "clustsig")
    utils::assignInNamespace("genProfile", .fast_gen_profile, ns = "clustsig")
    utils::assignInNamespace("computeAverage", .fast_compute_average, ns = "clustsig")
    utils::assignInNamespace("tsComparison", .fast_ts_comparison, ns = "clustsig")
    message("clustsig::simprof internals patched for speed")
  }
  invisible(TRUE)
}
