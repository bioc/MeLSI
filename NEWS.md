# MeLSI 1.0.1

- Release backport of the cumulative performance work developed in the 1.1.x
  devel series, brought to the Bioconductor release branch. Two fused C++ (Rcpp)
  kernels replace the R hot paths: the PERMANOVA pseudo-F statistic and the
  weak-learner metric optimizer (`src/melsi_fstat.cpp`). Together with the R-level
  optimizations (vectorized prefiltering, direct F-statistic in place of
  `vegan::adonis2`, diagonal-metric column scaling, diagonal-vector ensemble
  storage), this speeds up analysis substantially (up to ~100x on high-dimensional
  inputs) with no change to the statistical method.
- Reproducibility: results are bit-faithful to 1.0.0. The C++ optimizer draws its
  class pairs with R's own index sampler (`R_unif_index`), consuming the global RNG
  stream identically to the previous `sample()`-based loop; F-statistics reproduce
  to floating-point precision with no p-value changes. Validated against the
  published simulation tables (Type I error, scalability, parameter sensitivity,
  feature correlation).

# MeLSI 1.1.7

- Performance: the weak-learner metric optimizer (`optimize_weak_learner_robust`,
  the hot path of ensemble metric learning) is now a C++ kernel
  (`melsi_opt_weak_learner` in `src/melsi_fstat.cpp`, via Rcpp). Profiling showed
  this gradient-descent loop accounted for ~83% of total runtime; the C++ port
  runs it ~4.6x faster, for roughly ~2x end-to-end speedup, composing with the
  1.1.6 F-statistic kernel and `BPPARAM` parallelism.
- Reproducibility: results are bit-faithful to the R implementation. The C++
  loop draws its within/between-class pairs with R's own index sampler
  (`R_unif_index`), so it consumes the global RNG stream identically to the
  previous `sample()`-based loop. End-to-end F-statistics reproduce 1.1.6 to
  floating-point precision (~1e-15 relative) with no p-value changes.

# MeLSI 1.1.6

- Performance: the PERMANOVA pseudo-F statistic is now computed by a fused C++
  kernel (`src/melsi_fstat.cpp`, via Rcpp) that accumulates squared distances in
  a single pass over sample pairs, instead of materialising and squaring the
  full n x n distance matrix in R. This speeds up the F-statistic by ~4x and the
  overall analysis by roughly 1.3-2.2x (most for two-group / pairwise analyses),
  and composes multiplicatively with `BPPARAM` parallelism.
- Numerical note: results are statistically equivalent but not bit-identical to
  1.1.5. Computing squared distances directly (no sqrt round-trip) and summing in
  a different order changes F-statistics at the ~1e-12 relative level. Across
  extensive testing no permutation p-value changed; in rare boundary cases a
  p-value could differ by a single permutation step (1/(n_perms+1)). Set a
  larger `n_perms` if exact reproducibility of borderline p-values matters.

# MeLSI 1.1.5

- Dependencies: `vegan` moved from Imports to Suggests. The package computes the
  PERMANOVA F-statistic directly and no longer calls `vegan::adonis2()` at
  runtime, so `vegan` (and its transitive dependencies) is no longer installed
  for users; it is needed only to build the comparison section of the vignette.
- Cleanup: removed unused `stats::aov`, `stats::as.dist`, and `stats::t.test`
  imports and an unreachable t-test pre-filtering branch in the metric learner.
  Pre-filtering is handled by `apply_conservative_prefiltering()` as before.
- UX: permutation progress is now shown with a single progress bar instead of one
  message per permutation.
- UX: `melsi()` gains an optional `seed` argument for reproducible runs and now
  returns an object of class `"melsi"` with a `print()` method summarising the
  F-statistic, p-value, and top features. The returned object is still a list,
  so existing `$` access is unchanged.
- All statistical results (F-statistics, p-values, feature weights) are
  numerically identical to 1.1.4; these changes do not alter the method.

# MeLSI 1.1.4

- Documentation: document the `BPPARAM` argument to `melsi()` for parallel
  permutation testing. The argument has been available since the move to a
  `BiocParallel` backend; this release adds a worked example to the README and
  vignette and surfaces it in NEWS. Pass any `BiocParallelParam` object (for
  example `MulticoreParam(workers = 8)`) to distribute the permutation loop
  across cores; the default (`NULL`) runs permutations sequentially.

# MeLSI 1.1.3

- Performance: additional ~22% speedup via diagonal vector storage in ensemble
  aggregation (eliminates B x p^2 matrix allocation per permutation) and
  rejection sampling in gradient optimizer (eliminates setdiff() allocation
  per iteration). Combined with 1.1.2, total speedup is ~2.9x over 1.1.1.

# MeLSI 1.1.2

- Performance: 2.3x faster permutation testing via vectorized prefiltering, direct
  PERMANOVA F-statistic (replacing adonis2 overhead), and diagonal metric matrix
  optimization (replacing O(p^3) eigen decomposition with O(p) column scaling).
  Results are numerically identical; p-values unchanged.
- Suppress spurious NaN warnings from log2 fold-change on CLR-transformed data.

# MeLSI 0.99.9

- Lower minimum R dependency to R >= 4.5.0 to match Bioconductor 3.23 build system

# MeLSI 0.99.8

- Bump version to retrigger Bioconductor build report after CI configuration updates

# MeLSI 0.99.7

- Switch CI to grimbough/bioc-actions to match the Bioconductor Build System environment

# MeLSI 0.99.6

- Update CI to use R-devel and Bioconductor 3.23 targeting R >= 4.6.0

# MeLSI 0.99.5

- Bump version to retrigger build after dependency updates

# MeLSI 0.99.4

- Update minimum R dependency to R >= 4.6.0

# MeLSI 0.99.3

- Address remaining Bioconductor reviewer checklist items for package acceptance
- Enable progress messages for pairwise comparisons in multi-group analysis

# MeLSI 0.99.2

- Optimize vignette runtime to meet Bioconductor < 15 minute requirement
- Reduce example dataset size and permutation counts for faster vignette build

# MeLSI 0.99.1

- Resolve Bioconductor pre-acceptance review issues: warnings, license NOTE, non-standard
  directories, and build artifacts
- Add Matrix to dependencies

# MeLSI 0.99.0

## Initial Bioconductor Submission

- Initial submission of MeLSI (Metric Learning for Statistical Inference) to Bioconductor
- Novel machine learning method for microbiome beta diversity analysis
- Learns optimal distance metrics to improve statistical power in detecting group differences
- Comprehensive validation against standard methods (Bray-Curtis, Euclidean, Jaccard)
- Robust ensemble learning approach with conservative pre-filtering
- Validated on real microbiome datasets with proper Type I error control
- Provides feature importance weights for biological interpretability
- Includes helper functions for CLR transformation and visualization (VIP plots, PCoA)
- Full integration with Bioconductor ecosystem (phyloseq, microbiome packages)
