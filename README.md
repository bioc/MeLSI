# MeLSI: Metric Learning for Statistical Inference in Microbiome Analysis

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R](https://img.shields.io/badge/R-4.0+-blue.svg)](https://www.r-project.org/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17714848.svg)](https://doi.org/10.5281/zenodo.17714848)
[![BioC devel](https://bioconductor.org/shields/build/devel/bioc/MeLSI.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/MeLSI/)
[![BioC release](https://bioconductor.org/shields/availability/release/MeLSI.svg)](https://bioconductor.org/packages/MeLSI/)
[![Downloads](https://bioconductor.org/shields/downloads/release/MeLSI.svg)](https://bioconductor.org/packages/MeLSI/)

## Overview

**MeLSI** (Metric Learning for Statistical Inference) is a machine learning framework for microbiome beta diversity analysis that solves a fundamental problem: **traditional distance metrics are fixed and generic, missing biologically meaningful patterns in your data**. MeLSI learns adaptive distance metrics optimized for your specific dataset while maintaining rigorous statistical validity through permutation testing.

### The Problem

Current microbiome beta diversity analysis relies on fixed distance metrics like Bray-Curtis, Euclidean, or Jaccard that treat all microbial taxa equally. This "one-size-fits-all" approach often misses subtle but biologically important differences between groups.

### How MeLSI Works

MeLSI uses an innovative ensemble approach to learn optimal distance metrics for community composition analysis:

1. **Conservative Pre-filtering**: Selects taxa with highest variance using variance-based filtering, keeping 70% of features to maintain statistical power while reducing noise
2. **Bootstrap Sampling**: Creates multiple training sets by resampling your data
3. **Feature Subsampling**: Each learner uses a random subset of features (80%) to prevent overfitting
4. **Gradient Optimization**: Learns optimal weights for each feature subset using gradient descent
5. **Ensemble Averaging**: Combines 30 weak learners, weighting better performers more heavily
6. **Robust Distance Calculation**: Ensures numerical stability with eigenvalue decomposition
7. **Permutation Testing**: Validates significance using null distributions from permuted data (200 permutations for reliable p-values)

### Performance Results

Instead of using the same distance formula for every study, MeLSI learns what matters most for **your specific data**, resulting in:

- ✅ **Proper Type I Error Control** - No false positives on null data (p = 0.607 and 0.224)
- ✅ **Competitive Statistical Power** - Context-dependent performance with interpretable feature weights
- ✅ **Biological Interpretability** - Feature importance weights reveal which taxa drive group differences
- ✅ **Real-World Validation** - 9.1% improvement on Atlas1006, significant detection on DietSwap where traditional methods were marginal

## Requirements

### System Requirements
- **R version**: 4.0.0 or higher
- **Operating System**: Windows, macOS, or Linux

### Required R Packages
MeLSI depends on the following R packages (automatically installed):

**Core Dependencies:**
- `vegan` - Community ecology package for PERMANOVA calculations
- `ggplot2` - Graphics and visualization
- `stats` - Statistical functions (base R)
- `utils` - Utility functions (base R)

**Optional Dependencies:**
- `phyloseq` - Microbiome data structures (via BiocManager)
- `microbiome` - Microbiome analysis tools (via BiocManager)
- `BiocParallel` - For parallel permutation testing (see [Parallel permutation testing](#parallel-permutation-testing))

## Installation

MeLSI can be installed from Bioconductor:

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("MeLSI")
```

After installation, load the package:

```r
library(MeLSI)
```

## Running MeLSI

### Key Outputs

MeLSI automatically adapts to your data:
- **2 groups**: Standard pairwise analysis with feature importance
- **3+ groups**: Comprehensive omnibus and pairwise analysis with multiple testing correction

Results include:
- **F-statistic**: How well groups are separated (higher = better)
- **P-value**: Statistical significance (permutation-based, more reliable)
- **Feature importance weights**: Which taxa matter most - automatically displayed as VIP charts
- **Directionality information**: Indicates which group has higher abundance for each important taxon
- **Log2 fold-change values**: Quantitative differences between groups
- **Top features**: Displayed in console with their learned weights and directionality
- **Multiple testing correction**: Automatic FDR correction for multi-group comparisons
- **Diagnostics**: Quality metrics to ensure reliable results

### Input data

MeLSI requires microbiome count data formatted as:

**Feature abundance data frame**
- Formatted with features (taxa) as columns and samples as rows
- Contains raw count data or relative abundances
- Can include zeros (they will be handled appropriately)
- This can be a filepath to a tab-delimited file

**Group labels**
- Vector of group assignments for each sample
- Should be factors with meaningful group names
- Must match the order of samples in the abundance data

### Output structure

MeLSI returns a results object containing:

**Statistical Results**
- F-statistic: Test statistic for group differences
- P-value: Permutation-based p-value
- Feature weights: Learned importance weights for each taxon
- Directionality: Which group has higher abundance for each taxon
- Distance matrix: Learned distance matrix for PCoA plotting
- Null distribution: F-statistics from permuted data
- Diagnostics: Method performance metrics

## Quick Start (3 lines of code!)

The easiest way to use MeLSI:

```r
library(MeLSI)

# Your microbiome data: X = counts (samples × taxa), y = group labels
X_clr <- clr_transform(X)  # CLR transformation (recommended)
results <- melsi(X_clr, y)   # Run analysis - VIP plot auto-generated!
plot_pcoa(results, X_clr, y)  # PCoA plot
plot_vip(results)  # VIP plot with directionality (default)
plot_vip(results, directionality = FALSE)  # VIP plot without directionality

# That's it! Results include F-statistic, p-value, and feature importance
```

## Parallel permutation testing

The permutation test is the most expensive part of a MeLSI run, and it is
embarrassingly parallel. Pass a `BiocParallelParam` object via the `BPPARAM`
argument to distribute permutations across cores:

```r
library(MeLSI)
library(BiocParallel)

# Use 8 worker cores for the permutation loop (macOS/Linux)
results <- melsi(X_clr, y, n_perms = 999,
                 BPPARAM = MulticoreParam(workers = 8))

# On Windows, use SnowParam instead:
# results <- melsi(X_clr, y, n_perms = 999, BPPARAM = SnowParam(workers = 8))
```

When `BPPARAM` is left at its default (`NULL`), permutations run sequentially.
This works for omnibus, pairwise, and multi-group (with correction) analyses.

## Detailed Examples

### Pairwise Analysis (2 groups)
```r
# Load MeLSI package
library(MeLSI)

# Generate synthetic microbiome data with realistic taxonomic names
set.seed(42)
n_samples <- 60
n_taxa <- 50

# Create realistic bacterial species names
realistic_taxa <- c(
    "Bacteroides_vulgatus", "Bacteroides_thetaiotaomicron", "Bacteroides_fragilis",
    "Faecalibacterium_prausnitzii", "Ruminococcus_bromii", "Ruminococcus_gnavus",
    "Bifidobacterium_longum", "Lactobacillus_acidophilus", "Akkermansia_muciniphila",
    "Prevotella_copri", "Parabacteroides_distasonis", "Alistipes_putredinis",
    "Roseburia_intestinalis", "Coprococcus_comes", "Blautia_wexlerae",
    "Collinsella_aerofaciens", "Eggerthella_lenta", "Fusobacterium_nucleatum",
    "Veillonella_parvula", "Streptococcus_salivarius", "Enterococcus_faecalis",
    "Staphylococcus_epidermidis", "Propionibacterium_acnes", "Corynebacterium_striatum",
    "Escherichia_coli", "Klebsiella_pneumoniae", "Enterobacter_cloacae",
    "Clostridium_innocuum", "Clostridium_leptum", "Methanobrevibacter_smithii",
    "Desulfovibrio_piger", "Succinivibrio_dextrinosolvens", "Micrococcus_luteus",
    "Sarcina_ventriculi", "Slackia_equolifaciens", "Barnesiella_intestinihominis",
    "Odoribacter_splanchnicus", "Porphyromonas_gingivalis", "Bacteroides_ovatus",
    "Bacteroides_uniformis", "Bacteroides_caccae", "Bifidobacterium_adolescentis",
    "Bifidobacterium_bifidum", "Lactobacillus_casei", "Lactobacillus_plantarum",
    "Lactobacillus_rhamnosus", "Lactobacillus_reuteri", "Prevotella_oralis",
    "Dorea_longicatena", "Ruminococcus_bromii"
)

# Create base microbiome data (log-normal distribution)
X <- matrix(rlnorm(n_samples * n_taxa, meanlog = 2, sdlog = 1), 
            nrow = n_samples, ncol = n_taxa)

# IMPORTANT: Set realistic column names
colnames(X) <- realistic_taxa[1:n_taxa]

# Create group labels
y <- c(rep("Control", n_samples/2), rep("Treatment", n_samples/2))

# Add signal to first 10 taxa in Treatment group (simulate differential abundance)
X[31:60, 1:10] <- X[31:60, 1:10] * 1.5

# CLR transformation (recommended for microbiome data)
X_clr <- clr_transform(X)

# Run MeLSI analysis - VIP plot with directionality is auto-generated!
results <- melsi(X_clr, y, n_perms = 200, B = 30, show_progress = TRUE)

# Display results summary
cat("MeLSI Results:\n")
cat(sprintf("F-statistic: %.4f\n", results$F_observed))
cat(sprintf("P-value: %.4f\n", results$p_value))
cat(sprintf("Significant: %s\n", ifelse(results$p_value < 0.05, "Yes", "No")))

# Create PCoA plot
plot_pcoa(results, X_clr, y)

# Re-plot VIP with more features (directionality on by default)
plot_vip(results, top_n = 20)

# Plot VIP without directionality coloring
plot_vip(results, top_n = 20, directionality = FALSE)

# Access directionality information
cat("\nTop taxa with directionality:\n")
if (!is.null(results$directionality)) {
    dir_df <- data.frame(
        Taxon = names(results$directionality),
        Higher_in = results$directionality,
        Weight = results$feature_weights[names(results$directionality)],
        Log2FC = results$log2_fold_change
    )
    print(head(dir_df[order(-dir_df$Weight), ], 10))
}
```

### Multi-Group Analysis (3+ groups)
```r
# Create 4-group data with realistic taxonomic names
set.seed(42)
n_samples_per_group <- 30
n_taxa <- 50

# Use same realistic taxa as pairwise example
realistic_taxa <- c(
    "Bacteroides_vulgatus", "Bacteroides_thetaiotaomicron", "Bacteroides_fragilis",
    "Faecalibacterium_prausnitzii", "Ruminococcus_bromii", "Ruminococcus_gnavus",
    "Bifidobacterium_longum", "Lactobacillus_acidophilus", "Akkermansia_muciniphila",
    "Prevotella_copri", "Parabacteroides_distasonis", "Alistipes_putredinis",
    "Roseburia_intestinalis", "Coprococcus_comes", "Blautia_wexlerae",
    "Collinsella_aerofaciens", "Eggerthella_lenta", "Fusobacterium_nucleatum",
    "Veillonella_parvula", "Streptococcus_salivarius", "Enterococcus_faecalis",
    "Staphylococcus_epidermidis", "Propionibacterium_acnes", "Corynebacterium_striatum",
    "Escherichia_coli", "Klebsiella_pneumoniae", "Enterobacter_cloacae",
    "Clostridium_innocuum", "Clostridium_leptum", "Methanobrevibacter_smithii",
    "Desulfovibrio_piger", "Succinivibrio_dextrinosolvens", "Micrococcus_luteus",
    "Sarcina_ventriculi", "Slackia_equolifaciens", "Barnesiella_intestinihominis",
    "Odoribacter_splanchnicus", "Porphyromonas_gingivalis", "Bacteroides_ovatus",
    "Bacteroides_uniformis", "Bacteroides_caccae", "Bifidobacterium_adolescentis",
    "Bifidobacterium_bifidum", "Lactobacillus_casei", "Lactobacillus_plantarum",
    "Lactobacillus_rhamnosus", "Lactobacillus_reuteri", "Prevotella_oralis",
    "Dorea_longicatena", "Ruminococcus_bromii"
)

X <- matrix(rlnorm(n_samples_per_group * 4 * n_taxa, meanlog = 2, sdlog = 1), 
            nrow = n_samples_per_group * 4, ncol = n_taxa)

# IMPORTANT: Set realistic column names
colnames(X) <- realistic_taxa[1:n_taxa]

y <- rep(c("Control", "Mild", "Moderate", "Severe"), each = n_samples_per_group)

# Add progressive signal to simulate disease severity
X[31:60, 1:10] <- X[31:60, 1:10] * 1.4    # Mild
X[61:90, 1:10] <- X[61:90, 1:10] * 1.7    # Moderate  
X[91:120, 1:10] <- X[91:120, 1:10] * 2.0  # Severe

# CLR transformation
X_clr <- clr_transform(X)

# Run MeLSI analysis (automatically detects 4 groups and runs both omnibus and pairwise)
results <- melsi(X_clr, y, n_perms = 200, B = 30, show_progress = TRUE)

# Access omnibus results
cat("Omnibus F-statistic:", results$omnibus$F_observed, "\n")
cat("Omnibus p-value:", results$omnibus$p_value, "\n")

# Show global feature importance
cat("\nTop 5 globally important taxa:\n")
global_top <- head(sort(results$omnibus$feature_weights, decreasing = TRUE), 5)
for (i in 1:length(global_top)) {
    cat(sprintf("  %d. %-30s (weight: %.4f)\n", i, names(global_top)[i], global_top[i]))
}

# Access pairwise results
cat("\nSignificant pairwise comparisons:", nrow(results$pairwise$significant_pairs), "\n")
print(results$pairwise$summary_table)
```

## Key Advantages

- **Adaptive**: Learns what matters for your specific data, not generic patterns
- **Robust**: Ensemble approach prevents overfitting and improves generalization  
- **Interpretable**: Feature weights reveal which taxa drive group differences - automatically visualized with VIP charts
- **Directional**: Automatically determines which group has higher abundance for each important taxon
- **Validated**: Permutation testing ensures statistical reliability and proper Type I error control
- **Efficient**: Pre-filtering focuses computation on relevant features
- **User-friendly**: Automatic feature importance visualization with directionality coloring shows exactly which taxa matter most and how they differ

## Options

### Required parameters

- `X`: A matrix of feature abundances with samples as rows and features as columns
- `y`: A vector of group labels for each sample

### Analysis options

- `analysis_type` (default "auto"): Type of analysis - "auto", "pairwise", "omnibus", or "both"
- `n_perms` (default 200): Number of permutations for p-value calculation
- `B` (default 30): Number of weak learners in the ensemble
- `m_frac` (default 0.8): Fraction of features to use in each weak learner
- `show_progress` (default TRUE): Whether to display progress information
- `plot_vip` (default TRUE): Whether to automatically display Variable Importance Plot during analysis
- `correction_method` (default "BH"): Multiple testing correction method for multi-group analysis

### Plotting functions

- `plot_vip(melsi_results, top_n = 15, title = NULL, directionality = TRUE)`: Plot Variable Importance
  - `top_n`: Number of top features to display (default: 15)
  - `title`: Optional custom title
  - `directionality`: Whether to show directionality coloring (default: TRUE)
  
- `plot_pcoa(melsi_results, X, y, title = NULL)`: Plot Principal Coordinates Analysis
  - `X`: Original feature matrix (samples × taxa)
  - `y`: Group labels vector
  - `title`: Optional custom title

- `clr_transform(X, pseudocount = 1)`: Apply CLR transformation
  - `X`: Feature matrix (samples × taxa)
  - `pseudocount`: Small constant to add before log transformation (default: 1)

### Advanced options

- `learning_rate` (default 0.1): Learning rate for gradient descent optimization
- `filter_threshold` (default 0.1): Threshold for pre-filtering feature selection
- `max_iterations` (default 50): Maximum iterations for weak learner optimization

## Troubleshooting

**Question**: When I run MeLSI I see convergence warnings. Is this normal?
**Answer**: Some convergence warnings are normal, especially with small datasets. Check the diagnostics output for condition numbers and eigenvalue ranges.

**Question**: MeLSI is running slowly on my large dataset. How can I speed it up?
**Answer**: Reduce the number of permutations (`n_perms`) or ensemble size (`B`). Pre-filtering is automatically applied to improve efficiency.

**Question**: How do I interpret the learned metric weights?
**Answer**: Higher weights indicate taxa that contribute more to group separation. MeLSI automatically displays the top 5 most important features and generates VIP charts showing feature importance. For multi-group analysis, you get both global feature importance (omnibus) and comparison-specific importance (pairwise). Access all weights via `results$feature_weights` or `results$omnibus$feature_weights` and `results$pairwise$pairwise_results`.

**Question**: Should I use CLR transformation?
**Answer**: Yes, CLR transformation is recommended for microbiome data as it handles compositionality appropriately. MeLSI is designed to work with CLR-transformed data.

**Question**: Why does MeLSI sometimes give higher p-values than traditional methods?
**Answer**: This is actually a feature, not a bug! MeLSI is appropriately conservative and avoids false positives on borderline effects, while still detecting strong signals reliably. The interpretable feature weights provide biological insight even when p-values are similar to traditional methods.

**Question**: When should I choose MeLSI over traditional methods?
**Answer**: Prefer MeLSI for PERMANOVA workflows that need calibrated p-values plus interpretable taxa weights (e.g., diet interventions or subtle host phenotype comparisons). Traditional methods may be preferable when you need the fastest possible computation or when working with very large effect sizes and phylogenetic structure.

## Support

Check out the MeLSI examples in the `vignettes/` folder for an overview of analysis options and example runs. Users should start with the basic usage examples and then refer to the comprehensive testing suite as necessary.

If you have questions, please direct them to the [MeLSI Issues](https://github.com/NathanBresette/MeLSI/issues) page.

## Citation

If you use MeLSI in your research, please cite:

**Preprint:**
```bibtex
@article{bresette2025melsi,
  title={MeLSI: Metric Learning for Statistical Inference in Microbiome Community Composition Analysis},
  author={Bresette, Nathan and Ericsson, Aaron C. and Woods, Carter and Lin, Ai-Ling},
  year={2025},
  journal={bioRxiv},
  doi={10.64898/2025.12.04.692328}
}
```

**Software:**
```bibtex
@software{melsi_software,
  title={MeLSI: Metric Learning for Statistical Inference in Microbiome Analysis},
  author={Bresette, Nathan and Ericsson, Aaron C. and Woods, Carter and Lin, Ai-Ling},
  year={2025},
  url={https://github.com/NathanBresette/MeLSI},
  doi={10.5281/zenodo.17714848}
}
```

## Funding

This work was supported by NIH/NIA R56AG079586.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
