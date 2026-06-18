# Consolidated MeLSI Analysis Function
# Handles both pairwise (2 groups) and multi-group (3+ groups) analysis

# Tag a results list with the "melsi" S3 class without altering its contents,
# so existing `$` access keeps working and a tidy print method is available.
.as_melsi <- function(x) {
    if (!inherits(x, "melsi")) class(x) <- c("melsi", class(x))
    x
}

#' Print method for MeLSI results
#'
#' Concise summary of a \code{\link{melsi}} result: the F-statistic, p-value,
#' and the top features by learned weight.
#'
#' @param x An object of class \code{"melsi"} returned by \code{\link{melsi}}.
#' @param top_n Number of top features to display (default: 5).
#' @param ... Ignored; present for S3 method consistency.
#'
#' @return \code{x}, invisibly.
#'
#' @examples
#' test_data <- generate_test_data(n_samples = 40, n_taxa = 50, n_signal_taxa = 5)
#' X_clr <- clr_transform(test_data$counts)
#' results <- melsi(X_clr, test_data$metadata$Group, n_perms = 19, B = 10,
#'                  show_progress = FALSE, plot_vip = FALSE)
#' print(results)
#'
#' @export
print.melsi <- function(x, top_n = 5, ...) {
    if (!is.null(x$omnibus) || !is.null(x$pairwise)) {
        cat("MeLSI multi-group analysis\n")
        if (!is.null(x$omnibus)) {
            cat(sprintf("  Omnibus: F = %.4f, p = %.4f\n",
                        x$omnibus$F_observed, x$omnibus$p_value))
        }
        if (!is.null(x$pairwise$summary)) {
            cat("  Pairwise comparisons:", nrow(x$pairwise$summary), "\n")
        }
        return(invisible(x))
    }

    cat("MeLSI analysis\n")
    cat(sprintf("  F-statistic: %.4f\n", x$F_observed))
    cat(sprintf("  p-value:     %.4f\n", x$p_value))
    w <- x$feature_weights
    if (!is.null(w) && length(w) > 0) {
        k <- min(top_n, length(w))
        top <- sort(w, decreasing = TRUE)[seq_len(k)]
        cat(sprintf("  Top %d features by weight:\n", k))
        for (i in seq_len(k)) {
            cat(sprintf("    %s (%.4f)\n", names(top)[i], top[i]))
        }
    }
    invisible(x)
}

#' Run MeLSI Analysis
#'
#' Performs MeLSI (Metric Learning for Statistical Inference) analysis for microbiome data.
#' Automatically handles both pairwise comparisons (2 groups) and multi-group analysis (3+ groups).
#'
#' @param X A matrix of feature abundances with samples as rows and features as columns
#' @param y A vector of group labels for each sample
#' @param analysis_type Type of analysis to perform:
#'   - "auto" (default): Automatically choose based on number of groups
#'   - "pairwise": For 2 groups or all pairwise comparisons for 3+ groups
#'   - "omnibus": Global analysis for 3+ groups (requires at least 3 groups)
#'   - "both": Both omnibus and pairwise for 3+ groups
#' @param n_perms Number of permutations for p-value calculation (default: 200)
#' @param B Number of weak learners in the ensemble (default: 30)
#' @param m_frac Fraction of features to use in each weak learner (default: 0.8)
#' @param show_progress Whether to display progress information (default: TRUE)
#' @param plot_vip Whether to display Variable Importance Plot (default: TRUE)
#' @param correction_method Multiple testing correction method for pairwise comparisons (default: "BH")
#' @param BPPARAM A \code{\link[BiocParallel]{BiocParallelParam}} object specifying the
#'   parallel backend to use for permutation testing. If \code{NULL} (default),
#'   permutations run sequentially. Requires the \pkg{BiocParallel} package.
#' @param seed Optional integer used to set the random seed for reproducible
#'   results. If \code{NULL} (default), the current RNG state is used unchanged.
#'
#' @return An object of class \code{"melsi"} (a list, so existing \code{$}
#'         access is unchanged). For 2 groups or pairwise analysis it holds the
#'         F-statistic, p-value, and feature weights; for 3+ groups it contains
#'         the omnibus results, pairwise results, or both.
#'
#' @importFrom stats var dist p.adjust rpois setNames
#' @importFrom utils combn
#' @import ggplot2
#' @importFrom Rcpp evalCpp
#' @useDynLib MeLSI, .registration = TRUE
#'
#' @examples
#' # Generate test data
#' test_data <- generate_test_data(n_samples = 40, n_taxa = 50, n_signal_taxa = 5)
#' X <- test_data$counts
#' y <- test_data$metadata$Group
#'
#' # CLR transformation
#' X_clr <- clr_transform(X)
#'
#' # Run MeLSI analysis
#' results <- melsi(X_clr, y, n_perms = 19, B = 10, show_progress = FALSE)
#'
#' # Check results
#' stopifnot(is.list(results))
#' stopifnot("F_observed" %in% names(results))
#' stopifnot("p_value" %in% names(results))
#'
#' @export
melsi <- function(X, y, analysis_type = "auto", n_perms = 200, B = 30, m_frac = 0.8,
                 show_progress = TRUE, plot_vip = TRUE, correction_method = "BH",
                 BPPARAM = NULL, seed = NULL) {

    # Optional reproducibility: set the RNG seed when supplied (no-op otherwise).
    if (!is.null(seed)) set.seed(seed)

    # Validate that labels match the number of samples (rows of X). The C++
    # F-statistic kernel indexes the label vector by sample, so a length
    # mismatch must be caught here rather than reaching native code.
    if (length(y) != nrow(X)) {
        stop("length(y) (", length(y), ") must equal the number of rows (samples) in X (",
             nrow(X), ").")
    }

    # Validate input and ensure proper column names
    if (is.null(colnames(X)) || all(colnames(X) == "")) {
        colnames(X) <- paste0("Feature_", seq_len(ncol(X)))
        if (show_progress) {
            warning("Input data has no column names. Using generic feature names.")
        }
    }

    groups <- unique(y)
    n_groups <- length(groups)
    
    if (n_groups < 2) {
        stop("MeLSI requires at least 2 groups.")
    }
    
    # Determine analysis type
    if (analysis_type == "auto") {
        if (n_groups == 2) {
            analysis_type <- "pairwise"
        } else {
            analysis_type <- "both"
        }
    }
    
    # Validate analysis type for number of groups
    if (n_groups == 2 && analysis_type %in% c("omnibus", "both")) {
        if (show_progress) {
            message("Only 2 groups detected. Running pairwise analysis.")
        }
        analysis_type <- "pairwise"
    }
    
    if (n_groups >= 3 && analysis_type == "omnibus") {
        # Omnibus only for 3+ groups
        return(.as_melsi(run_omnibus_analysis(X, y, n_perms, B, m_frac, show_progress, plot_vip, BPPARAM)))
    }

    if (analysis_type == "pairwise") {
        # Pairwise analysis
        if (n_groups == 2) {
            # Standard pairwise
            return(.as_melsi(run_pairwise_analysis(X, y, n_perms, B, m_frac, show_progress, plot_vip, BPPARAM)))
        } else {
            # All pairwise comparisons
            return(.as_melsi(run_all_pairwise_analysis(X, y, n_perms, B, m_frac, show_progress,
                                           plot_vip, correction_method, BPPARAM)))
        }
    }

    if (analysis_type == "both") {
        # Both omnibus and pairwise
        if (show_progress) {
            message("=== MeLSI Multi-Group Analysis ===")
            message("Groups: ", paste(groups, collapse = ", "))
            message("Sample sizes: ", paste(names(table(y)), ":", table(y), collapse = ", "))
        }

        results <- list()

        # Run omnibus analysis
        if (show_progress) message("Running omnibus analysis...")
        results$omnibus <- run_omnibus_analysis(X, y, n_perms, B, m_frac, show_progress, plot_vip, BPPARAM)

        # Run pairwise analysis
        if (show_progress) message("Running pairwise analysis...")
        results$pairwise <- run_all_pairwise_analysis(X, y, n_perms, B, m_frac, show_progress,
                                                     plot_vip, correction_method, BPPARAM)

        return(.as_melsi(results))
    }
    
    stop("Invalid analysis_type. Use 'auto', 'pairwise', 'omnibus', or 'both'.")
}

# Helper function: Run standard pairwise analysis (2 groups)
run_pairwise_analysis <- function(X, y, n_perms, B, m_frac, show_progress, plot_vip,
                                  BPPARAM = NULL) {

    if (show_progress) {
        message("--- Starting MeLSI Analysis ---")
    }

    # 1. Learn metric on observed data (with conservative pre-filtering)
    if (show_progress) {
        message("Learning metric on observed data...")
    }

    # Apply conservative pre-filtering
    X_filtered <- apply_conservative_prefiltering(X, y, filter_frac = 0.7)

    M_observed <- learn_melsi_metric_robust(X_filtered, y, B = B, m_frac = m_frac)

    dist_observed <- calculate_mahalanobis_dist_robust(X_filtered, M_observed)
    F_observed <- .melsi_F_scaled(.melsi_scale_mahal(X_filtered, M_observed), y)

    # 2. Generate null distribution with CONSISTENT pre-filtering
    if (show_progress) {
        message("Generating null distribution with ", n_perms, " permutations...")
    }

    .perm_fn_pairwise <- function(p, X, y, B, m_frac) {
        y_permuted <- sample(y)
        X_filtered_perm <- apply_conservative_prefiltering(X, y_permuted, filter_frac = 0.7)
        M_permuted <- learn_melsi_metric_robust(X_filtered_perm, y_permuted, B = B, m_frac = m_frac)
        .melsi_F_scaled(.melsi_scale_mahal(X_filtered_perm, M_permuted), y_permuted)
    }

    use_parallel <- !is.null(BPPARAM) && requireNamespace("BiocParallel", quietly = TRUE)
    if (use_parallel) {
        if (show_progress) message("Using BiocParallel for permutation testing...")
        F_null <- BiocParallel::bplapply(seq_len(n_perms), .perm_fn_pairwise,
                                         X = X, y = y, B = B, m_frac = m_frac,
                                         BPPARAM = BPPARAM)
        F_null <- unlist(F_null)
    } else {
        F_null <- numeric(n_perms)
        pb <- if (show_progress) utils::txtProgressBar(min = 0, max = n_perms, style = 3) else NULL
        for (p in seq_len(n_perms)) {
            F_null[p] <- .perm_fn_pairwise(p, X, y, B, m_frac)
            if (!is.null(pb)) utils::setTxtProgressBar(pb, p)
        }
        if (!is.null(pb)) close(pb)
    }

    # 3. Calculate p-value
    p_value <- (sum(F_null >= F_observed) + 1) / (n_perms + 1)
    
    # 4. Extract feature weights and calculate directionality
    feature_weights <- diag(M_observed)
    names(feature_weights) <- colnames(X_filtered)
    
    # Calculate directionality (which group has higher abundance)
    groups <- unique(y)
    directionality_info <- NULL
    mean_abundances <- NULL
    log2_fold_change <- NULL
    
    # Always calculate directionality for 2 groups
    if (length(groups) == 2) {
        group1_idx <- which(y == groups[1])
        group2_idx <- which(y == groups[2])
        
        # Ensure we have valid indices
        if (length(group1_idx) > 0 && length(group2_idx) > 0) {
        # Calculate mean abundances for each group
        mean_group1 <- colMeans(X_filtered[group1_idx, , drop = FALSE])
        mean_group2 <- colMeans(X_filtered[group2_idx, , drop = FALSE])
        
            # Determine directionality - ensure it's always a character vector with names
        directionality_info <- ifelse(mean_group1 > mean_group2, 
                                          paste0("Higher in ", as.character(groups[1])), 
                                          paste0("Higher in ", as.character(groups[2])))
        names(directionality_info) <- colnames(X_filtered)
            
            # Verify directionality was created correctly
            if (is.null(directionality_info) || length(directionality_info) == 0) {
                warning("Failed to calculate directionality. Creating default values.")
                directionality_info <- rep("Unknown", ncol(X_filtered))
                names(directionality_info) <- colnames(X_filtered)
            }
        
        # Calculate fold change and log2 fold change
        # Add small epsilon to both to avoid division issues and ensure positive values for log2
        fold_change <- (mean_group1 + 1e-10) / (mean_group2 + 1e-10)
        log2_fold_change <- suppressWarnings(log2(fold_change))
        log2_fold_change[!is.finite(log2_fold_change)] <- 0
        names(log2_fold_change) <- colnames(X_filtered)
        
        # Store mean abundances
        mean_abundances <- list(
            group1 = mean_group1,
            group2 = mean_group2,
            group1_name = as.character(groups[1]),
            group2_name = as.character(groups[2])
        )
        } else {
            warning("Invalid group indices. Cannot calculate directionality.")
        }
    }
    
    if (show_progress) {
        message("Analysis completed!")
        message("F-statistic: ", round(F_observed, 4))
        message("P-value: ", round(p_value, 4))
        message("\nFeature importance: Access results$feature_weights to see which taxa")
        message("contributed most to the learned distance metric.")
        if (!is.null(directionality_info)) {
            message("Directionality: Access results$directionality to see which group has higher abundance.")
        }
        
        # Show top 5 features with directionality if available
        if (length(feature_weights) >= 5) {
            top_5_idx <- order(feature_weights, decreasing = TRUE)[seq_len(5)]
            message("\nTop 5 most important features:")
            for (i in seq_len(5)) {
                idx <- top_5_idx[i]
                feature_name <- names(feature_weights)[idx]
                if (is.null(feature_name) || feature_name == "" || is.na(feature_name)) {
                    feature_name <- paste0("Feature_", idx)
                }
                direction_text <- ""
                if (!is.null(directionality_info) && !is.null(directionality_info[idx])) {
                    direction_text <- paste0(" [", directionality_info[idx], "]")
                }
                message(sprintf("  %d. %s (weight: %.4f)%s", i, feature_name, feature_weights[idx], direction_text))
            }
        }
    }
    
    # 5. Create VIP plot if requested
    if (plot_vip && length(feature_weights) > 0) {
        tryCatch({
            plot_feature_importance(feature_weights, directionality = directionality_info)
        }, error = function(e) {
            if (show_progress) {
                warning("Could not generate VIP plot: ", e$message)
            }
        })
    }
    
    # Return results - directionality should always be included for 2-group analysis
    # (will be NULL for multi-group, but should be a named vector for 2 groups)
    return(list(
        F_observed = F_observed,
        p_value = p_value,
        F_null = F_null,
        feature_weights = feature_weights,
        directionality = directionality_info,  # Named vector: "Higher in [group]" for each feature
        mean_abundances = mean_abundances,
        log2_fold_change = log2_fold_change,
        metric_matrix = M_observed,
        distance_matrix = dist_observed,  # reuse already-computed dist
        diagnostics = list(
            n_features_used = ncol(X_filtered),
            n_permutations = n_perms,
            ensemble_size = B
        )
    ))
}

# Helper function: Run omnibus analysis (3+ groups)
run_omnibus_analysis <- function(X, y, n_perms, B, m_frac, show_progress, plot_vip,
                                  BPPARAM = NULL) {
    
    groups <- unique(y)
    n_groups <- length(groups)
    
    if (show_progress) {
        message("--- Starting MeLSI Omnibus Analysis ---")
        message("Groups: ", paste(groups, collapse = ", "))
        message("Sample sizes: ", paste(table(y), collapse = " "))
    }
    
    # 1. Apply conservative pre-filtering
    if (show_progress) {
        message("Applying conservative pre-filtering...")
    }
    
    X_filtered <- apply_conservative_prefiltering_multi(X, y, filter_frac = 0.7)
    
    # 2. Learn global metric optimized for all group pairs
    if (show_progress) {
        message("Learning global metric for all group pairs...")
    }
    
    M_observed <- learn_melsi_metric_omnibus(X_filtered, y, B = B, m_frac = m_frac)
    
    # 3. Calculate omnibus F-statistic
    dist_observed <- calculate_mahalanobis_dist_robust(X_filtered, M_observed)
    F_observed <- .melsi_F_scaled(.melsi_scale_mahal(X_filtered, M_observed), y)
    
    # 4. Generate null distribution
    if (show_progress) {
        message("Generating null distribution with ", n_perms, " permutations...")
    }
    
    .perm_fn_omnibus <- function(p, X, y, B, m_frac) {
        y_permuted <- sample(y)
        X_filtered_perm <- apply_conservative_prefiltering_multi(X, y_permuted, filter_frac = 0.7)
        M_permuted <- learn_melsi_metric_omnibus(X_filtered_perm, y_permuted, B = B, m_frac = m_frac)
        .melsi_F_scaled(.melsi_scale_mahal(X_filtered_perm, M_permuted), y_permuted)
    }

    use_parallel <- !is.null(BPPARAM) && requireNamespace("BiocParallel", quietly = TRUE)
    if (use_parallel) {
        if (show_progress) message("Using BiocParallel for permutation testing...")
        F_null <- BiocParallel::bplapply(seq_len(n_perms), .perm_fn_omnibus,
                                         X = X, y = y, B = B, m_frac = m_frac,
                                         BPPARAM = BPPARAM)
        F_null <- unlist(F_null)
    } else {
        F_null <- numeric(n_perms)
        pb <- if (show_progress) utils::txtProgressBar(min = 0, max = n_perms, style = 3) else NULL
        for (p in seq_len(n_perms)) {
            F_null[p] <- .perm_fn_omnibus(p, X, y, B, m_frac)
            if (!is.null(pb)) utils::setTxtProgressBar(pb, p)
        }
        if (!is.null(pb)) close(pb)
    }

    # 5. Calculate p-value
    p_value <- (sum(F_null >= F_observed) + 1) / (n_perms + 1)
    
    # 6. Extract feature weights and calculate directionality (highest mean group)
    feature_weights <- diag(M_observed)
    names(feature_weights) <- colnames(X_filtered)
    
    # Calculate which group has highest mean abundance for each feature
    directionality_info <- NULL
    mean_abundances <- NULL
    
    if (n_groups >= 2) {
        # Calculate mean abundances for each group
        mean_by_group <- list()
        for (g in groups) {
            group_idx <- which(y == g)
            mean_by_group[[as.character(g)]] <- colMeans(X_filtered[group_idx, , drop = FALSE])
        }
        
        # Determine which group has highest mean for each feature
        mean_matrix <- do.call(rbind, mean_by_group)
        max_group_idx <- apply(mean_matrix, 2, which.max)
        directionality_info <- as.character(groups[max_group_idx])
        names(directionality_info) <- colnames(X_filtered)
        
        # Store mean abundances
        mean_abundances <- mean_by_group
        names(mean_abundances) <- as.character(groups)
    }
    
    if (show_progress) {
        message("Omnibus analysis completed!")
        message("F-statistic: ", round(F_observed, 4))
        message("P-value: ", round(p_value, 4))
        message("\nGlobal feature importance: Access results$feature_weights to see which taxa")
        message("contributed most to overall group separation.")
        if (!is.null(directionality_info)) {
            message("Directionality: Access results$directionality to see which group has highest mean abundance.")
        }
        
        # Show top 5 features with directionality if available
        if (length(feature_weights) >= 5) {
            top_5_idx <- order(feature_weights, decreasing = TRUE)[seq_len(5)]
            message("\nTop 5 globally important features:")
            for (i in seq_len(5)) {
                idx <- top_5_idx[i]
                feature_name <- names(feature_weights)[idx]
                if (is.null(feature_name) || feature_name == "" || is.na(feature_name)) {
                    feature_name <- paste0("Feature_", idx)
                }
                direction_text <- ""
                if (!is.null(directionality_info) && !is.null(directionality_info[idx])) {
                    direction_text <- paste0(" [Highest in ", directionality_info[idx], "]")
                }
                message(sprintf("  %d. %s (weight: %.4f)%s", i, feature_name, feature_weights[idx], direction_text))
            }
        }
    }
    
    # 7. Create VIP plot if requested
    if (plot_vip && length(feature_weights) > 0) {
        tryCatch({
            plot_feature_importance(feature_weights, 
                                  main_title = "Global Feature Importance (Omnibus)",
                                  directionality = directionality_info)
        }, error = function(e) {
            if (show_progress) {
                warning("Could not generate VIP plot: ", e$message)
            }
        })
    }
    
    return(list(
        F_observed = F_observed,
        p_value = p_value,
        F_null = F_null,
        feature_weights = feature_weights,
        directionality = directionality_info,
        mean_abundances = mean_abundances,
        metric_matrix = M_observed,
        distance_matrix = dist_observed,  # reuse already-computed dist
        group_info = list(
            groups = groups,
            n_groups = n_groups,
            sample_sizes = table(y)
        ),
        diagnostics = list(
            n_features_used = ncol(X_filtered),
            n_permutations = n_perms,
            ensemble_size = B
        )
    ))
}

# Helper function: Run all pairwise comparisons (3+ groups)
run_all_pairwise_analysis <- function(X, y, n_perms, B, m_frac, show_progress, plot_vip,
                                      correction_method, BPPARAM = NULL) {
    
    groups <- unique(y)
    n_groups <- length(groups)
    
    if (show_progress) {
        message("--- Starting MeLSI Pairwise Analysis ---")
        message("Groups: ", paste(groups, collapse = ", "))
        message("Total comparisons: ", choose(n_groups, 2))
    }
    
    # Generate all pairwise combinations
    pairwise_combinations <- combn(groups, 2, simplify = FALSE)
    n_comparisons <- length(pairwise_combinations)
    
    # Store results
    pairwise_results <- list()
    summary_data <- data.frame(
        Group1 = character(n_comparisons),
        Group2 = character(n_comparisons),
        F_statistic = numeric(n_comparisons),
        P_value = numeric(n_comparisons),
        stringsAsFactors = FALSE
    )
    
    # Run pairwise comparisons
    for (i in seq_len(n_comparisons)) {
        pair <- pairwise_combinations[[i]]
        group1 <- pair[1]
        group2 <- pair[2]
        
        if (show_progress) {
            message(sprintf("Comparison %d/%d: %s vs %s", i, n_comparisons, group1, group2))
        }
        
        # Subset data for this pair
        pair_indices <- y %in% pair
        X_pair <- X[pair_indices, , drop = FALSE]
        y_pair <- y[pair_indices]
        
        # Run MeLSI for this pair
        pair_result <- tryCatch({
            run_pairwise_analysis(X_pair, y_pair, n_perms = n_perms, B = B, m_frac = m_frac,
                                 show_progress = show_progress, plot_vip = FALSE,
                                 BPPARAM = BPPARAM)
        }, error = function(e) {
            if (show_progress) {
                warning("Analysis failed for ", group1, " vs ", group2, ": ", e$message)
            }
            return(NULL)
        })
        
        if (!is.null(pair_result)) {
            # Store results
            comparison_name <- paste(group1, group2, sep = "_vs_")
            pairwise_results[[comparison_name]] <- pair_result
            
            # Update summary table
            summary_data$Group1[i] <- group1
            summary_data$Group2[i] <- group2
            summary_data$F_statistic[i] <- pair_result$F_observed
            summary_data$P_value[i] <- pair_result$p_value
            
            if (show_progress) {
                message(sprintf("  F-statistic: %.4f, P-value: %.4f", 
                           pair_result$F_observed, pair_result$p_value))
            }
        }
    }
    
    # Apply multiple testing correction
    summary_data$P_value_corrected <- p.adjust(summary_data$P_value, method = correction_method)
    summary_data$Significant <- summary_data$P_value_corrected < 0.05
    
    # Identify significant pairs
    significant_pairs <- summary_data[summary_data$Significant, ]
    
    if (show_progress) {
        message("\n=== Pairwise Analysis Summary ===")
        message("Multiple testing correction: ", correction_method)
        message("Significant pairs (corrected p < 0.05): ", nrow(significant_pairs))
        
        if (nrow(significant_pairs) > 0) {
            message("\nSignificant comparisons:")
            for (i in seq_len(nrow(significant_pairs))) {
                row <- significant_pairs[i, ]
                message(sprintf("  %s vs %s: F=%.3f, p=%.4f (corrected: %.4f)",
                           row$Group1, row$Group2, row$F_statistic, 
                           row$P_value, row$P_value_corrected))
            }
        }
    }
    
    return(list(
        pairwise_results = pairwise_results,
        summary_table = summary_data,
        significant_pairs = significant_pairs,
        correction_method = correction_method
    ))
}

# All helper functions from original melsi_robust.R
# (Conservative pre-filtering, metric learning, distance calculation, etc.)

# Conservative pre-filtering function
apply_conservative_prefiltering <- function(X, y, filter_frac = 0.7) {
    # Keep more features, use less aggressive filtering
    classes <- unique(y)
    if (length(classes) != 2 || ncol(X) <= 10) {
        return(X)
    }
    
    class1_indices <- which(y == classes[1])
    class2_indices <- which(y == classes[2])
    
    # Calculate feature importance with more conservative approach
    mean_group1    <- colMeans(X[class1_indices, , drop = FALSE])
    mean_group2    <- colMeans(X[class2_indices, , drop = FALSE])
    var_group1     <- apply(X[class1_indices, , drop = FALSE], 2, var)
    var_group2     <- apply(X[class2_indices, , drop = FALSE], 2, var)
    var_combined   <- sqrt(pmax(var_group1 + var_group2, 0))
    feature_importance <- abs(mean_group1 - mean_group2) / (var_combined + 1e-10)
    feature_importance[!is.finite(feature_importance)] <- 0
    
    # Keep more features (70% instead of 50%)
    n_keep <- max(10, floor(ncol(X) * filter_frac))
    top_features <- order(feature_importance, decreasing = TRUE)[seq_len(n_keep)]
    
    # Ensure column names are preserved
    filtered_X <- X[, top_features, drop = FALSE]
    colnames(filtered_X) <- colnames(X)[top_features]
    
    return(filtered_X)
}

# Helper function: Conservative pre-filtering for multi-group data
apply_conservative_prefiltering_multi <- function(X, y, filter_frac = 0.7) {
    classes <- unique(y)
    if (length(classes) < 2 || ncol(X) <= 10) {
        return(X)
    }
    
    # Calculate feature importance using vectorized one-way ANOVA F-statistic
    grand_mean <- colMeans(X)
    n  <- nrow(X)
    k  <- length(classes)
    SS_between <- numeric(ncol(X))
    SS_within  <- numeric(ncol(X))
    for (g in classes) {
        idx <- which(y == g)
        n_g <- length(idx)
        gm  <- colMeans(X[idx, , drop = FALSE])
        SS_between <- SS_between + n_g * (gm - grand_mean)^2
        X_c <- sweep(X[idx, , drop = FALSE], 2, gm, "-")
        SS_within  <- SS_within + colSums(X_c^2)
    }
    feature_importance <- (SS_between / (k - 1)) / (SS_within / (n - k) + 1e-10)
    feature_importance[!is.finite(feature_importance)] <- 0
    
    # Keep top features
    n_keep <- max(10, floor(ncol(X) * filter_frac))
    top_features <- order(feature_importance, decreasing = TRUE)[seq_len(n_keep)]
    
    # Ensure column names are preserved
    filtered_X <- X[, top_features, drop = FALSE]
    colnames(filtered_X) <- colnames(X)[top_features]
    
    return(filtered_X)
}

# Helper function: PERMANOVA F-statistic from an already-scaled data matrix,
# via the fused C++ kernel. Avoids building and squaring the full n x n distance
# matrix; used in the hot paths where only the F-statistic (not the distance
# matrix) is needed.
.melsi_F_scaled <- function(Xs, labels) {
    # Guard the native boundary: the kernel indexes labels by sample.
    if (length(labels) != nrow(Xs)) {
        stop("Internal error: label length does not match number of samples.")
    }
    g <- match(labels, unique(labels)) - 1L
    melsi_permanova_f(t(Xs), as.integer(g), length(unique(labels)))
}

# Apply the Mahalanobis column scaling (matches calculate_mahalanobis_dist_robust)
# without forming a distance object.
.melsi_scale_mahal <- function(X, M) {
    w <- 1 / sqrt(pmax(diag(M), 1e-6))
    w[!is.finite(w)] <- 1e-3
    sweep(X, 2, w, "*")
}

# Helper function: Calculate PERMANOVA F-statistic (direct formula, avoids adonis2 overhead)
calculate_permanova_F <- function(dist_obj, labels) {
    n  <- attr(dist_obj, "Size")
    groups <- unique(labels)
    k  <- length(groups)
    D2 <- as.matrix(dist_obj)^2
    SS_total  <- sum(D2) / (2 * n)
    SS_within <- 0
    for (g in groups) {
        idx <- which(labels == g)
        SS_within <- SS_within + sum(D2[idx, idx]) / (2 * length(idx))
    }
    SS_between <- SS_total - SS_within
    (SS_between / (k - 1)) / (SS_within / (n - k))
}

# Helper function: Robust Mahalanobis Distance
# M is always diagonal in this implementation (only diag(M) is ever modified),
# so we skip the O(p^3) eigen decomposition and scale columns directly.
calculate_mahalanobis_dist_robust <- function(X, M) {
    w <- 1 / sqrt(pmax(diag(M), 1e-6))
    w[!is.finite(w)] <- 1e-3
    return(dist(sweep(X, 2, w, "*"), method = "euclidean"))
}

# Helper function: Optimize weak learner
optimize_weak_learner_robust <- function(X, y, n_iterations = 50, learning_rate = 0.1) {
    n_features <- ncol(X)

    # Diagonal-metric gradient descent in C++ (melsi_opt_weak_learner). The C++
    # loop draws its within/between pairs with R's own index sampler, so it
    # consumes the global RNG identically to the previous pure-R loop and returns
    # a bit-identical diagonal. Group ids match the R reference: g == 0 is
    # unique(y)[1] (class 1), g == 1 is unique(y)[2] (class 2); the early-stopping
    # F statistic uses all groups. Guards (k < 2, class size < 2) are handled in
    # C++ and return the identity diagonal, matching the old return(diag(...)).
    g <- match(y, unique(y)) - 1L
    k <- length(unique(y))
    diag_M <- melsi_opt_weak_learner(X, as.integer(g), k, n_iterations, learning_rate)
    diag(pmax(diag_M, 0.01))
}

# Helper function: Ensemble metric learning with bootstrap and feature subsampling
# Shared logic for both robust (pairwise) and omnibus (multi-group) metric learning
.learn_ensemble_metric <- function(X, y, B, m_frac, optimizer_fn) {
    n_samples <- nrow(X)
    n_features <- ncol(X)
    m <- max(2, floor(n_features * m_frac))

    learned_matrices <- vector("list", B)
    valid_count <- 0
    f_stats <- numeric(B)

    for (b in seq_len(B)) {
        boot_indices <- sample(seq_len(n_samples), n_samples, replace = TRUE)
        if (length(unique(y[boot_indices])) < 2) next

        X_boot <- X[boot_indices, , drop = FALSE]
        y_boot <- y[boot_indices]

        feature_indices <- sample(seq_len(n_features), m, replace = FALSE)
        X_subset <- X_boot[, feature_indices, drop = FALSE]

        M_weak <- optimizer_fn(X_subset, y_boot)

        tryCatch({
            f_stat <- .melsi_F_scaled(.melsi_scale_mahal(X_subset, M_weak), y_boot)

            if (is.finite(f_stat) && f_stat > 0) {
                valid_count <- valid_count + 1
                # Store only the diagonal (M is always diagonal); identity (1) elsewhere
                diag_full <- rep(1, n_features)
                diag_full[feature_indices] <- diag(M_weak)
                learned_matrices[[valid_count]] <- diag_full
                f_stats[valid_count] <- f_stat
            }
        }, error = function(e) {
            # Skip this weak learner if it fails
        })

        if (valid_count >= B) break
    }

    if (valid_count == 0) {
        warning("No valid weak learners found. Returning identity matrix.")
        return(diag(n_features))
    }

    weights <- f_stats[seq_len(valid_count)]
    weights <- weights / sum(weights)

    # learned_matrices holds diagonal vectors; weighted sum directly
    diag_ensemble <- numeric(n_features)
    for (i in seq_len(valid_count)) {
        diag_ensemble <- diag_ensemble + weights[i] * learned_matrices[[i]]
    }
    return(diag(pmax(diag_ensemble, 1e-6)))
}

# Helper function: Learn MeLSI metric
# Pre-filtering is handled upstream by apply_conservative_prefiltering(), so the
# metric learner operates directly on the already-filtered feature matrix.
learn_melsi_metric_robust <- function(X, y, B = 20, m_frac = 0.7) {
    .learn_ensemble_metric(X, y, B, m_frac, optimize_weak_learner_robust)
}

# Helper function: Learn omnibus metric for multi-group data
learn_melsi_metric_omnibus <- function(X, y, B = 30, m_frac = 0.8) {
    .learn_ensemble_metric(X, y, B, m_frac, optimize_weak_learner_omnibus)
}

# Helper function: Optimize weak learner for omnibus analysis
optimize_weak_learner_omnibus <- function(X, y, n_iterations = 50, learning_rate = 0.1) {
    n_samples <- nrow(X)
    n_features <- ncol(X)
    
    # Start with identity matrix
    M <- diag(n_features)
    
    # Get all group pairs
    groups <- unique(y)
    n_groups <- length(groups)
    
    if (n_groups < 2) return(M)
    
    # Track convergence
    prev_f_stat <- -Inf
    stagnation_count <- 0
    group_pairs <- combn(groups, 2, simplify = FALSE)

    for (iter in seq_len(n_iterations)) {
        # Sample from all group pairs (balanced sampling)
        
        # Randomly select a group pair to optimize for this iteration
        selected_pair <- sample(group_pairs, 1)[[1]]
        group1_indices <- which(y == selected_pair[1])
        group2_indices <- which(y == selected_pair[2])
        
        if (length(group1_indices) < 2 || length(group2_indices) < 2) next
        
        # Sample from selected pair (same as pairwise optimization)
        i1 <- sample(group1_indices, 1)
        j1 <- group1_indices[group1_indices != i1][[sample.int(length(group1_indices) - 1L, 1)]]
        i2 <- sample(group2_indices, 1)
        j2 <- group2_indices[group2_indices != i2][[sample.int(length(group2_indices) - 1L, 1)]]
        
        # Compute differences
        diff1 <- X[i1, ] - X[j1, ]  # Within group 1
        diff2 <- X[i2, ] - X[j2, ]  # Within group 2
        diff3 <- X[i1, ] - X[i2, ]  # Between groups
        
        # Gradient calculation
        grad_between <- diff3^2
        grad_within <- -(diff1^2 + diff2^2) / 2
        total_gradient <- grad_between + grad_within
        
        # Adaptive learning rate
        current_learning_rate <- learning_rate * (1 / (1 + iter * 0.1))
        
        # Update
        diag(M) <- diag(M) + current_learning_rate * total_gradient
        diag(M) <- pmax(diag(M), 0.01)
        
        # Early stopping check
        if (iter %% 20 == 0) {
            current_f_stat <- tryCatch({
                .melsi_F_scaled(sweep(X, 2, sqrt(pmax(diag(M), 0)), "*"), y)
            }, error = function(e) 0)
            
            if (current_f_stat <= prev_f_stat) {
                stagnation_count <- stagnation_count + 1
                if (stagnation_count >= 5) break
            } else {
                stagnation_count <- 0
            }
            prev_f_stat <- current_f_stat
        }
    }
    
    return(M)
}

#' Plot Feature Importance from MeLSI Analysis
#'
#' Creates a barplot showing the top features ranked by their learned weights
#'
#' @param feature_weights Named vector of feature weights
#' @param top_n Number of top features to display (default: 8)
#' @param main_title Optional title for the plot
#' @param directionality Optional named vector indicating which group has higher abundance for each feature
#'
#' @return A ggplot2 object (invisibly)
#'
#' @examples
#' # Generate test data and run MeLSI
#' test_data <- generate_test_data(n_samples = 30, n_taxa = 20, n_signal_taxa = 5)
#' X <- test_data$counts
#' y <- test_data$metadata$Group
#' X_clr <- clr_transform(X)
#' results <- melsi(X_clr, y, n_perms = 19, B = 10, show_progress = FALSE)
#' 
#' # Plot feature importance
#' plot_feature_importance(results$feature_weights, top_n = 10)
#'
#' @export
plot_feature_importance <- function(feature_weights, top_n = 8, main_title = NULL, directionality = NULL) {
    # Validate input
    if (length(feature_weights) == 0) {
        stop("feature_weights is empty")
    }
    
    # Sort features by weight
    sorted_weights <- sort(feature_weights, decreasing = TRUE)
    
    # Take top N features
    n_display <- min(top_n, length(sorted_weights))
    top_weights <- sorted_weights[seq_len(n_display)]
    
    # Get feature names
    feature_names <- names(top_weights)
    if (is.null(feature_names)) {
        feature_names <- paste0("Feature_", seq_len(n_display))
    }
    
    # Truncate long names for better display (more aggressive truncation)
    feature_names <- ifelse(nchar(feature_names) > 20, 
                           paste0(substr(feature_names, 1, 17), "..."),
                           feature_names)
    
    # Ensure unique feature names (in case truncation created duplicates)
    # Add suffix to duplicates
    if (any(duplicated(feature_names))) {
        for (i in seq_along(feature_names)) {
            if (sum(feature_names == feature_names[i]) > 1) {
                # Find all occurrences of this name
                dup_indices <- which(feature_names == feature_names[i])
                # Keep first as-is, add numbers to others
                if (length(dup_indices) > 1) {
                    for (j in seq_len(length(dup_indices))[-1]) {
                        feature_names[dup_indices[j]] <- paste0(feature_names[dup_indices[j]], "_", j)
                    }
                }
            }
        }
    }
    
    # Extract directionality for top features if provided
    directionality_colors <- NULL
    directionality_labels <- NULL
    if (!is.null(directionality) && length(directionality) > 0) {
        # Match directionality to top features
        directionality_labels <- directionality[names(top_weights)]
        
        # Create color mapping based on directionality
        # Extract unique group names from directionality strings
        unique_groups <- unique(directionality_labels)
        if (length(unique_groups) == 2) {
            # Two groups - use two colors
            directionality_colors <- ifelse(directionality_labels == unique_groups[1], 
                                           "#E63946",  # Red for group 1
                                           "#457B9D")  # Blue for group 2
        } else if (length(unique_groups) >= 3) {
            # Three or more groups - assign distinct colors
            # Use a color palette that works well for 3+ groups
            group_colors <- c("#E63946", "#457B9D", "#F77F00", "#06A77D", "#7209B7", "#A8DADC")
            # Cycle through colors if more than 6 groups
            color_map <- stats::setNames(
                group_colors[seq_len(min(length(unique_groups), length(group_colors)))],
                unique_groups[seq_len(min(length(unique_groups), length(group_colors)))]
            )
            # If more groups than colors, cycle
            if (length(unique_groups) > length(group_colors)) {
                remaining <- unique_groups[(length(group_colors) + 1):length(unique_groups)]
                extra_colors <- rep(group_colors, length.out = length(remaining))
                color_map <- c(color_map, stats::setNames(extra_colors, remaining))
            }
            directionality_colors <- color_map[directionality_labels]
        } else {
            # Single group or missing - use default color
            directionality_colors <- rep("steelblue", length(directionality_labels))
        }
    }
    
    # Create data frame for ggplot
    plot_data <- data.frame(
        Feature = factor(feature_names, levels = rev(feature_names)),  # Reverse for top-to-bottom ordering
        Weight = top_weights
    )
    
    # Add directionality info if available
    if (!is.null(directionality_labels)) {
        plot_data$Directionality <- directionality_labels
        plot_data$Color <- directionality_colors
    }
    
    # Create title
    title_text <- if (!is.null(main_title)) main_title else paste0("Top ", n_display, " Features by Importance")
    
    # Create ggplot with or without directionality coloring
    if (!is.null(directionality_labels) && !is.null(directionality_colors)) {
        # Create proper color mapping for scale_fill_manual
        # Map each unique group to its corresponding color
        unique_groups <- unique(directionality_labels)
        # Get the color for the first occurrence of each unique group
        color_map <- stats::setNames(
            vapply(unique_groups, function(g) {
                idx <- which(directionality_labels == g)[1]
                directionality_colors[idx]
            }, character(1)),
            unique_groups
        )
        
        # Plot with directionality colors
        p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Weight, y = Feature, fill = Directionality)) +
            ggplot2::geom_col() +
            ggplot2::scale_fill_manual(values = color_map,
                                      name = "Higher in") +
            ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", Weight)), 
                              hjust = -0.1, size = 3) +
            ggplot2::labs(
                title = title_text,
                x = "Feature Weight",
                y = ""
            ) +
            ggplot2::theme_bw() +
            ggplot2::theme(
                plot.title = ggplot2::element_text(size = 12, hjust = 0.5),
                axis.text.y = ggplot2::element_text(size = 10),
                axis.text.x = ggplot2::element_text(size = 9),
                axis.title.x = ggplot2::element_text(size = 11),
                panel.grid.minor = ggplot2::element_blank(),
                plot.margin = ggplot2::margin(10, 30, 10, 10),
                legend.position = "right"
            ) +
            ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.1)))
    } else {
        # Plot without directionality (original behavior)
        p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Weight, y = Feature)) +
            ggplot2::geom_col(fill = "steelblue") +
            ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", Weight)), 
                              hjust = -0.1, size = 3) +
            ggplot2::labs(
                title = title_text,
                x = "Feature Weight",
                y = ""
            ) +
            ggplot2::theme_bw() +
            ggplot2::theme(
                plot.title = ggplot2::element_text(size = 12, hjust = 0.5),
                axis.text.y = ggplot2::element_text(size = 10),
                axis.text.x = ggplot2::element_text(size = 9),
                axis.title.x = ggplot2::element_text(size = 11),
                panel.grid.minor = ggplot2::element_blank(),
                plot.margin = ggplot2::margin(10, 30, 10, 10)
            ) +
            ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.1)))
    }
    
    # Return the plot object for further customization if needed    
    return(invisible(p))
}

#' Plot VIP from MeLSI Results (User-Friendly Wrapper)
#'
#' Simplified function to plot Variable Importance (VIP) directly from MeLSI results.
#' Automatically extracts feature weights and optionally includes directionality information.
#'
#' @param melsi_results Results object from melsi() function
#' @param top_n Number of top features to display (default: 15)
#' @param title Optional custom title for the plot
#' @param directionality Whether to include directionality coloring (default: TRUE)
#'
#' @return A ggplot2 object (invisibly)
#'
#' @examples
#' # Generate test data and run MeLSI
#' test_data <- generate_test_data(n_samples = 30, n_taxa = 20, n_signal_taxa = 5)
#' X <- test_data$counts
#' y <- test_data$metadata$Group
#' X_clr <- clr_transform(X)
#' results <- melsi(X_clr, y, n_perms = 19, B = 10, show_progress = FALSE)
#' 
#' # Plot VIP with directionality (default)
#' plot_vip(results, top_n = 10)
#'
#' @export
plot_vip <- function(melsi_results, top_n = 15, title = NULL, directionality = TRUE) {
    
    # Check if results is valid
    if (is.null(melsi_results) || !is.list(melsi_results)) {
        stop("melsi_results must be a valid MeLSI results object")
    }
    
    # Extract feature weights
    if (is.null(melsi_results$feature_weights)) {
        stop("No feature weights found in results. Make sure you ran melsi() successfully.")
    }
    
    feature_weights <- melsi_results$feature_weights
    
    # Extract directionality if requested
    directionality_data <- NULL
    if (directionality) {
        directionality_data <- melsi_results$directionality
    }
    
    # Call the underlying plot function
    plot_feature_importance(
        feature_weights = feature_weights,
        top_n = top_n,
        main_title = title,
        directionality = directionality_data
    )
}

#' Plot PCoA from MeLSI Results
#'
#' Creates a Principal Coordinates Analysis (PCoA) plot using the learned MeLSI distance matrix.
#'
#' @param melsi_results Results object from melsi() function
#' @param X Original feature matrix (samples x taxa)
#' @param y Group labels vector
#' @param title Optional custom title for the plot (default: "PCoA using MeLSI Distance")
#'
#' @return A ggplot2 object (invisibly)
#'
#' @importFrom stats cmdscale
#'
#' @examples
#' # Generate test data and run MeLSI
#' test_data <- generate_test_data(n_samples = 30, n_taxa = 20, n_signal_taxa = 5)
#' X <- test_data$counts
#' y <- test_data$metadata$Group
#' X_clr <- clr_transform(X)
#' results <- melsi(X_clr, y, n_perms = 19, B = 10, show_progress = FALSE)
#' 
#' # Plot PCoA
#' plot_pcoa(results, X_clr, y)
#'
#' @export
plot_pcoa <- function(melsi_results, X, y, title = "PCoA using MeLSI Distance") {
    
    # Check if results has distance matrix
    if (is.null(melsi_results$distance_matrix)) {
        stop("No distance matrix found in results.")
    }
    
    dist_matrix <- melsi_results$distance_matrix
    
    # Run PCoA (classical MDS)
    pcoa_result <- cmdscale(dist_matrix, k = 2, eig = TRUE)
    
    # Calculate variance explained
    var_explained <- pcoa_result$eig / sum(abs(pcoa_result$eig)) * 100
    
    # Create plot data
    plot_data <- data.frame(
        PC1 = pcoa_result$points[, 1],
        PC2 = pcoa_result$points[, 2],
        Group = as.factor(y)
    )
    
    # Create ggplot
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = PC1, y = PC2, color = Group)) +
        ggplot2::geom_point(size = 3, alpha = 0.7) +
        ggplot2::stat_ellipse(level = 0.95, linetype = 2) +
        ggplot2::labs(
            title = title,
            x = sprintf("PCoA1 (%.1f%%)", var_explained[1]),
            y = sprintf("PCoA2 (%.1f%%)", var_explained[2])
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
            legend.position = "right"
        )
    
    return(invisible(p))
}

#' CLR Transformation for Microbiome Data
#'
#' Applies centered log-ratio (CLR) transformation to microbiome count data.
#' This transformation is recommended for microbiome data before running MeLSI.
#'
#' @param X Feature matrix (samples x taxa) with raw counts or relative abundances
#' @param pseudocount Small constant to add before log transformation (default: 1)
#'
#' @return CLR-transformed matrix with preserved column names
#'
#' @examples
#' # Generate synthetic data
#' test_data <- generate_test_data(n_samples = 20, n_taxa = 30, n_signal_taxa = 5)
#' X <- test_data$counts
#' 
#' # Transform microbiome data
#' X_clr <- clr_transform(X)
#' 
#' # Verify transformation
#' stopifnot(is.matrix(X_clr))
#' stopifnot(nrow(X_clr) == nrow(X))
#'
#' @export
clr_transform <- function(X, pseudocount = 1) {
    
    # Add pseudocount and take log
    X_log <- log(X + pseudocount)
    
    # Center by row means (CLR transformation)
    X_clr <- X_log - rowMeans(X_log)
    
    # Preserve column names
    colnames(X_clr) <- colnames(X)
    rownames(X_clr) <- rownames(X)
    
    return(X_clr)
}
