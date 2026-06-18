# Test file for MeLSI package
# Comprehensive tests for all MeLSI functionality

library(testthat)
library(MeLSI)

test_that("generate_test_data works correctly", {
  # Test basic generation
  test_data <- generate_test_data(n_samples = 20, n_taxa = 30, n_signal_taxa = 5)
  
  expect_true(is.list(test_data))
  expect_true("counts" %in% names(test_data))
  expect_true("metadata" %in% names(test_data))
  expect_true("signal_taxa" %in% names(test_data))
  
  expect_equal(nrow(test_data$counts), 20)
  expect_equal(ncol(test_data$counts), 30)
  expect_equal(nrow(test_data$metadata), 20)
  expect_equal(length(test_data$signal_taxa), 5)
  
  # Test counts are non-negative
  expect_true(all(test_data$counts >= 0))
  
  # Test groups are balanced
  expect_true("Group" %in% names(test_data$metadata))
})

test_that("clr_transform works correctly", {
  # Generate test data
  test_data <- generate_test_data(n_samples = 20, n_taxa = 30, n_signal_taxa = 5)
  X <- test_data$counts
  
  # Test CLR transformation
  X_clr <- clr_transform(X)
  
  expect_true(is.matrix(X_clr))
  expect_equal(nrow(X_clr), nrow(X))
  expect_equal(ncol(X_clr), ncol(X))
  expect_equal(colnames(X_clr), colnames(X))
  
  # Test CLR properties (row means should be ~0)
  row_means <- rowMeans(X_clr)
  expect_true(all(abs(row_means) < 1e-10))
  
  # Test with custom pseudocount
  X_clr2 <- clr_transform(X, pseudocount = 0.5)
  expect_true(is.matrix(X_clr2))
  expect_equal(dim(X_clr2), dim(X))
})

test_that("melsi basic functionality works", {
  # Generate test data
  test_data <- generate_test_data(n_samples = 40, n_taxa = 50, n_signal_taxa = 5)
  X <- test_data$counts
  y <- test_data$metadata$Group
  
  # CLR transformation
  X_clr <- clr_transform(X)
  
  # Test basic MeLSI run
  expect_no_error({
    results <- melsi(
      X_clr, y, 
      n_perms = 19, 
      B = 10, 
      show_progress = FALSE
    )
  })
  
  # Test results structure
  expect_true(is.list(results))
  expect_true("F_observed" %in% names(results))
  expect_true("p_value" %in% names(results))
  expect_true("feature_weights" %in% names(results))
  expect_true("distance_matrix" %in% names(results))
  
  # Test F-statistic is positive
  expect_true(results$F_observed > 0)
  expect_true(is.numeric(results$F_observed))
  
  # Test p-value is between 0 and 1
  expect_true(results$p_value >= 0 && results$p_value <= 1)
  expect_true(is.numeric(results$p_value))
  
  # Test feature weights (may be fewer due to pre-filtering)
  expect_true(is.numeric(results$feature_weights))
  expect_true(length(results$feature_weights) > 0)
  expect_true(length(results$feature_weights) <= ncol(X_clr))
  
  # Test distance matrix (may be dist object or matrix)
  expect_true(is.matrix(results$distance_matrix) || inherits(results$distance_matrix, "dist"))
  if (is.matrix(results$distance_matrix)) {
    expect_equal(nrow(results$distance_matrix), nrow(X_clr))
    expect_equal(ncol(results$distance_matrix), nrow(X_clr))
  } else {
    # dist object
    expect_equal(attr(results$distance_matrix, "Size"), nrow(X_clr))
  }
})

test_that("plot_vip works correctly", {
  # Generate test data and run MeLSI
  test_data <- generate_test_data(n_samples = 30, n_taxa = 20, n_signal_taxa = 5)
  X <- test_data$counts
  y <- test_data$metadata$Group
  X_clr <- clr_transform(X)
  results <- melsi(X_clr, y, n_perms = 19, B = 10, show_progress = FALSE)
  
  # Test plot_vip with directionality
  expect_no_error({
    p1 <- plot_vip(results, top_n = 10, directionality = TRUE)
  })
  expect_true(inherits(p1, "ggplot"))
  
  # Test plot_vip without directionality
  expect_no_error({
    p2 <- plot_vip(results, top_n = 10, directionality = FALSE)
  })
  expect_true(inherits(p2, "ggplot"))
  
  # Test with custom title
  expect_no_error({
    p3 <- plot_vip(results, title = "Test VIP Plot")
  })
})

test_that("plot_pcoa works correctly", {
  # Generate test data and run MeLSI
  test_data <- generate_test_data(n_samples = 30, n_taxa = 20, n_signal_taxa = 5)
  X <- test_data$counts
  y <- test_data$metadata$Group
  X_clr <- clr_transform(X)
  results <- melsi(X_clr, y, n_perms = 19, B = 10, show_progress = FALSE)
  
  # Test plot_pcoa
  expect_no_error({
    p <- plot_pcoa(results, X_clr, y)
  })
  expect_true(inherits(p, "ggplot"))
  
  # Test with custom title
  expect_no_error({
    p2 <- plot_pcoa(results, X_clr, y, title = "Test PCoA Plot")
  })
})

test_that("plot_feature_importance works correctly", {
  # Generate test data and run MeLSI
  test_data <- generate_test_data(n_samples = 30, n_taxa = 20, n_signal_taxa = 5)
  X <- test_data$counts
  y <- test_data$metadata$Group
  X_clr <- clr_transform(X)
  results <- melsi(X_clr, y, n_perms = 19, B = 10, show_progress = FALSE)
  
  # Test plot_feature_importance
  expect_no_error({
    p <- plot_feature_importance(results$feature_weights, top_n = 10)
  })
  expect_true(inherits(p, "ggplot"))
  
  # Test with directionality
  directionality <- if (!is.null(results$directionality)) results$directionality else NULL
  expect_no_error({
    p2 <- plot_feature_importance(results$feature_weights, top_n = 10, 
                                   directionality = directionality)
  })
})

test_that("melsi handles edge cases", {
  # Test with minimal data
  X_minimal <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, ncol = 2)
  y_minimal <- factor(c("A", "B", "A"))
  X_minimal_clr <- clr_transform(X_minimal)
  
  expect_no_error({
    results <- melsi(X_minimal_clr, y_minimal, n_perms = 9, B = 5, 
                     show_progress = FALSE)
  })
  
  expect_true(results$F_observed > 0)
  expect_true(results$p_value >= 0 && results$p_value <= 1)
})

test_that("melsi input validation works", {
  test_data <- generate_test_data(n_samples = 20, n_taxa = 15, n_signal_taxa = 3)
  X <- test_data$counts
  y <- test_data$metadata$Group
  X_clr <- clr_transform(X)
  
  # Test error with single group
  expect_error({
    melsi(X_clr, rep("A", nrow(X_clr)), n_perms = 9, B = 5, show_progress = FALSE)
  }, "at least 2 groups")
  
  # Test error with mismatched dimensions
  expect_error({
    melsi(X_clr, y[1:(length(y)-1)], n_perms = 9, B = 5, show_progress = FALSE)
  })
})
