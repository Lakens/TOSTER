# Test script comparing TOSTER::perm_t_test with MKinfer::perm.t.test
# Author: Auto-generated for TOSTER package validation
# Purpose: Ensure equivalent results between the two implementations for standard alternatives

library(testthat)
library(MKinfer)
library(TOSTER)
# Source your perm_t_test function or load TOSTER
# source("perm_t_test.R")  # Uncomment if sourcing directly
# library(TOSTER)          # Uncomment if using installed package

# =============================================================================
# IMPORTANT NOTES ON EXPECTED DIFFERENCES:
# =============================================================================
# 1. MKinfer::perm.t.test computes p-values as (b+1)/(R+1) following Phipson & Smyth (2010)
#    where b = number of permutation statistics >= observed (for greater alternative)
#    Your implementation may use b/R or (b+1)/(R+1) - verify and adjust tolerance accordingly
#
# 2. CI methods may differ - MKinfer uses percentile bootstrap on the permutation distribution
#    Your implementation should match this approach
#
# 3. For exact tests (when R=NULL or R >= max_perms), results should be essentially identical
#
# 4. For Monte Carlo tests, set a seed and use the same R for comparable results
#    Even then, small differences are expected due to implementation details
# =============================================================================

# Helper function to compare two htest objects
compare_htest <- function(result1, result2, tol = 0.05,
                          check_ci = TRUE, check_pval = TRUE) {

  # Compare estimates
  expect_equal(unname(result1$estimate), unname(result2$estimate),
               tolerance = 1e-10,
               info = "Estimates should match exactly")

  # Compare test statistics (if both report them)
  if (!is.null(result1$statistic) && !is.null(result2$statistic)) {
    expect_equal(unname(result1$statistic), unname(result2$statistic),
                 tolerance = 1e-10,
                 info = "Test statistics should match exactly")
  }

  # Compare p-values with tolerance for Monte Carlo variation
  if (check_pval) {
    expect_equal(result1$p.value, result2$p.value,
                 tolerance = tol,
                 info = paste("P-values should be similar (tolerance:", tol, ")"))
  }

  # Compare confidence intervals with tolerance
  if (check_ci && !is.null(result1$conf.int) && !is.null(result2$conf.int)) {
    expect_equal(as.numeric(result1$conf.int), as.numeric(result2$conf.int),
                 tolerance = tol,
                 info = paste("CIs should be similar (tolerance:", tol, ")"))
  }

  invisible(TRUE)
}

# =============================================================================
# Test Context: One-sample tests
# =============================================================================
test_that("One-sample permutation t-test matches MKinfer", {

  set.seed(12345)
  x <- rnorm(8, mean = 0.5, sd = 1)
  x

  # Test 1: Two-sided, mu = 0, exact permutation (small sample)
  # With n=15, there are 2^15 = 32768 permutations
  result_toster <- perm_t_test(x, mu = 0, alternative = "two.sided", R = 999)
  result_toster
  result_mkinfer <- perm.t.test(x, mu = 0, alternative = "two.sided", R = 999)
  result_mkinfer
  # For Monte Carlo, allow some tolerance
  compare_htest(result_toster, result_mkinfer, tol = 0.1)

  # Test 2: One-sided greater
  set.seed(12345)
  result_toster <- perm_t_test(x, mu = 0, alternative = "greater", R = 999)
  set.seed(12345)
  result_mkinfer <- perm.t.test(x, mu = 0, alternative = "greater", R = 999)

  compare_htest(result_toster, result_mkinfer, tol = 0.1)

  # Test 3: One-sided less
  set.seed(12345)
  result_toster <- perm_t_test(x, mu = 0, alternative = "less", R = 999)
  set.seed(12345)
  result_mkinfer <- perm.t.test(x, mu = 0, alternative = "less", R = 999)

  compare_htest(result_toster, result_mkinfer, tol = 0.1)

  # Test 4: Small sample - exact permutation
  set.seed(42)
  x_small <- rnorm(8, mean = 1, sd = 1)  # 2^8 = 256 permutations - should be exact

  result_toster <- perm_t_test(x_small, mu = 0, alternative = "two.sided", R = NULL)
  result_mkinfer <- perm.t.test(x_small, mu = 0, alternative = "two.sided", R = 10000)

  # With exact permutation, p-values should be very close
  compare_htest(result_toster, result_mkinfer, tol = 0.02)

})

# =============================================================================
# Test Context: Two-sample independent tests
# =============================================================================
test_that("Two-sample independent permutation t-test matches MKinfer", {

  set.seed(123)
  x <- rnorm(12, mean = 0, sd = 1)
  y <- rnorm(10, mean = 0.5, sd = 1.2)

  # Test 1: Two-sided, Welch (var.equal = FALSE)
  set.seed(999)
  result_toster <- perm_t_test(x, y, alternative = "two.sided", var.equal = FALSE, R = 999)
  set.seed(999)
  result_mkinfer <- perm.t.test(x, y, alternative = "two.sided", var.equal = FALSE, R = 999)

  compare_htest(result_toster, result_mkinfer, tol = 0.1)

  # Test 2: Two-sided, pooled variance (var.equal = TRUE)
  set.seed(999)
  result_toster <- perm_t_test(x, y, alternative = "two.sided", var.equal = TRUE, R = 999)
  set.seed(999)
  result_mkinfer <- perm.t.test(x, y, alternative = "two.sided", var.equal = TRUE, R = 999)

  compare_htest(result_toster, result_mkinfer, tol = 0.1)

  # Test 3: One-sided greater
  set.seed(999)
  result_toster <- perm_t_test(x, y, alternative = "greater", var.equal = FALSE, R = 999)
  set.seed(999)
  result_mkinfer <- perm.t.test(x, y, alternative = "greater", var.equal = FALSE, R = 999)

  compare_htest(result_toster, result_mkinfer, tol = 0.1)

  # Test 4: One-sided less
  set.seed(999)
  result_toster <- perm_t_test(x, y, alternative = "less", var.equal = FALSE, R = 999)
  set.seed(999)
  result_mkinfer <- perm.t.test(x, y, alternative = "less", var.equal = FALSE, R = 999)

  compare_htest(result_toster, result_mkinfer, tol = 0.1)

  # Test 5: Small samples - exact permutation
  set.seed(42)
  x_small <- rnorm(5, mean = 0, sd = 1)
  y_small <- rnorm(5, mean = 1, sd = 1)
  # choose(10, 5) = 252 permutations

  result_toster <- perm_t_test(x_small, y_small, alternative = "two.sided", R = NULL)
  result_mkinfer <- perm.t.test(x_small, y_small, alternative = "two.sided", R = 10000)

  # With exact permutation, should be very close
  compare_htest(result_toster, result_mkinfer, tol = 0.02)

})

# =============================================================================
# Test Context: Paired samples tests
# =============================================================================
test_that("Paired permutation t-test matches MKinfer", {

  set.seed(456)
  x <- rnorm(12, mean = 5, sd = 1)
  y <- x + rnorm(12, mean = 0.3, sd = 0.5)  # Paired with small positive effect

  # Test 1: Two-sided
  set.seed(777)
  result_toster <- perm_t_test(x, y, paired = TRUE, alternative = "two.sided", R = 999)
  set.seed(777)
  result_mkinfer <- perm.t.test(x, y, paired = TRUE, alternative = "two.sided", R = 999)

  compare_htest(result_toster, result_mkinfer, tol = 0.1)

  # Test 2: One-sided (expect y > x, so x - y should be negative, testing "less")
  set.seed(777)
  result_toster <- perm_t_test(x, y, paired = TRUE, alternative = "less", R = 999)
  set.seed(777)
  result_mkinfer <- perm.t.test(x, y, paired = TRUE, alternative = "less", R = 999)

  compare_htest(result_toster, result_mkinfer, tol = 0.1)

  # Test 3: Small paired sample - exact
  set.seed(42)
  x_small <- rnorm(8, mean = 10, sd = 2)
  y_small <- x_small + rnorm(8, mean = 1, sd = 1)
  # 2^8 = 256 permutations

  result_toster <- perm_t_test(x_small, y_small, paired = TRUE, alternative = "two.sided", R = NULL)
  result_mkinfer <- perm.t.test(x_small, y_small, paired = TRUE, alternative = "two.sided", R = 10000)

  compare_htest(result_toster, result_mkinfer, tol = 0.02)

})

# =============================================================================
# Test Context: Formula interface
# =============================================================================
test_that("Formula interface produces equivalent results to x, y interface", {

  data(sleep)

  # Using formula
  set.seed(111)
  result_formula <- perm_t_test(extra ~ group, data = sleep, R = 999)

  # Using x, y vectors
  x <- sleep$extra[sleep$group == 1]
  y <- sleep$extra[sleep$group == 2]
  set.seed(111)
  result_xy <- perm_t_test(x, y, R = 999)

  # Should produce identical results
  expect_equal(result_formula$statistic, result_xy$statistic, tolerance = 1e-10)
  expect_equal(result_formula$p.value, result_xy$p.value, tolerance = 1e-10)
  expect_equal(as.numeric(result_formula$conf.int), as.numeric(result_xy$conf.int), tolerance = 1e-10)

})

# =============================================================================
# Test Context: Edge cases and special situations
# =============================================================================
test_that("Edge cases are handled correctly", {

  # Test 1: Equal groups (no difference)
  set.seed(42)
  x <- rnorm(20, mean = 0, sd = 1)
  y <- rnorm(20, mean = 0, sd = 1)

  result <- perm_t_test(x, y, R = 999)
  expect_true(result$p.value > 0.05)  # Should not reject when no difference

  # Test 2: Large difference (should reject)
  set.seed(42)
  x <- rnorm(20, mean = 0, sd = 1)
  y <- rnorm(20, mean = 3, sd = 1)

  result <- perm_t_test(x, y, R = 999)
  expect_true(result$p.value < 0.05)  # Should reject with large difference

  # Test 3: Non-zero mu in one-sample
  set.seed(42)
  x <- rnorm(20, mean = 5, sd = 1)

  result <- perm_t_test(x, mu = 5, R = 999)
  expect_true(result$p.value > 0.05)  # Should not reject when mu matches true mean

  # Test 4: Handling of NA values
  x_na <- c(1, 2, NA, 4, 5, 6, 7, 8, 9, 10)
  result <- perm_t_test(x_na, mu = 0, R = 499)
  expect_equal(length(result$estimate), 1)  # Should run without error

})

# =============================================================================
# Test Context: P-value computation method verification
# =============================================================================
test_that("P-value computation follows Phipson & Smyth (2010) recommendation", {

  # The p-value should never be exactly 0 when using (b+1)/(R+1)
  # This is the recommended formula to avoid zero p-values

  set.seed(42)
  x <- rnorm(10, mean = 5, sd = 0.1)  # Very extreme data

  result <- perm_t_test(x, mu = 0, R = 99)

  # P-value should be at least 1/R_used (your implementation uses max(pval, 1/R_used))
  expect_true(result$p.value >= 1/result$R.used,
              info = "P-value should never be less than 1/R to avoid zero p-values")

  # Also check that p-value is not exactly zero
  expect_true(result$p.value > 0,
              info = "P-value should never be exactly zero")

})

# =============================================================================
# Test Context: Symmetric vs non-symmetric two-sided p-value
# =============================================================================
test_that("Symmetric and non-symmetric p-value options work correctly", {

  set.seed(42)
  x <- c(1, 2, 3, 4, 5, 6, 7, 20)  # Skewed data

  # Symmetric (default): p = mean(|T*| >= |T|)
  result_sym <- perm_t_test(x, mu = 0, alternative = "two.sided",
                            symmetric = TRUE, R = 999)


  # Non-symmetric: p = 2 * min(P(T* <= T), P(T* > T))
  result_nonsym <- perm_t_test(x, mu = 0, alternative = "two.sided",
                               symmetric = FALSE, R = 999)

  # Both should return valid p-values
  expect_true(result_sym$p.value >= 0 && result_sym$p.value <= 1)
  expect_true(result_nonsym$p.value >= 0 && result_nonsym$p.value <= 1)

  # They may differ, especially with skewed data
  # Just verify both work and produce reasonable results

})

# =============================================================================
# Test Context: Degrees of freedom
# =============================================================================
test_that("Degrees of freedom match standard t-test", {

  set.seed(42)
  x <- rnorm(15, mean = 0, sd = 1)
  y <- rnorm(12, mean = 0, sd = 1.5)

  # Compare df with standard t.test
  # One-sample
  result_perm_1s <- perm_t_test(x, mu = 0, R = 99)
  result_t_1s <- t.test(x, mu = 0)
  expect_equal(result_perm_1s$parameter, result_t_1s$parameter, tolerance = 1e-10)

  # Two-sample Welch
  result_perm_2s <- perm_t_test(x, y, var.equal = FALSE, R = 99)
  result_t_2s <- t.test(x, y, var.equal = FALSE)
  expect_equal(unname(result_perm_2s$parameter), unname(result_t_2s$parameter), tolerance = 1e-10)

  # Two-sample pooled
  result_perm_pool <- perm_t_test(x, y, var.equal = TRUE, R = 99)
  result_t_pool <- t.test(x, y, var.equal = TRUE)
  expect_equal(unname(result_perm_pool$parameter), unname(result_t_pool$parameter), tolerance = 1e-10)

})

# =============================================================================
# Test Context: Test statistics match parametric t-test
# =============================================================================
test_that("Test statistics match standard t-test", {

  set.seed(42)
  x <- rnorm(20, mean = 0.5, sd = 1)
  y <- rnorm(18, mean = 0, sd = 1.2)

  # One-sample
  result_perm <- perm_t_test(x, mu = 0, R = 99)
  result_t <- t.test(x, mu = 0)
  expect_equal(unname(result_perm$statistic), unname(result_t$statistic), tolerance = 1e-10)

  # Two-sample
  result_perm <- perm_t_test(x, y, R = 99)
  result_t <- t.test(x, y)
  expect_equal(unname(result_perm$statistic), unname(result_t$statistic), tolerance = 1e-10)

  # Paired
  set.seed(42)
  x_p <- rnorm(15, 10, 2)
  y_p <- x_p + rnorm(15, 0.5, 1)

  result_perm <- perm_t_test(x_p, y_p, paired = TRUE, R = 99)
  result_t <- t.test(x_p, y_p, paired = TRUE)
  expect_equal(unname(result_perm$statistic), unname(result_t$statistic), tolerance = 1e-10)

})

# =============================================================================
# Test Context: Exact permutation consistency
# =============================================================================
test_that("Exact permutation gives deterministic results", {

  # Very small sample for exact enumeration
  set.seed(42)
  x <- rnorm(6)
  y <- rnorm(6)

  # Run twice - should get identical results with exact permutation
  result1 <- perm_t_test(x, y, R = NULL)  # Exact

result2 <- perm_t_test(x, y, R = NULL)  # Exact

  expect_equal(result1$p.value, result2$p.value, tolerance = 1e-15)
  expect_equal(as.numeric(result1$conf.int), as.numeric(result2$conf.int), tolerance = 1e-15)

})

# =============================================================================
# Test Context: Monte Carlo reproducibility with seed
# =============================================================================
test_that("Monte Carlo results are reproducible with same seed", {

  set.seed(42)
  x <- rnorm(30)
  y <- rnorm(30)

  # Run with same seed - should get identical results
  set.seed(123)
  result1 <- perm_t_test(x, y, R = 999)

  set.seed(123)
  result2 <- perm_t_test(x, y, R = 999)

  expect_equal(result1$p.value, result2$p.value, tolerance = 1e-15)
  expect_equal(as.numeric(result1$conf.int), as.numeric(result2$conf.int), tolerance = 1e-15)

})

# =============================================================================
# Test Context: Return structure matches htest class requirements
# =============================================================================
test_that("Return object has correct htest structure", {

  set.seed(42)
  x <- rnorm(20)
  y <- rnorm(20)

  result <- perm_t_test(x, y, R = 99)

  # Check class
  expect_s3_class(result, "htest")

  # Check required components
  expect_true("statistic" %in% names(result))
  expect_true("parameter" %in% names(result))
  expect_true("p.value" %in% names(result))
  expect_true("conf.int" %in% names(result))
  expect_true("estimate" %in% names(result))
  expect_true("null.value" %in% names(result))
  expect_true("alternative" %in% names(result))
  expect_true("method" %in% names(result))
  expect_true("data.name" %in% names(result))

  # Check that CI has conf.level attribute
  expect_true(!is.null(attr(result$conf.int, "conf.level")))

})

# =============================================================================
# Test Context: Comparison with specific MKinfer behavior
# =============================================================================
test_that("Specific comparison scenarios with MKinfer", {

  # Use the sleep dataset which is classic for t-tests
  data(sleep)

  # Independent samples test
  x <- sleep$extra[sleep$group == 1]
  y <- sleep$extra[sleep$group == 2]

  set.seed(2024)
  result_toster <- perm_t_test(x, y, alternative = "two.sided", R = 1999)
  set.seed(2024)
  result_mkinfer <- perm.t.test(x, y, alternative = "two.sided", R = 1999)

  # Print comparison for diagnostic purposes
  cat("\n=== Sleep data comparison ===\n")
  cat("TOSTER estimate:", result_toster$estimate, "\n")
  cat("MKinfer estimate:", result_mkinfer$estimate, "\n")
  cat("TOSTER t-stat:", result_toster$statistic, "\n")
  cat("MKinfer t-stat:", result_mkinfer$statistic, "\n")
  cat("TOSTER p-value:", result_toster$p.value, "\n")
  cat("MKinfer p-value:", result_mkinfer$p.value, "\n")
  cat("TOSTER CI:", result_toster$conf.int, "\n")
  cat("MKinfer CI:", result_mkinfer$conf.int, "\n")

  # Estimates and statistics should match exactly
  expect_equal(unname(result_toster$estimate), unname(result_mkinfer$estimate),
               tolerance = 1e-10)
  expect_equal(unname(result_toster$statistic), unname(result_mkinfer$statistic),
               tolerance = 1e-10)

  # P-values and CIs with some tolerance for Monte Carlo
  expect_equal(result_toster$p.value, result_mkinfer$p.value, tolerance = 0.05)

})

# =============================================================================
# Run all tests
# =============================================================================
cat("\n========================================\n")
cat("Running perm_t_test vs MKinfer comparison tests\n
")
cat("========================================\n\n")

# Run tests
test_results <- test_file(test_path = ".", reporter = "summary")

cat("\n========================================\n")
cat("Test Summary Complete\n")
cat("========================================\n")
