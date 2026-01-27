# Unit tests for hodges_lehmann function
# Run with: testthat::test_file("test-hodges_lehmann.R")

context("hodges_lehmann")

# =============================================================================
# Test Data Setup
# =============================================================================

set.seed(42)
x_sample <- rnorm(15, mean = 5, sd = 2)
y_sample <- rnorm(12, mean = 6, sd = 2)
paired_x <- c(5.1, 4.8, 6.2, 5.7, 6.0, 5.5, 4.9, 5.8)
paired_y <- c(5.6, 5.2, 6.7, 6.1, 6.5, 5.8, 5.3, 6.2)

# =============================================================================
# Basic Structure and Class Tests
# =============================================================================

test_that("hodges_lehmann returns htest object with all expected components", {
  set.seed(123)
  result <- hodges_lehmann(x_sample, y_sample, R = 199)

  expect_s3_class(result, "htest")

  # Check all required htest components
  expect_true("statistic" %in% names(result))
  expect_true("p.value" %in% names(result))
  expect_true("conf.int" %in% names(result))
  expect_true("estimate" %in% names(result))
  expect_true("null.value" %in% names(result))
  expect_true("alternative" %in% names(result))
  expect_true("method" %in% names(result))
  expect_true("data.name" %in% names(result))

  # Check hodges_lehmann specific components
  expect_true("call" %in% names(result))
  expect_true("R" %in% names(result))
  expect_true("R.used" %in% names(result))
})

test_that("hodges_lehmann keeps permutation distribution when keep_perm = TRUE", {
  set.seed(123)
  result_keep <- hodges_lehmann(x_sample, y_sample, R = 99, keep_perm = TRUE)
  result_no_keep <- hodges_lehmann(x_sample, y_sample, R = 99, keep_perm = FALSE)

  expect_true("perm.stat" %in% names(result_keep))
  expect_true("perm.eff" %in% names(result_keep))
  expect_false("perm.stat" %in% names(result_no_keep))
  expect_false("perm.eff" %in% names(result_no_keep))
})

# =============================================================================
# Hodges-Lehmann Estimator Tests
# =============================================================================

test_that("HL1 estimator matches wilcox.test pseudomedian for one-sample", {
  set.seed(123)
  x <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)

  # Get wilcox.test estimate
  wt_result <- wilcox.test(x, conf.int = TRUE)

  # Get hodges_lehmann estimate (asymptotic)
  hl_result <- hodges_lehmann(x, mu = 0)

  # Estimates should match

  expect_equal(as.numeric(hl_result$estimate),
               as.numeric(wt_result$estimate),
               tolerance = 1e-10)
})

test_that("HL2 estimator matches wilcox.test estimate for two-sample", {
  set.seed(123)
  x <- c(1, 2, 3, 4, 5)
  y <- c(6, 7, 8, 9, 10)

  # Get wilcox.test estimate
  wt_result <- wilcox.test(x, y, conf.int = TRUE)

  # Get hodges_lehmann estimate (asymptotic)
  hl_result <- hodges_lehmann(x, y)

  # Estimates should match exactly (both compute x - y)
  expect_equal(as.numeric(hl_result$estimate),
               as.numeric(wt_result$estimate),
               tolerance = 1e-10)
})

test_that("HL1 estimator is median of Walsh averages", {
  x <- c(1, 3, 5, 7)

  # Manual calculation of Walsh averages (i <= j)
  # Pairs: (1,1)/2=1, (1,3)/2=2, (1,5)/2=3, (1,7)/2=4,
  #        (3,3)/2=3, (3,5)/2=4, (3,7)/2=5,
  #        (5,5)/2=5, (5,7)/2=6,
  #        (7,7)/2=7
  walsh <- c(1, 2, 3, 4, 3, 4, 5, 5, 6, 7)
  expected <- median(walsh)

  result <- hodges_lehmann(x, mu = 0)

  expect_equal(as.numeric(result$estimate), expected, tolerance = 1e-10)
})

test_that("HL2 estimator is median of pairwise differences", {
  x <- c(1, 2)
  y <- c(5, 6)

  # Manual calculation: x_i - y_j (matches wilcox.test sign convention)
  # 1-5=-4, 1-6=-5, 2-5=-3, 2-6=-4
  diffs <- c(-4, -5, -3, -4)
  expected <- median(diffs)

  result <- hodges_lehmann(x, y)

  expect_equal(as.numeric(result$estimate), expected, tolerance = 1e-10)
})

# =============================================================================
# One-Sample Tests
# =============================================================================

test_that("one-sample hodges_lehmann works correctly (asymptotic)", {
  set.seed(123)
  result <- hodges_lehmann(x_sample, mu = 5)

  expect_s3_class(result, "htest")
  expect_true(grepl("One Sample", result$method))
  expect_true(grepl("Asymptotic", result$method))
  expect_equal(result$null.value, c("location" = 5))
  expect_named(result$estimate, "pseudomedian")

  # P-value should be in valid range
  expect_gte(result$p.value, 0)
  expect_lte(result$p.value, 1)

  # Confidence interval should have correct structure
  expect_length(result$conf.int, 2)
})

test_that("one-sample hodges_lehmann works correctly (permutation)", {
  set.seed(123)
  result <- hodges_lehmann(x_sample, mu = 5, R = 199)

  expect_s3_class(result, "htest")
  expect_true(grepl("One Sample", result$method))
  expect_true(grepl("Randomization|Permutation", result$method))
})

test_that("one-sample hodges_lehmann handles mu = 0 correctly", {
  set.seed(123)
  # Generate data with mean noticeably different from 0
  x_nonzero <- rnorm(20, mean = 3, sd = 1)
  result <- hodges_lehmann(x_nonzero, mu = 0, R = 499)

  expect_s3_class(result, "htest")
  # Should have small p-value since location is clearly not 0
  expect_lt(result$p.value, 0.05)
})

# =============================================================================
# Two-Sample Tests
# =============================================================================

test_that("two-sample hodges_lehmann works with asymptotic method", {
  set.seed(123)
  result <- hodges_lehmann(x_sample, y_sample)

  expect_s3_class(result, "htest")
  expect_true(grepl("Two Sample", result$method))
  expect_true(grepl("Asymptotic", result$method))
  expect_named(result$estimate, "difference in location")
})

test_that("two-sample hodges_lehmann works with permutation method", {
  set.seed(123)
  result <- hodges_lehmann(x_sample, y_sample, R = 199)

  expect_s3_class(result, "htest")
  expect_true(grepl("Two Sample", result$method))
  expect_true(grepl("Randomization|Permutation", result$method))
})

test_that("two-sample hodges_lehmann works with formula interface", {
  set.seed(123)
  # Using built-in sleep data
  result <- hodges_lehmann(extra ~ group, data = sleep, R = 199)

  expect_s3_class(result, "htest")
  expect_equal(result$data.name, "extra by group")
})

# =============================================================================
# Paired Tests
# =============================================================================

test_that("paired hodges_lehmann works correctly", {
  set.seed(123)
  result <- hodges_lehmann(paired_x, paired_y, paired = TRUE, R = 199)

  expect_s3_class(result, "htest")
  expect_true(grepl("Paired", result$method))
  expect_named(result$estimate, "pseudomedian of differences")
})

test_that("paired hodges_lehmann detects significant difference", {
  set.seed(123)
  # These paired samples have a consistent positive difference
  before <- c(5, 6, 7, 8, 9, 10, 11, 12)
  after <- c(5.8, 7.1, 7.9, 9.2, 9.8, 11.1, 11.9, 13.2)

  result <- hodges_lehmann(before, after, paired = TRUE,
                           alternative = "less", R = 499)

  expect_s3_class(result, "htest")
  # before < after consistently, so "less" should be significant
  expect_lt(result$p.value, 0.05)
})

# =============================================================================
# Alternative Hypotheses
# =============================================================================

test_that("alternative = 'two.sided' works correctly", {
  set.seed(123)
  result <- hodges_lehmann(x_sample, y_sample, alternative = "two.sided", R = 199)

  expect_equal(result$alternative, "two.sided")
  expect_length(result$conf.int, 2)
  expect_true(all(is.finite(result$conf.int)))
})

test_that("alternative = 'less' works correctly", {
  set.seed(123)
  result <- hodges_lehmann(x_sample, y_sample, alternative = "less", R = 199)

  expect_equal(result$alternative, "less")
  expect_equal(result$conf.int[1], -Inf)
  expect_true(is.finite(result$conf.int[2]))
})

test_that("alternative = 'greater' works correctly", {
  set.seed(123)
  result <- hodges_lehmann(x_sample, y_sample, alternative = "greater", R = 199)

  expect_equal(result$alternative, "greater")
  expect_true(is.finite(result$conf.int[1]))
  expect_equal(result$conf.int[2], Inf)
})

test_that("alternative = 'equivalence' works correctly", {
  set.seed(123)
  result <- hodges_lehmann(x_sample, y_sample,
                           alternative = "equivalence",
                           mu = c(-2, 2), R = 199)

  expect_equal(result$alternative, "equivalence")
  expect_length(result$null.value, 2)
  expect_equal(as.numeric(result$null.value), c(-2, 2))

  # Confidence level should be 1 - 2*alpha = 0.90 for equivalence
  expect_equal(attr(result$conf.int, "conf.level"), 0.90)
})

test_that("alternative = 'minimal.effect' works correctly", {
  set.seed(123)
  result <- hodges_lehmann(x_sample, y_sample,
                           alternative = "minimal.effect",
                           mu = c(-2, 2), R = 199)

  expect_equal(result$alternative, "minimal.effect")
  expect_length(result$null.value, 2)
  expect_equal(as.numeric(result$null.value), c(-2, 2))
})

test_that("equivalence with single mu value creates symmetric bounds", {
  set.seed(123)
  result <- hodges_lehmann(x_sample, y_sample,
                           alternative = "equivalence",
                           mu = 2, R = 199)

  expect_equal(as.numeric(result$null.value), c(-2, 2))
})

# =============================================================================
# Scale Estimator Tests
# =============================================================================

test_that("scale = 'S1' works correctly", {
  set.seed(123)
  result <- hodges_lehmann(x_sample, y_sample, R = 199, scale = "S1")

  expect_s3_class(result, "htest")
  expect_gte(result$p.value, 0)
  expect_lte(result$p.value, 1)
})

test_that("scale = 'S2' works correctly", {
  set.seed(123)
  result <- hodges_lehmann(x_sample, y_sample, R = 199, scale = "S2")

  expect_s3_class(result, "htest")
  expect_gte(result$p.value, 0)
  expect_lte(result$p.value, 1)
})

test_that("different scale estimators can produce different results", {
  set.seed(123)
  result_s1 <- hodges_lehmann(x_sample, y_sample, R = 199, scale = "S1")
  set.seed(123)
  result_s2 <- hodges_lehmann(x_sample, y_sample, R = 199, scale = "S2")

  # Estimates should be the same
  expect_equal(result_s1$estimate, result_s2$estimate)

  # P-values may differ
  expect_gte(result_s1$p.value, 0)
  expect_gte(result_s2$p.value, 0)
})

# =============================================================================
# P-value Method Tests
# =============================================================================

test_that("p_method = 'plusone' produces valid p-values", {
  set.seed(123)
  result <- hodges_lehmann(x_sample, y_sample, R = 199, p_method = "plusone")

  expect_gt(result$p.value, 0)  # Should never be exactly 0
  expect_lte(result$p.value, 1)
})

test_that("p_method = 'exact' produces valid p-values", {
  set.seed(123)
  result <- hodges_lehmann(x_sample, y_sample, R = 199, p_method = "exact")

  expect_gte(result$p.value, 0)
  expect_lte(result$p.value, 1)
})

test_that("p_method defaults to 'exact' for exact permutation", {
  set.seed(123)
  small_x <- c(1, 2, 3, 4)
  small_y <- c(5, 6, 7)

  # Exact permutation (choose(7, 4) = 35)
  expect_message(
    result <- hodges_lehmann(small_x, small_y, R = 1000),
    "exact permutations"
  )

  # P-value could be exactly 0 with "exact" method
  expect_gte(result$p.value, 0)
})

# =============================================================================
# Exact Permutation Tests
# =============================================================================

test_that("exact permutation is computed for small two-sample", {
  set.seed(123)
  small_x <- c(1, 2, 3, 4)
  small_y <- c(5, 6, 7)

  # choose(7, 4) = 35 permutations
  expect_message(
    result <- hodges_lehmann(small_x, small_y, R = 1000),
    "exact permutations"
  )

  expect_s3_class(result, "htest")
  expect_true(grepl("Exact", result$method))
  expect_equal(result$R.used, choose(7, 4))
})

test_that("exact permutation for one-sample small samples", {
  set.seed(123)
  small_x <- c(1, 2, 3, 4, 5)

  # 2^5 = 32 sign permutations
  expect_message(
    result <- hodges_lehmann(small_x, mu = 0, R = 1000),
    "exact permutations"
  )

  expect_true(grepl("Exact", result$method))
  expect_equal(result$R.used, 2^5)
})

test_that("Randomization is used when R is specified and smaller than max perms", {
  set.seed(123)
  # Large enough sample that R = 199 triggers Randomization
  result <- hodges_lehmann(x_sample, y_sample, R = 199)

  expect_true(grepl("Randomization", result$method))
  expect_equal(result$R, 199)
})

# =============================================================================
# Asymptotic Test
# =============================================================================

test_that("asymptotic test is used when R = NULL", {
  set.seed(123)
  result <- hodges_lehmann(x_sample, y_sample, R = NULL)

  expect_true(grepl("Asymptotic", result$method))
  expect_null(result$R)
  expect_true(is.na(result$R.used))
})

test_that("asymptotic test produces valid p-values", {
  set.seed(123)
  result <- hodges_lehmann(x_sample, y_sample)

  expect_gte(result$p.value, 0)
  expect_lte(result$p.value, 1)
})

test_that("asymptotic CIs are finite for two-sided test", {
  set.seed(123)
  result <- hodges_lehmann(x_sample, y_sample, alternative = "two.sided")

  expect_true(all(is.finite(result$conf.int)))
})

# =============================================================================
# Missing Value Handling
# =============================================================================

test_that("hodges_lehmann handles NA values correctly", {
  set.seed(123)
  x_na <- c(x_sample, NA, NA)
  y_na <- c(NA, y_sample, NA)

  result <- hodges_lehmann(x_na, y_na, R = 199)

  expect_s3_class(result, "htest")
})

test_that("paired test handles NA values correctly", {
  set.seed(123)
  x_na <- c(paired_x, NA)
  y_na <- c(paired_y, NA)

  result <- hodges_lehmann(x_na, y_na, paired = TRUE, R = 199)

  expect_s3_class(result, "htest")
})

# =============================================================================
# Error Handling
# =============================================================================

test_that("error for non-numeric x", {
  expect_error(hodges_lehmann(c("a", "b", "c")),
               "'x' must be numeric")
})

test_that("error for non-numeric y", {
  expect_error(hodges_lehmann(x_sample, c("a", "b", "c")),
               "'y' must be numeric")
})

test_that("error for invalid alpha", {
  expect_error(hodges_lehmann(x_sample, y_sample, alpha = -0.1),
               "'alpha' must be a single number between 0 and 1")
  expect_error(hodges_lehmann(x_sample, y_sample, alpha = 1.5),
               "'alpha' must be a single number between 0 and 1")
})

test_that("error for invalid R", {
  expect_error(hodges_lehmann(x_sample, y_sample, R = 0),
               "'R' must be NULL .* or a positive integer")
  expect_error(hodges_lehmann(x_sample, y_sample, R = -10),
               "'R' must be NULL .* or a positive integer")
})

test_that("error for mu = 0 with equivalence/minimal.effect", {
  expect_error(hodges_lehmann(x_sample, y_sample, alternative = "equivalence", mu = 0),
               "For equivalence or minimal.effect testing")
})

test_that("error for wrong mu length with standard alternatives", {
  expect_error(hodges_lehmann(x_sample, y_sample, alternative = "two.sided", mu = c(1, 2)),
               "'mu' must be a single finite number for this alternative")
})

test_that("error for unequal lengths in paired test", {
  expect_error(hodges_lehmann(c(1, 2, 3), c(4, 5), paired = TRUE),
               "'x' and 'y' must have the same length for paired test")
})

test_that("error for insufficient sample size", {
  expect_error(hodges_lehmann(c(1), mu = 0),
               "need at least 2 observations")
})

test_that("error for insufficient sample size in each group", {
  expect_error(hodges_lehmann(c(1), c(2, 3, 4)),
               "need at least 2 observations in each group")
})

test_that("error for formula with wrong number of groups", {
  df <- data.frame(value = 1:9, group = rep(c("a", "b", "c"), 3))
  expect_error(hodges_lehmann(value ~ group, data = df),
               "grouping factor must have exactly 2 levels")
})

test_that("error for incorrect formula", {
  expect_error(hodges_lehmann(~ group, data = sleep),
               "'formula' missing or incorrect")
})

# =============================================================================
# Reproducibility Tests
# =============================================================================

test_that("results are reproducible with set.seed for permutation test", {
  set.seed(42)
  result1 <- hodges_lehmann(x_sample, y_sample, R = 199)
  set.seed(42)
  result2 <- hodges_lehmann(x_sample, y_sample, R = 199)

  expect_equal(result1$p.value, result2$p.value)
  expect_equal(result1$conf.int, result2$conf.int)
  expect_equal(result1$statistic, result2$statistic)
})

test_that("asymptotic test is deterministic", {
  result1 <- hodges_lehmann(x_sample, y_sample)
  result2 <- hodges_lehmann(x_sample, y_sample)

  expect_equal(result1$p.value, result2$p.value)
  expect_equal(result1$conf.int, result2$conf.int)
  expect_equal(result1$statistic, result2$statistic)
})

# =============================================================================
# Confidence Interval Tests
# =============================================================================

test_that("confidence level attribute is correct", {
  result_two <- hodges_lehmann(x_sample, y_sample, alternative = "two.sided",
                               alpha = 0.05, R = 199)
  expect_equal(attr(result_two$conf.int, "conf.level"), 0.95)

  result_less <- hodges_lehmann(x_sample, y_sample, alternative = "less",
                                alpha = 0.05, R = 199)
  expect_equal(attr(result_less$conf.int, "conf.level"), 0.95)

  result_equiv <- hodges_lehmann(x_sample, y_sample, alternative = "equivalence",
                                 mu = c(-3, 3), alpha = 0.05, R = 199)
  expect_equal(attr(result_equiv$conf.int, "conf.level"), 0.90)
})

test_that("different alpha produces appropriate confidence intervals", {
  set.seed(123)

  result_95 <- hodges_lehmann(x_sample, y_sample, alpha = 0.05, R = 299)
  result_99 <- hodges_lehmann(x_sample, y_sample, alpha = 0.01, R = 299)

  # 99% CI (alpha=0.01) should be wider than 95% CI (alpha=0.05)
  width_95 <- result_95$conf.int[2] - result_95$conf.int[1]
  width_99 <- result_99$conf.int[2] - result_99$conf.int[1]

  expect_gt(width_99, width_95)
})

# =============================================================================
# Method String Tests
# =============================================================================

test_that("method string reflects test type correctly", {
  # Two-sample asymptotic
  result1 <- hodges_lehmann(x_sample, y_sample)
  expect_match(result1$method, "Asymptotic")
  expect_match(result1$method, "Two Sample")

  # Two-sample permutation
  result2 <- hodges_lehmann(x_sample, y_sample, R = 199)
  expect_match(result2$method, "Randomization|Permutation")
  expect_match(result2$method, "Two Sample")

  # One-sample
  result3 <- hodges_lehmann(x_sample, mu = 5)
  expect_match(result3$method, "One Sample")

  # Paired
  result4 <- hodges_lehmann(paired_x, paired_y, paired = TRUE, R = 199)
  expect_match(result4$method, "Paired")
})

# =============================================================================
# Permutation Distribution Tests
# =============================================================================

test_that("permutation distribution has correct length", {
  set.seed(123)
  R <- 199
  result <- hodges_lehmann(x_sample, y_sample, R = R, keep_perm = TRUE)

  expect_length(result$perm.stat, result$R.used)
  expect_length(result$perm.eff, result$R.used)
})

test_that("permutation distribution is numeric", {
  set.seed(123)
  result <- hodges_lehmann(x_sample, y_sample, R = 199, keep_perm = TRUE)

  expect_type(result$perm.stat, "double")
  expect_type(result$perm.eff, "double")
  expect_true(all(is.finite(result$perm.stat)))
  expect_true(all(is.finite(result$perm.eff)))
})

# =============================================================================
# Edge Cases
# =============================================================================

test_that("handles samples with equal values", {
  set.seed(123)
  x_const <- rep(5, 10)

  # Should work but may produce warnings about scale
  expect_warning(
    result <- hodges_lehmann(x_const, mu = 5, R = 199),
    "Scale estimate"
  )
})

test_that("handles samples with different variances", {
  set.seed(123)
  x_low_var <- c(4, 5, 6)
  y_high_var <- c(2, 5, 8)  # Same median, different variance

  result <- hodges_lehmann(x_low_var, y_high_var, R = 199)

  expect_s3_class(result, "htest")
})

test_that("handles very small samples for paired test", {
  set.seed(123)
  x_tiny <- c(1, 2, 3)
  y_tiny <- c(2.1, 2.8, 4.2)

  # 2^3 = 8 exact permutations
  expect_message(
    result <- hodges_lehmann(x_tiny, y_tiny, paired = TRUE, R = 1000),
    "exact permutations"
  )

  expect_s3_class(result, "htest")
})

# =============================================================================
# Comparison with wilcox.test Direction Consistency
# =============================================================================

test_that("direction of effect matches wilcox.test", {
  set.seed(123)
  # Clear difference: x < y
  x_low <- c(1, 2, 3, 4, 5)
  y_high <- c(10, 11, 12, 13, 14)

  hl_result <- hodges_lehmann(x_low, y_high)
  wt_result <- wilcox.test(x_low, y_high, conf.int = TRUE)

  # Estimates should match exactly (both compute x - y)
  expect_equal(as.numeric(hl_result$estimate),
               as.numeric(wt_result$estimate),
               tolerance = 1e-10)
})

test_that("one-sided tests give appropriate p-values", {
  set.seed(123)
  # x clearly less than y
  x_low <- c(1, 2, 3, 4, 5)
  y_high <- c(10, 11, 12, 13, 14)

  result_less <- hodges_lehmann(x_low, y_high, alternative = "less", R = 499)
  result_greater <- hodges_lehmann(x_low, y_high, alternative = "greater", R = 499)

  # hodges_lehmann computes x - y, which is negative when x < y
  # "less" tests H1: location shift < 0, should be highly significant
  expect_lt(result_less$p.value, 0.05)
  # "greater" tests H1: location shift > 0, should not be significant
  expect_gt(result_greater$p.value, 0.5)
})

# =============================================================================
# Print Method Test
# =============================================================================

test_that("print method works without error", {
  result <- hodges_lehmann(x_sample, y_sample, R = 99)

  # Should print without error
  expect_output(print(result))
})
