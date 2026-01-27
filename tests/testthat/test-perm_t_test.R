# Unit tests for perm_t_test function
# Run with: testthat::test_file("test-perm_t_test.R")

context("perm_t_test")


# Test Data Setup -----------------


set.seed(42)
x_sample <- rnorm(15, mean = 5, sd = 2)
y_sample <- rnorm(12, mean = 6, sd = 2)
paired_x <- c(5.1, 4.8, 6.2, 5.7, 6.0, 5.5, 4.9, 5.8)
paired_y <- c(5.6, 5.2, 6.7, 6.1, 6.5, 5.8, 5.3, 6.2)


# Basic Structure and Class Tests -----------------


test_that("perm_t_test returns htest object with all expected components", {
  set.seed(123)
  result <- perm_t_test(x_sample, y_sample, R = 199)

  expect_s3_class(result, "htest")

  # Check all required htest components
  expect_true("statistic" %in% names(result))
  expect_true("parameter" %in% names(result))
  expect_true("p.value" %in% names(result))
  expect_true("conf.int" %in% names(result))
  expect_true("estimate" %in% names(result))
  expect_true("null.value" %in% names(result))
  expect_true("alternative" %in% names(result))
  expect_true("method" %in% names(result))
  expect_true("data.name" %in% names(result))

  # Check additional perm_t_test specific components
  expect_true("stderr" %in% names(result))
  expect_true("call" %in% names(result))
  expect_true("R" %in% names(result))
  expect_true("R.used" %in% names(result))
})

test_that("perm_t_test keeps permutation distribution when keep_perm = TRUE", {
  set.seed(123)
  result_keep <- perm_t_test(x_sample, y_sample, R = 99, keep_perm = TRUE)
  result_no_keep <- perm_t_test(x_sample, y_sample, R = 99, keep_perm = FALSE)

  expect_true("perm.stat" %in% names(result_keep))
  expect_true("perm.eff" %in% names(result_keep))
  expect_false("perm.stat" %in% names(result_no_keep))
  expect_false("perm.eff" %in% names(result_no_keep))
})


# One-Sample Tests -----------------


test_that("one-sample perm_t_test works correctly", {
  set.seed(123)
  result <- perm_t_test(x_sample, mu = 5, R = 199)

  expect_s3_class(result, "htest")
  expect_true(grepl("One Sample", result$method))
  expect_equal(result$null.value, c("mean" = 5))
  expect_named(result$estimate, "mean of x")

  # P-value should be in valid range

  expect_gte(result$p.value, 0)
  expect_lte(result$p.value, 1)

  # Confidence interval should contain the estimate
  expect_length(result$conf.int, 2)
})

test_that("one-sample perm_t_test handles mu = 0 correctly", {
  set.seed(123)
  # Generate data with mean noticeably different from 0
  x_nonzero <- rnorm(20, mean = 3, sd = 1)
  result <- perm_t_test(x_nonzero, mu = 0, R = 199)

  expect_s3_class(result, "htest")
  # Should have small p-value since mean is clearly not 0
  expect_lt(result$p.value, 0.05)
})

# Two-Sample Tests-----------------

test_that("two-sample perm_t_test works with default settings", {
  set.seed(123)
  result <- perm_t_test(x_sample, y_sample, R = 199)

  expect_s3_class(result, "htest")
  expect_true(grepl("Two Sample", result$method))
  expect_true(grepl("Welch", result$method))  # Default is var.equal = FALSE
  expect_length(result$estimate, 2)
  expect_named(result$estimate, c("mean of x", "mean of y"))
})

test_that("two-sample perm_t_test works with var.equal = TRUE", {
  set.seed(123)
  result <- perm_t_test(x_sample, y_sample, var.equal = TRUE, R = 199)

  expect_s3_class(result, "htest")
  expect_true(grepl("Two Sample", result$method))
  expect_false(grepl("Welch", result$method))
})

test_that("two-sample perm_t_test works with formula interface", {
  set.seed(123)
  # Using built-in sleep data
  result <- perm_t_test(extra ~ group, data = sleep, R = 199)

  expect_s3_class(result, "htest")
  expect_equal(result$data.name, "extra by group")
})


# Paired Tests-----------------

test_that("paired perm_t_test works correctly", {
  set.seed(123)
  result <- perm_t_test(paired_x, paired_y, paired = TRUE, R = 199)

  expect_s3_class(result, "htest")
  expect_true(grepl("Paired", result$method))
  expect_named(result$estimate, "mean of the differences")
})

test_that("paired perm_t_test detects significant difference", {
  set.seed(123)
  # These paired samples have a consistent positive difference with some variability
  before <- c(5, 6, 7, 8, 9, 10, 11, 12)
  after <- c(5.8, 7.1, 7.9, 9.2, 9.8, 11.1, 11.9, 13.2)  # Mean diff ~1, with variability

  result <- perm_t_test(before, after, paired = TRUE, alternative = "less", R = 199)

  expect_s3_class(result, "htest")
  # before < after consistently, so "less" should be significant
  expect_lt(result$p.value, 0.05)
})


# Alternative Hypotheses -----------------
test_that("alternative = 'two.sided' works correctly",
{
  set.seed(123)
  result <- perm_t_test(x_sample, y_sample, alternative = "two.sided", R = 199)

  expect_equal(result$alternative, "two.sided")
  expect_length(result$conf.int, 2)
  expect_true(all(is.finite(result$conf.int)))
})

test_that("alternative = 'less' works correctly", {
  set.seed(123)
  result <- perm_t_test(x_sample, y_sample, alternative = "less", R = 199)

  expect_equal(result$alternative, "less")
  expect_equal(result$conf.int[1], -Inf)
  expect_true(is.finite(result$conf.int[2]))
})

test_that("alternative = 'greater' works correctly", {
  set.seed(123)
  result <- perm_t_test(x_sample, y_sample, alternative = "greater", R = 199)

  expect_equal(result$alternative, "greater")
  expect_true(is.finite(result$conf.int[1]))
  expect_equal(result$conf.int[2], Inf)
})

test_that("alternative = 'equivalence' works correctly", {
  set.seed(123)
  result <- perm_t_test(x_sample, y_sample,
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
  result <- perm_t_test(x_sample, y_sample,
                        alternative = "minimal.effect",
                        mu = c(-2, 2), R = 199)

  expect_equal(result$alternative, "minimal.effect")
  expect_length(result$null.value, 2)
  expect_equal(as.numeric(result$null.value), c(-2, 2))
})

test_that("equivalence with single mu value creates symmetric bounds", {
  set.seed(123)
  result <- perm_t_test(x_sample, y_sample,
                        alternative = "equivalence",
                        mu = 2, R = 199)

  expect_equal(as.numeric(result$null.value), c(-2, 2))
})


# Trimming (Yuen's Test) -----------------

test_that("trimmed t-test (tr > 0) works correctly", {
  set.seed(123)
  result <- perm_t_test(x_sample, y_sample, tr = 0.1, R = 199)

  expect_s3_class(result, "htest")
  expect_true(grepl("Yuen", result$method))
  expect_named(result$estimate, c("trimmed mean of x", "trimmed mean of y"))
})

test_that("one-sample trimmed t-test works", {
  set.seed(123)
  result <- perm_t_test(x_sample, mu = 5, tr = 0.2, R = 199)

  expect_true(grepl("Yuen", result$method))
  expect_named(result$estimate, "trimmed mean of x")
})

test_that("paired trimmed t-test works", {
  set.seed(123)
  result <- perm_t_test(paired_x, paired_y, paired = TRUE, tr = 0.1, R = 199)

  expect_true(grepl("Yuen", result$method))
  expect_true(grepl("Paired", result$method))
})


# P-value Method Tests -----------------

test_that("p_method = 'plusone' produces valid p-values", {
  set.seed(123)
  result <- perm_t_test(x_sample, y_sample, R = 199, p_method = "plusone")

  expect_gt(result$p.value, 0)  # Should never be exactly 0
  expect_lte(result$p.value, 1)
})

test_that("p_method = 'exact' produces valid p-values", {
  set.seed(123)
  result <- perm_t_test(x_sample, y_sample, R = 199, p_method = "exact")

  expect_gte(result$p.value, 0)
  expect_lte(result$p.value, 1)
})

test_that("different p_methods can produce different results", {
  set.seed(123)
  result_plus <- perm_t_test(x_sample, y_sample, R = 99, p_method = "plusone")
  set.seed(123)
  result_exact <- perm_t_test(x_sample, y_sample, R = 99, p_method = "exact")

  # They may or may not differ depending on the data, but both should be valid
  expect_gte(result_plus$p.value, 0)
  expect_gte(result_exact$p.value, 0)
})


# Studentized vs Non-Studentized Tests-----------------

test_that("perm_se = TRUE (studentized) produces valid results", {
  set.seed(123)
  result <- perm_t_test(x_sample, y_sample, perm_se = TRUE, R = 199)

  expect_s3_class(result, "htest")
  expect_gte(result$p.value, 0)
  expect_lte(result$p.value, 1)
})

test_that("perm_se = FALSE (non-studentized) produces valid results", {
  set.seed(123)
  result <- perm_t_test(x_sample, y_sample, perm_se = FALSE, R = 199)

  expect_s3_class(result, "htest")
  expect_gte(result$p.value, 0)
  expect_lte(result$p.value, 1)
})

test_that("studentized and non-studentized can produce different results", {
  set.seed(123)
  result_stud <- perm_t_test(x_sample, y_sample, perm_se = TRUE, R = 199)
  set.seed(123)
  result_non <- perm_t_test(x_sample, y_sample, perm_se = FALSE, R = 199)

  # t-statistics should be the same (calculated from original data)
  expect_equal(result_stud$statistic, result_non$statistic)

  # But p-values might differ (permutation distributions differ)
  # Both should be valid
  expect_gte(result_stud$p.value, 0)
  expect_gte(result_non$p.value, 0)
})


# Symmetric vs Equal-Tail Two-Sided Tests-----------------

test_that("symmetric = TRUE works for two-sided test", {
  set.seed(123)
  result <- perm_t_test(x_sample, y_sample,
                        alternative = "two.sided",
                        symmetric = TRUE, R = 199)

  expect_s3_class(result, "htest")
  expect_equal(result$alternative, "two.sided")
})

test_that("symmetric = FALSE works for two-sided test", {
  set.seed(123)
  result <- perm_t_test(x_sample, y_sample,
                        alternative = "two.sided",
                        symmetric = FALSE, R = 199)

  expect_s3_class(result, "htest")
  expect_equal(result$alternative, "two.sided")
})


# Exact Permutation Tests-----------------

test_that("exact permutation is computed for small samples", {
  set.seed(123)
  small_x <- c(1, 2, 3, 4)
  small_y <- c(5, 6, 7)

  # With R = NULL or large R, should compute exact permutations
  # choose(7, 4) = 35 permutations
  expect_message(
    result <- perm_t_test(small_x, small_y, R = NULL),
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
    result <- perm_t_test(small_x, mu = 0, R = NULL),
    "exact permutations"
  )

  expect_true(grepl("Exact", result$method))
  expect_equal(result$R.used, 2^5)
})

test_that("Randomization is used when R is specified and smaller than max perms", {
  set.seed(123)
  # Large enough sample that R = 199 triggers Randomization
  result <- perm_t_test(x_sample, y_sample, R = 199)

  expect_true(grepl("Randomization", result$method))
  expect_equal(result$R, 199)
})


# Missing Value Handling -----------------

test_that("perm_t_test handles NA values correctly", {
  set.seed(123)
  x_na <- c(x_sample, NA, NA)
  y_na <- c(NA, y_sample, NA)

  result <- perm_t_test(x_na, y_na, R = 199)

  expect_s3_class(result, "htest")
  # Should work with reduced sample
})
test_that("paired test handles NA values correctly", {
  set.seed(123)
  x_na <- c(paired_x, NA)
  y_na <- c(paired_y, NA)

  result <- perm_t_test(x_na, y_na, paired = TRUE, R = 199)

  expect_s3_class(result, "htest")
})


# Error Handling -----------------

test_that("error for invalid alpha", {
  expect_error(perm_t_test(x_sample, y_sample, alpha = -0.1),
               "'alpha' must be a single number between 0 and 1")
  expect_error(perm_t_test(x_sample, y_sample, alpha = 1.5),
               "'alpha' must be a single number between 0 and 1")
  expect_error(perm_t_test(x_sample, y_sample, alpha = c(0.05, 0.10)),
               "'alpha' must be a single number between 0 and 1")
})

test_that("error for invalid tr", {
  expect_error(perm_t_test(x_sample, y_sample, tr = -0.1),
               "'tr' must be a single number between 0 and 0.5")
  expect_error(perm_t_test(x_sample, y_sample, tr = 0.5),
               "'tr' must be a single number between 0 and 0.5")
  expect_error(perm_t_test(x_sample, y_sample, tr = 0.6),
               "'tr' must be a single number between 0 and 0.5")
})

test_that("error for invalid R", {
  expect_error(perm_t_test(x_sample, y_sample, R = 0),
               "'R' must be NULL .* or a positive integer")
  expect_error(perm_t_test(x_sample, y_sample, R = -10),
               "'R' must be NULL .* or a positive integer")
})

test_that("error for invalid perm_se", {
  expect_error(perm_t_test(x_sample, y_sample, perm_se = "yes"),
               "'perm_se' must be TRUE or FALSE")
  expect_error(perm_t_test(x_sample, y_sample, perm_se = 1),
               "'perm_se' must be TRUE or FALSE")
})

test_that("error for mu = 0 with equivalence/minimal.effect", {
  expect_error(perm_t_test(x_sample, y_sample, alternative = "equivalence", mu = 0),
               "For equivalence or minimal.effect testing")
})

test_that("error for wrong mu length with standard alternatives", {
  expect_error(perm_t_test(x_sample, y_sample, alternative = "two.sided", mu = c(1, 2)),
               "'mu' must be a single finite number for this alternative")
})

test_that("error for missing y in paired test", {
  expect_error(perm_t_test(x_sample, paired = TRUE),
               "'y' is missing for paired test")
})

test_that("error for insufficient sample size", {
  expect_error(perm_t_test(c(1), mu = 0),
               "not enough 'x' observations")
})

test_that("error for sample too small for trimming", {
  expect_error(perm_t_test(c(1, 2, 3), mu = 0, tr = 0.4),
               "Sample size too small for specified trimming")
})

test_that("error for formula with wrong number of groups", {
  df <- data.frame(value = 1:9, group = rep(c("a", "b", "c"), 3))
  expect_error(perm_t_test(value ~ group, data = df),
               "grouping factor must have exactly 2 levels")
})

test_that("error for incorrect formula", {
  expect_error(perm_t_test(~ group, data = sleep),
               "'formula' missing or incorrect")
})


# Reproducibility Tests -----------------

test_that("results are reproducible with set.seed", {
  set.seed(42)
  result1 <- perm_t_test(x_sample, y_sample, R = 199)
  set.seed(42)
  result2 <- perm_t_test(x_sample, y_sample, R = 199)

  expect_equal(result1$p.value, result2$p.value)
  expect_equal(result1$conf.int, result2$conf.int)
  expect_equal(result1$statistic, result2$statistic)
})


# Specific Value Tests (Regression Tests) -----------------

test_that("t-statistic matches expected calculation", {
  set.seed(123)
  x <- c(1, 2, 3, 4, 5)
  y <- c(6, 7, 8, 9, 10)

  result <- perm_t_test(x, y, R = 199)

  # Calculate expected t-statistic manually
  mx <- mean(x)
  my <- mean(y)
  vx <- var(x)
  vy <- var(y)
  stderr <- sqrt(vx/5 + vy/5)
  expected_t <- (mx - my) / stderr

  expect_equal(as.numeric(result$statistic), expected_t, tolerance = 1e-10)
})

test_that("one-sample t-statistic matches expected calculation", {
  x <- c(1, 2, 3, 4, 5)
  mu <- 2

  result <- perm_t_test(x, mu = mu, R = 199)

  # Calculate expected t-statistic manually
  mx <- mean(x)
  stderr <- sd(x) / sqrt(5)
  expected_t <- (mx - mu) / stderr

  expect_equal(as.numeric(result$statistic), expected_t, tolerance = 1e-10)
})

test_that("degrees of freedom are correct for two-sample Welch test", {
  x <- c(1, 2, 3, 4, 5)
  y <- c(6, 7, 8, 9, 10)

  result <- perm_t_test(x, y, var.equal = FALSE, R = 199)

  # Welch df calculation
  vx <- var(x)
  vy <- var(y)
  nx <- 5
  ny <- 5
  stderrx <- sqrt(vx/nx)
  stderry <- sqrt(vy/ny)
  stderr <- sqrt(stderrx^2 + stderry^2)
  expected_df <- stderr^4 / (stderrx^4/(nx-1) + stderry^4/(ny-1))

  expect_equal(as.numeric(result$parameter), expected_df, tolerance = 1e-10)
})

test_that("degrees of freedom are correct for two-sample pooled test", {
  x <- c(1, 2, 3, 4, 5)
  y <- c(6, 7, 8, 9, 10)

  result <- perm_t_test(x, y, var.equal = TRUE, R = 199)

  expected_df <- 5 + 5 - 2
  expect_equal(as.numeric(result$parameter), expected_df)
})


# Confidence Interval Tests-----------------

test_that("confidence intervals have correct coverage property conceptually", {
  set.seed(123)
  # Generate data where true difference is 0
  x_null <- rnorm(20, mean = 5, sd = 2)
  y_null <- rnorm(20, mean = 5, sd = 2)

  result <- perm_t_test(x_null, y_null, R = 499)

  # CI should contain 0 more often than not when null is true
  # This is a single test, so we just check the structure
  expect_length(result$conf.int, 2)
  expect_lt(result$conf.int[1], result$conf.int[2])
})

test_that("confidence level attribute is correct", {
  result_two <- perm_t_test(x_sample, y_sample, alternative = "two.sided",
                             alpha = 0.05, R = 199)
  expect_equal(attr(result_two$conf.int, "conf.level"), 0.95)

  result_less <- perm_t_test(x_sample, y_sample, alternative = "less",
                              alpha = 0.05, R = 199)
  expect_equal(attr(result_less$conf.int, "conf.level"), 0.95)

  result_equiv <- perm_t_test(x_sample, y_sample, alternative = "equivalence",
                               mu = c(-3, 3), alpha = 0.05, R = 199)
  expect_equal(attr(result_equiv$conf.int, "conf.level"), 0.90)
})


# Method String Tests-----------------

test_that("method string reflects test type correctly", {
  # Two-sample Welch
  result1 <- perm_t_test(x_sample, y_sample, var.equal = FALSE, R = 199)
  expect_match(result1$method, "Welch")
  expect_match(result1$method, "Two Sample")

  # Two-sample pooled
  result2 <- perm_t_test(x_sample, y_sample, var.equal = TRUE, R = 199)
  expect_false(grepl("Welch", result2$method))
  expect_match(result2$method, "Two Sample")

  # One-sample
  result3 <- perm_t_test(x_sample, mu = 5, R = 199)
  expect_match(result3$method, "One Sample")

  # Paired
  result4 <- perm_t_test(paired_x, paired_y, paired = TRUE, R = 199)
  expect_match(result4$method, "Paired")

  # Yuen (trimmed)
  result5 <- perm_t_test(x_sample, y_sample, tr = 0.1, R = 199)
  expect_match(result5$method, "Yuen")
})


# Permutation Distribution Tests-----------------

test_that("permutation distribution has correct length", {
  set.seed(123)
  R <- 199
  result <- perm_t_test(x_sample, y_sample, R = R, keep_perm = TRUE)

  expect_length(result$perm.stat, result$R.used)
  expect_length(result$perm.eff, result$R.used)
})

test_that("permutation distribution is numeric", {
  set.seed(123)
  result <- perm_t_test(x_sample, y_sample, R = 199, keep_perm = TRUE)

  expect_type(result$perm.stat, "double")
  expect_type(result$perm.eff, "double")
  expect_true(all(is.finite(result$perm.stat)))
  expect_true(all(is.finite(result$perm.eff)))
})


# Edge Cases -----------------

test_that("handles equal samples", {
  set.seed(123)
  x_eq <- c(1, 2, 3, 4, 5)

  result <- perm_t_test(x_eq, x_eq, R = 199)

  expect_s3_class(result, "htest")
  expect_equal(as.numeric(result$statistic), 0)
})

test_that("handles samples with equal means but different variances", {
  set.seed(123)
  x_same_mean <- c(4, 5, 6)
  y_same_mean <- c(2, 5, 8)  # Same mean, different variance

  result_welch <- perm_t_test(x_same_mean, y_same_mean, var.equal = FALSE, R = 199)
  result_pool <- perm_t_test(x_same_mean, y_same_mean, var.equal = TRUE, R = 199)

  # Both should have t-statistic of 0
  expect_equal(as.numeric(result_welch$statistic), 0)
  expect_equal(as.numeric(result_pool$statistic), 0)
})

test_that("handles very small samples for paired test", {
  set.seed(123)
  x_tiny <- c(1, 2, 3)
  y_tiny <- c(2.1, 2.8, 4.2)  # Add variability to differences

  # 2^3 = 8 exact permutations
  expect_message(
    result <- perm_t_test(x_tiny, y_tiny, paired = TRUE, R = NULL),
    "exact permutations"
  )

  expect_s3_class(result, "htest")
})


# Comparison with Standard t.test (Direction Consistency) -----------------

test_that("direction of effect matches t.test", {
  set.seed(123)
  # Clear difference: x < y
  x_low <- c(1, 2, 3, 4, 5)
  y_high <- c(10, 11, 12, 13, 14)

  perm_result <- perm_t_test(x_low, y_high, R = 199)
  t_result <- t.test(x_low, y_high)

  # Signs should match (compare numeric values without names)

  expect_equal(sign(as.numeric(perm_result$statistic)),
               sign(as.numeric(t_result$statistic)))

  # Estimates should be very close
  expect_equal(as.numeric(perm_result$estimate), as.numeric(t_result$estimate),
               tolerance = 1e-10)
})

test_that("one-sided tests give appropriate p-values", {
  set.seed(123)
  # x clearly less than y
  x_low <- c(1, 2, 3, 4, 5)
  y_high <- c(10, 11, 12, 13, 14)

  result_less <- perm_t_test(x_low, y_high, alternative = "less", R = 199)
  result_greater <- perm_t_test(x_low, y_high, alternative = "greater", R = 199)

  # "less" should be highly significant
  expect_lt(result_less$p.value, 0.05)
  # "greater" should not be significant
  expect_gt(result_greater$p.value, 0.5)
})


# Alpha Level Tests -----------------

test_that("different alpha levels produce appropriate confidence intervals", {
  set.seed(123)

  result_05 <- perm_t_test(x_sample, y_sample, alpha = 0.05, R = 299)
  result_01 <- perm_t_test(x_sample, y_sample, alpha = 0.01, R = 299)

  # 99% CI should be wider than 95% CI
  width_95 <- result_05$conf.int[2] - result_05$conf.int[1]
  width_99 <- result_01$conf.int[2] - result_01$conf.int[1]

  expect_gt(width_99, width_95)
})


# Print Method Test-----------------

test_that("print method works without error", {
  result <- perm_t_test(x_sample, y_sample, R = 99)

  # Should print without error
  expect_output(print(result))
})
