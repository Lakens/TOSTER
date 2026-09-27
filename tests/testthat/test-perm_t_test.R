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
  skip_on_cran()
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
  skip_on_cran()
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
  skip_on_cran()
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
  skip_on_cran()
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
  skip_on_cran()
  set.seed(123)
  result <- perm_t_test(x_sample, y_sample, R = 199)

  expect_s3_class(result, "htest")
  expect_true(grepl("Two Sample", result$method))
  expect_true(grepl("Welch", result$method))  # Default is var.equal = FALSE
  expect_length(result$estimate, 3)
  expect_named(result$estimate, c("mean of group x", "mean of group y", "mean difference (x - y)"))
})

test_that("two-sample perm_t_test works with var.equal = TRUE", {
  skip_on_cran()
  set.seed(123)
  result <- perm_t_test(x_sample, y_sample, var.equal = TRUE, R = 199)

  expect_s3_class(result, "htest")
  expect_true(grepl("Two Sample", result$method))
  expect_false(grepl("Welch", result$method))
})

test_that("two-sample perm_t_test works with formula interface", {
  skip_on_cran()
  set.seed(123)
  # Using built-in sleep data
  result <- perm_t_test(extra ~ group, data = sleep, R = 199)

  expect_s3_class(result, "htest")
  expect_equal(result$data.name, "extra by group")
})


# Paired Tests-----------------

test_that("paired perm_t_test works correctly", {
  skip_on_cran()
  set.seed(123)
  result <- perm_t_test(paired_x, paired_y, paired = TRUE, R = 199)

  expect_s3_class(result, "htest")
  expect_true(grepl("Paired", result$method))
  expect_named(result$estimate, "mean of the differences (z = x - y)")
})

test_that("paired perm_t_test detects significant difference", {
  skip_on_cran()
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
  skip_on_cran()
  set.seed(123)
  result <- perm_t_test(x_sample, y_sample, alternative = "two.sided", R = 199)

  expect_equal(result$alternative, "two.sided")
  expect_length(result$conf.int, 2)
  expect_true(all(is.finite(result$conf.int)))
})

test_that("alternative = 'less' works correctly", {
  skip_on_cran()
  set.seed(123)
  result <- perm_t_test(x_sample, y_sample, alternative = "less", R = 199)

  expect_equal(result$alternative, "less")
  expect_equal(result$conf.int[1], -Inf)
  expect_true(is.finite(result$conf.int[2]))
})

test_that("alternative = 'greater' works correctly", {
  skip_on_cran()
  set.seed(123)
  result <- perm_t_test(x_sample, y_sample, alternative = "greater", R = 199)

  expect_equal(result$alternative, "greater")
  expect_true(is.finite(result$conf.int[1]))
  expect_equal(result$conf.int[2], Inf)
})

test_that("alternative = 'equivalence' works correctly", {
  skip_on_cran()
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
  skip_on_cran()
  set.seed(123)
  result <- perm_t_test(x_sample, y_sample,
                        alternative = "minimal.effect",
                        mu = c(-2, 2), R = 199)

  expect_equal(result$alternative, "minimal.effect")
  expect_length(result$null.value, 2)
  expect_equal(as.numeric(result$null.value), c(-2, 2))
})

test_that("equivalence with single mu value creates symmetric bounds", {
  skip_on_cran()
  set.seed(123)
  result <- perm_t_test(x_sample, y_sample,
                        alternative = "equivalence",
                        mu = 2, R = 199)

  expect_equal(as.numeric(result$null.value), c(-2, 2))
})


# Trimming (Yuen's Test) -----------------

test_that("trimmed t-test (tr > 0) works correctly", {
  skip_on_cran()
  set.seed(123)
  result <- perm_t_test(x_sample, y_sample, tr = 0.1, R = 199)

  expect_s3_class(result, "htest")
  expect_true(grepl("Yuen", result$method))
  expect_equal(length(result$estimate), 3)
  expect_equal(names(result$estimate)[1:2],
               c("trimmed mean of x", "trimmed mean of y"))
  expect_true(grepl("trimmed mean difference", names(result$estimate)[3]))
})

test_that("one-sample trimmed t-test works", {
  skip_on_cran()
  set.seed(123)
  result <- perm_t_test(x_sample, mu = 5, tr = 0.2, R = 199)

  expect_true(grepl("Yuen", result$method))
  expect_true(grepl("trimmed mean of x", names(result$estimate)))
})

test_that("paired trimmed t-test works", {
  skip_on_cran()
  set.seed(123)
  result <- perm_t_test(paired_x, paired_y, paired = TRUE, tr = 0.1, R = 199)

  expect_true(grepl("Yuen", result$method))
  expect_true(grepl("Paired", result$method))
})


# P-value Method Tests -----------------

test_that("p_method = 'plusone' produces valid p-values", {
  skip_on_cran()
  set.seed(123)
  result <- perm_t_test(x_sample, y_sample, R = 199, p_method = "plusone")

  expect_gt(result$p.value, 0)  # Should never be exactly 0
  expect_lte(result$p.value, 1)
})

test_that("p_method = 'exact' produces valid p-values", {
  skip_on_cran()
  set.seed(123)
  result <- perm_t_test(x_sample, y_sample, R = 199, p_method = "exact")

  expect_gte(result$p.value, 0)
  expect_lte(result$p.value, 1)
})

test_that("different p_methods can produce different results", {
  skip_on_cran()
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
  skip_on_cran()
  set.seed(123)
  result <- perm_t_test(x_sample, y_sample, perm_se = TRUE, R = 199)

  expect_s3_class(result, "htest")
  expect_gte(result$p.value, 0)
  expect_lte(result$p.value, 1)
})

test_that("perm_se = FALSE (non-studentized) produces valid results", {
  skip_on_cran()
  set.seed(123)
  result <- perm_t_test(x_sample, y_sample, perm_se = FALSE, R = 199)

  expect_s3_class(result, "htest")
  expect_gte(result$p.value, 0)
  expect_lte(result$p.value, 1)
})

test_that("studentized and non-studentized can produce different results", {
  skip_on_cran()
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
  skip_on_cran()
  set.seed(123)
  result <- perm_t_test(x_sample, y_sample,
                        alternative = "two.sided",
                        symmetric = TRUE, R = 199)

  expect_s3_class(result, "htest")
  expect_equal(result$alternative, "two.sided")
})

test_that("symmetric = FALSE works for two-sided test", {
  skip_on_cran()
  set.seed(123)
  result <- perm_t_test(x_sample, y_sample,
                        alternative = "two.sided",
                        symmetric = FALSE, R = 199)

  expect_s3_class(result, "htest")
  expect_equal(result$alternative, "two.sided")
})


# Exact Permutation Tests-----------------

test_that("exact permutation is computed for small samples", {
  skip_on_cran()
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
  skip_on_cran()
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
  skip_on_cran()
  set.seed(123)
  # Large enough sample that R = 199 triggers Randomization
  result <- perm_t_test(x_sample, y_sample, R = 199)

  expect_true(grepl("Randomization", result$method))
  expect_equal(result$R, 199)
})


# Missing Value Handling -----------------

test_that("perm_t_test handles NA values correctly", {
  skip_on_cran()
  set.seed(123)
  x_na <- c(x_sample, NA, NA)
  y_na <- c(NA, y_sample, NA)

  result <- perm_t_test(x_na, y_na, R = 199)

  expect_s3_class(result, "htest")
  # Should work with reduced sample
})
test_that("paired test handles NA values correctly", {
  skip_on_cran()
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
  skip_on_cran()
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
  skip_on_cran()
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
  skip_on_cran()
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
  skip_on_cran()
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
  skip_on_cran()
  x <- c(1, 2, 3, 4, 5)
  y <- c(6, 7, 8, 9, 10)

  result <- perm_t_test(x, y, var.equal = TRUE, R = 199)

  expected_df <- 5 + 5 - 2
  expect_equal(as.numeric(result$parameter), expected_df)
})


# Confidence Interval Tests-----------------

test_that("confidence intervals have correct coverage property conceptually", {
  skip_on_cran()
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
  skip_on_cran()
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
  skip_on_cran()
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
  skip_on_cran()
  set.seed(123)
  R <- 199
  result <- perm_t_test(x_sample, y_sample, R = R, keep_perm = TRUE)

  expect_length(result$perm.stat, result$R.used)
  expect_length(result$perm.eff, result$R.used)
})

test_that("permutation distribution is numeric", {
  skip_on_cran()
  set.seed(123)
  result <- perm_t_test(x_sample, y_sample, R = 199, keep_perm = TRUE)

  expect_type(result$perm.stat, "double")
  expect_type(result$perm.eff, "double")
  expect_true(all(is.finite(result$perm.stat)))
  expect_true(all(is.finite(result$perm.eff)))
})


# Edge Cases -----------------

test_that("handles equal samples", {
  skip_on_cran()
  set.seed(123)
  x_eq <- c(1, 2, 3, 4, 5)

  result <- perm_t_test(x_eq, x_eq, R = 199)

  expect_s3_class(result, "htest")
  expect_equal(as.numeric(result$statistic), 0)
})

test_that("handles samples with equal means but different variances", {
  skip_on_cran()
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
  skip_on_cran()
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
  skip_on_cran()
  set.seed(123)
  # Clear difference: x < y
  x_low <- c(1, 2, 3, 4, 5)
  y_high <- c(10, 11, 12, 13, 14)

  perm_result <- perm_t_test(x_low, y_high, R = 199)
  t_result <- simple_htest(x_low, y_high, test = "t")

  # Signs should match (compare numeric values without names)

  expect_equal(sign(as.numeric(perm_result$statistic)),
               sign(as.numeric(t_result$statistic)))

  # Estimates should be very close
  expect_equal(as.numeric(perm_result$estimate), as.numeric(t_result$estimate),
               tolerance = 1e-10)
})

test_that("one-sided tests give appropriate p-values", {
  skip_on_cran()
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
  skip_on_cran()
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
  skip_on_cran()
  result <- perm_t_test(x_sample, y_sample, R = 99)

  # Should print without error
  expect_output(print(result))
})


# CI / p-value duality (#120) -----------------

# Re-run the test at a new null value using the same permutations
# (same seed for randomized R, identical enumeration for exact)
perm_p_at <- function(args, mu, seed = 1) {
  args$mu <- mu
  set.seed(seed)
  suppressMessages(do.call(perm_t_test, args))$p.value
}

# Check that mu just inside the CI gives p > alpha and just outside gives p <= alpha
expect_ci_duality <- function(args, seed = 1) {
  alpha <- if (is.null(args$alpha)) 0.05 else args$alpha
  set.seed(seed)
  res <- suppressMessages(do.call(perm_t_test, args))
  ci <- res$conf.int
  eps <- 1e-6 * res$stderr

  if (is.finite(ci[1])) {
    expect_gt(perm_p_at(args, ci[1] + eps, seed), alpha)
    expect_lte(perm_p_at(args, ci[1] - eps, seed), alpha)
  }
  if (is.finite(ci[2])) {
    expect_gt(perm_p_at(args, ci[2] - eps, seed), alpha)
    expect_lte(perm_p_at(args, ci[2] + eps, seed), alpha)
  }
  invisible(res)
}

set.seed(2024)
bf_x <- rnorm(8, mean = 1, sd = 4)
bf_y <- rnorm(20, mean = 0, sd = 1)
small_x <- rnorm(5, mean = 1, sd = 3)
small_y <- rnorm(9, mean = 0, sd = 1)

test_that("CI agrees with p-value: two-sample Welch (Behrens-Fisher), randomized", {
  skip_on_cran()
  for (alt in c("two.sided", "less", "greater")) {
    for (sym in c(TRUE, FALSE)) {
      expect_ci_duality(list(x = bf_x, y = bf_y, alternative = alt,
                             symmetric = sym, R = 999))
    }
  }
})

test_that("CI agrees with p-value: two-sample exact enumeration, both p_methods", {
  skip_on_cran()
  for (pm in c("exact", "plusone")) {
    for (alt in c("two.sided", "less", "greater")) {
      expect_ci_duality(list(x = small_x, y = small_y, alternative = alt,
                             p_method = pm))
      expect_ci_duality(list(x = small_x, y = small_y, alternative = alt,
                             symmetric = FALSE, p_method = pm))
    }
  }
})

test_that("CI agrees with p-value: var.equal, trimming, and perm_se = FALSE", {
  skip_on_cran()
  for (alt in c("two.sided", "less", "greater")) {
    expect_ci_duality(list(x = bf_x, y = bf_y, alternative = alt,
                           var.equal = TRUE, R = 499))
    expect_ci_duality(list(x = bf_x, y = bf_y, alternative = alt,
                           tr = 0.2, R = 499))
    expect_ci_duality(list(x = bf_x, y = bf_y, alternative = alt,
                           perm_se = FALSE, R = 499))
    expect_ci_duality(list(x = bf_x, y = bf_y, alternative = alt,
                           var.equal = TRUE, tr = 0.2, R = 499,
                           p_method = "exact"))
  }
})

test_that("CI agrees with p-value: one-sample and paired", {
  skip_on_cran()
  set.seed(7)
  x1 <- rexp(10) - 0.3
  for (alt in c("two.sided", "less", "greater")) {
    for (sym in c(TRUE, FALSE)) {
      # exact sign-flip enumeration (1024 permutations)
      expect_ci_duality(list(x = x1, alternative = alt, symmetric = sym))
      # randomized sign flips with trimming
      expect_ci_duality(list(x = c(x1, x1 * 1.5), alternative = alt,
                             symmetric = sym, tr = 0.1, R = 999))
      # paired
      expect_ci_duality(list(x = paired_x, y = paired_y, paired = TRUE,
                             alternative = alt, symmetric = sym))
    }
  }
})

test_that("CI agrees with p-value for alpha other than 0.05", {
  skip_on_cran()
  expect_ci_duality(list(x = bf_x, y = bf_y, alpha = 0.1, R = 999))
  expect_ci_duality(list(x = bf_x, y = bf_y, alpha = 0.01, R = 999,
                         symmetric = FALSE))
})

test_that("equivalence and minimal effect decisions agree with the 1 - 2*alpha CI", {
  skip_on_cran()
  base <- list(x = bf_x, y = bf_y, R = 999)
  set.seed(1)
  ci <- suppressMessages(do.call(perm_t_test, c(base, list(
    alternative = "equivalence", mu = c(-10, 10)))))$conf.int
  expect_equal(attr(ci, "conf.level"), 0.90)
  eps <- 1e-6

  run <- function(alt, mu) {
    set.seed(1)
    suppressMessages(do.call(perm_t_test, c(base, list(
      alternative = alt, mu = mu))))
  }

  # CI is identical regardless of the bounds
  expect_equal(run("equivalence", c(-1, 1))$conf.int, ci)
  expect_equal(run("minimal.effect", c(-1, 1))$conf.int, ci)

  # Equivalence: significant iff CI lies within the bounds
  expect_lte(run("equivalence", c(ci[1] - eps, ci[2] + eps))$p.value, 0.05)
  expect_gt(run("equivalence", c(ci[1] + eps, ci[2] + eps))$p.value, 0.05)
  expect_gt(run("equivalence", c(ci[1] - eps, ci[2] - eps))$p.value, 0.05)

  # Minimal effect: significant iff CI lies entirely outside the bounds
  expect_lte(run("minimal.effect", c(ci[2] + eps, ci[2] + 5))$p.value, 0.05)
  expect_gt(run("minimal.effect", c(ci[2] - eps, ci[2] + 5))$p.value, 0.05)
  expect_lte(run("minimal.effect", c(ci[1] - 5, ci[1] - eps))$p.value, 0.05)
  expect_gt(run("minimal.effect", c(ci[1] - 5, ci[1] + eps))$p.value, 0.05)
})

test_that("symmetric two-sided CI is centered on the estimate", {
  skip_on_cran()
  set.seed(1)
  res <- suppressMessages(perm_t_test(bf_x, bf_y, R = 999, symmetric = TRUE))
  est <- unname(res$estimate[3])
  expect_equal(est - res$conf.int[1], res$conf.int[2] - est)
})

test_that("CI is infinite when too few permutations to reject", {
  # one-sample n = 4: 16 sign flips; with plusone the smallest one-sided
  # p-value is 1/17 > 0.05 and the equal-tail p-value is at least 2/17
  x4 <- c(1.2, 0.8, 2.1, 1.5)
  res <- suppressMessages(perm_t_test(x4, symmetric = FALSE, p_method = "plusone"))
  expect_equal(res$R.used, 16)
  expect_equal(res$conf.int[1], -Inf)
  expect_equal(res$conf.int[2], Inf)

  res_less <- suppressMessages(perm_t_test(x4, alternative = "less",
                                           p_method = "plusone"))
  expect_equal(res_less$conf.int[2], Inf)

  # with exact (b/R) counting, p = 0 is attainable, so the CI is finite
  res_exact <- suppressMessages(perm_t_test(x4, symmetric = FALSE, p_method = "exact"))
  expect_true(all(is.finite(res_exact$conf.int)))
})

test_that("perm_crit returns the order statistic matching the p-value rule", {
  tstat <- c(-2, -1, 0, 1, 2, 3, 4, 5, 6, 7)
  # exact: reject when b/10 <= 0.2, i.e. b <= 2
  # upper: 3rd largest; lower: 3rd smallest; abs: 3rd largest of |T|
  # (critical values are widened by perm_tol() to match perm_count())
  expect_equal(perm_crit(tstat, 0.2, "exact", "upper"), 5 + perm_tol(5))
  expect_equal(perm_crit(tstat, 0.2, "exact", "lower"), 0 - perm_tol(0))
  expect_equal(perm_crit(tstat, 0.2, "exact", "abs"), 5 + perm_tol(5))
  # plusone: reject when (b+1)/11 <= 0.2, i.e. b <= 1
  expect_equal(perm_crit(tstat, 0.2, "plusone", "upper"), 6 + perm_tol(6))
  expect_equal(perm_crit(tstat, 0.2, "plusone", "lower"), -1 - perm_tol(-1))
  # exact: b = 0 gives p = 0, so the extreme order statistic is the critical value
  expect_equal(perm_crit(tstat, 0.05, "exact", "upper"), 7 + perm_tol(7))
  expect_equal(perm_crit(tstat, 0.05, "exact", "lower"), -2 - perm_tol(-2))
  # plusone: (0+1)/11 > 0.05, so the test can never reject
  expect_equal(perm_crit(tstat, 0.05, "plusone", "upper"), Inf)
  expect_equal(perm_crit(tstat, 0.05, "plusone", "lower"), -Inf)
  expect_equal(perm_crit(tstat, 0.05, "plusone", "abs"), Inf)
})

# Floating point ties in permutation counts ----

test_that("perm_count treats values within floating point error as ties", {
  t_obs <- 0.1 + 0.2          # 0.30000000000000004
  TSTAT <- c(-0.1, 0.3, 0.3, 2) # 0.3 is mathematically equal but slightly smaller
  expect_false(all(TSTAT[2:3] >= t_obs))
  expect_equal(perm_count(TSTAT, t_obs, "ge"), 3)
  expect_equal(perm_count(TSTAT, t_obs, "le"), 3)
  expect_equal(perm_count(TSTAT, -t_obs, "abs"), 3)
  # genuinely different values are unaffected
  expect_equal(perm_count(c(0.2999, 0.3001), 0.3, "ge"), 1)
  expect_equal(perm_count(c(0.2999, 0.3001), 0.3, "le"), 1)
})

test_that("tied data: every mathematically tied permutation is counted", {
  skip_on_cran()
  # Ordinal data give many permutations whose statistics tie exactly in
  # theory but differ in the last bits in floating point
  x <- c(3, 5, 3, 3, 4)
  y <- c(4, 3, 3, 2, 4)
  res <- suppressMessages(perm_t_test(x, y, p_method = "exact"))
  near <- abs(res$perm.stat - res$statistic) < 1e-9
  expect_gt(sum(near), 1)
  b_expected <- sum(abs(res$perm.stat) >= abs(res$statistic) - 1e-9)
  expect_equal(res$p.value, b_expected / res$R.used)
})

test_that("CI agrees with p-value under Behrens-Fisher (regression for #120)", {
  skip_on_cran()
  # Small group with large variance: the old percentile CI of raw permuted
  # differences disagreed with the studentized p-value for these data
  set.seed(1)
  x <- rnorm(8, 0.9, 4)
  y <- rnorm(30, 0, 1)
  set.seed(1)
  res <- suppressMessages(perm_t_test(x, y, R = 999))
  old_ci <- quantile(res$perm.eff, c(0.025, 0.975), names = FALSE)

  p_sig <- res$p.value <= 0.05
  expect_false(p_sig)
  # old interval excluded 0 despite p > 0.05
  expect_true(old_ci[1] > 0 || old_ci[2] < 0)
  # new interval agrees with the p-value
  expect_equal(res$conf.int[1] > 0 || res$conf.int[2] < 0, p_sig)
})
