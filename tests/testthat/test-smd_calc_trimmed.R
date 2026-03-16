# Tests for trimmed SMD functionality in smd_calc and boot_smd_calc
# Following the specification in references/working/trimmed_smd_spec.md

# Helper to access internal functions
trim_rescale <- TOSTER:::trim_rescale
winvar <- TOSTER:::winvar
trim_h <- TOSTER:::trim_h

# --- 11.1 Parameter Validation ---

test_that("tr parameter validation works", {
  x <- rnorm(30)
  y <- rnorm(30)

  # Invalid tr values
  expect_error(smd_calc(x = x, y = y, tr = -0.1), "tr")
  expect_error(smd_calc(x = x, y = y, tr = 0.5), "tr")
  expect_error(smd_calc(x = x, y = y, tr = "a"), "tr")
  expect_error(smd_calc(x = x, y = y, tr = c(0.1, 0.2)), "tr")

  # tr = 0 should work (default)
  expect_no_error(smd_calc(x = x, y = y, tr = 0))

  # Incompatible combinations
  expect_error(smd_calc(x = x, y = y, paired = TRUE,
                        rm_correction = TRUE, tr = 0.2),
               "rm.*not.*supported.*trimming|not currently supported with trimming")
  expect_error(smd_calc(x = x, y = y, tr = 0.2, smd_ci = "goulet"),
               "goulet.*not.*supported.*trimming|not supported with trimming")

  # Same for boot
  expect_error(boot_smd_calc(x = x, y = y, tr = -0.1, R = 99), "tr")
  expect_error(boot_smd_calc(x = x, y = y, paired = TRUE,
                              rm_correction = TRUE, tr = 0.2, R = 99),
               "rm.*not.*supported.*trimming|not currently supported with trimming")
})


# --- 11.2 Backward Compatibility (tr = 0) ---

test_that("tr = 0 gives identical results to no-tr call", {
  set.seed(42)
  x <- rnorm(30, mean = 5, sd = 2)
  y <- rnorm(30, mean = 7, sd = 2)

  # Two-sample
  res_default <- smd_calc(x = x, y = y)
  res_tr0 <- smd_calc(x = x, y = y, tr = 0)
  expect_equal(res_default$estimate, res_tr0$estimate)
  expect_equal(res_default$conf.int, res_tr0$conf.int)
  expect_equal(res_default$stderr, res_tr0$stderr)

  # Paired
  res_default_p <- smd_calc(x = x, y = y, paired = TRUE)
  res_tr0_p <- smd_calc(x = x, y = y, paired = TRUE, tr = 0)
  expect_equal(res_default_p$estimate, res_tr0_p$estimate)
  expect_equal(res_default_p$conf.int, res_tr0_p$conf.int)

  # One-sample
  res_default_o <- smd_calc(x = x)
  res_tr0_o <- smd_calc(x = x, tr = 0)
  expect_equal(res_default_o$estimate, res_tr0_o$estimate)
  expect_equal(res_default_o$conf.int, res_tr0_o$conf.int)
})


# --- 11.3 Known Values Under Normality ---

test_that("trimmed SMD approximates untrimmed under normality for large n", {
  set.seed(123)
  x <- rnorm(200, mean = 0, sd = 1)
  y <- rnorm(200, mean = 0.5, sd = 1)

  res_std <- smd_calc(x = x, y = y, bias_correction = FALSE, smd_ci = "z")
  res_trim <- smd_calc(x = x, y = y, bias_correction = FALSE, tr = 0.2, smd_ci = "z")

  # Under normality with large n, these should be in the same ballpark
  expect_equal(unname(res_std$estimate), unname(res_trim$estimate), tolerance = 0.15)
})


# --- 11.4 Rescaling Constant Accuracy ---

test_that("trim_rescale returns correct known values", {
  # c(0) = 1
  expect_equal(trim_rescale(0), 1)

  # c(0.2) should be approximately 0.642
  expect_equal(trim_rescale(0.2), 0.642, tolerance = 0.001)

  # c(0.1) -- compute expected value
  a <- qnorm(0.9)
  expected <- sqrt(1 - 0.2 + 0.2 * a^2 - 2 * a * dnorm(a))
  expect_equal(trim_rescale(0.1), expected)
})


# --- 11.5 Robustness to Outliers ---

test_that("trimmed SMD is robust to outliers", {
  set.seed(456)
  x <- rnorm(30, mean = 0, sd = 1)
  y <- rnorm(30, mean = 0.5, sd = 1)

  # Add outliers
  x_out <- c(x, 50, -50)  # extreme outliers
  y_out <- c(y, 55, -45)

  res_clean <- smd_calc(x = x, y = y, bias_correction = FALSE, smd_ci = "z")
  res_contaminated <- smd_calc(x = x_out, y = y_out, bias_correction = FALSE, smd_ci = "z")
  res_robust <- smd_calc(x = x_out, y = y_out, bias_correction = FALSE, tr = 0.2, smd_ci = "z")

  # Standard SMD should be heavily affected by outliers
  expect_true(abs(unname(res_contaminated$estimate) - unname(res_clean$estimate)) > 0.1)

  # Trimmed SMD should be much closer to the clean data result
  expect_true(abs(unname(res_robust$estimate) - unname(res_clean$estimate)) <
              abs(unname(res_contaminated$estimate) - unname(res_clean$estimate)))
})


# --- 11.6 Degrees of Freedom Modification ---

test_that("degrees of freedom are correctly adjusted for trimming", {
  set.seed(789)
  x <- rnorm(30)
  y <- rnorm(30)

  # For tr = 0.2, g = floor(0.2 * 30) = 6, h = 30 - 12 = 18
  # df for pooled two-sample = h1 + h2 - 2 = 18 + 18 - 2 = 34
  res <- smd_calc(x = x, y = y, tr = 0.2, smd_ci = "t", var.equal = TRUE,
                  alternative = "two.sided", null.value = 0, test_method = "t")
  # The df should be 34 (not 58 as in the untrimmed case)
  expect_equal(res$parameter[["df"]], 34)

  # For one-sample: df = h - 1 = 18 - 1 = 17
  res_one <- smd_calc(x = x, tr = 0.2, smd_ci = "t",
                      alternative = "two.sided", null.value = 0, test_method = "t")
  expect_equal(res_one$parameter[["df"]], 17)
})


# --- 11.7 All Denominator Options with Trimming ---

test_that("all supported denom options work with tr > 0", {
  set.seed(101)
  x <- rnorm(30, mean = 0, sd = 1)
  y <- rnorm(30, mean = 0.8, sd = 1.5)

  # Two-sample denominators
  expect_no_error(smd_calc(x = x, y = y, denom = "pooled", tr = 0.1))
  expect_no_error(smd_calc(x = x, y = y, denom = "avg", tr = 0.1))
  expect_no_error(smd_calc(x = x, y = y, denom = "glass1", tr = 0.1))
  expect_no_error(smd_calc(x = x, y = y, denom = "glass2", tr = 0.1))

  # Paired denominators
  expect_no_error(smd_calc(x = x, y = y, paired = TRUE, denom = "z", tr = 0.1))
  expect_no_error(smd_calc(x = x, y = y, paired = TRUE, denom = "glass1", tr = 0.1))

  # One-sample
  expect_no_error(smd_calc(x = x, tr = 0.1))
})


# --- 11.8 Bootstrap Pass-Through ---

test_that("boot_smd_calc correctly passes tr to internal calls", {
  set.seed(202)
  x <- rnorm(20, mean = 0, sd = 1)
  y <- rnorm(20, mean = 0.5, sd = 1)

  # Should run without error
  res <- boot_smd_calc(x = x, y = y, tr = 0.2, R = 199, boot_ci = "perc")
  expect_s3_class(res, "htest")
  expect_true(!is.na(res$estimate))
  expect_true(all(!is.na(res$conf.int)))

  # Estimate should match smd_calc with same trimming
  res_point <- smd_calc(x = x, y = y, tr = 0.2, smd_ci = "z")
  expect_equal(unname(res$estimate), unname(res_point$estimate))
})


# --- 11.9 CI Method Compatibility ---

test_that("CI methods work with trimming", {
  set.seed(303)
  x <- rnorm(40)
  y <- rnorm(40, mean = 0.5)

  # These should work
  expect_no_error(smd_calc(x = x, y = y, tr = 0.2, smd_ci = "nct"))
  expect_no_error(smd_calc(x = x, y = y, tr = 0.2, smd_ci = "t"))
  expect_no_error(smd_calc(x = x, y = y, tr = 0.2, smd_ci = "z"))

  # This should fail
  expect_error(smd_calc(x = x, y = y, tr = 0.2, smd_ci = "goulet"),
               "goulet.*not.*supported|not supported with trimming")
})


# --- 11.10 Hypothesis Testing with Trimming ---

test_that("hypothesis testing works with trimmed SMD", {
  set.seed(404)
  x <- rnorm(30)
  y <- rnorm(30, mean = 1)

  # Two-sided test
  res <- smd_calc(x = x, y = y, tr = 0.2,
                  alternative = "two.sided", null.value = 0)
  expect_true(!is.na(res$p.value))
  expect_true(res$p.value >= 0 && res$p.value <= 1)

  # Equivalence test
  res_eq <- smd_calc(x = x, y = y, tr = 0.2,
                     alternative = "equivalence", null.value = c(-2, 2))
  expect_true(!is.na(res_eq$p.value))

  # Bootstrap hypothesis test
  res_boot <- boot_smd_calc(x = x, y = y, tr = 0.2, R = 199,
                            alternative = "two.sided", null.value = 0)
  expect_true(!is.na(res_boot$p.value))
})


# --- 11.11 Edge Case: Small Sample with Trimming ---

test_that("trimming with small samples is handled correctly", {
  x <- rnorm(5)
  y <- rnorm(5)

  # tr = 0.2 with n = 5: g = floor(0.2*5) = 1, h = 3, which is OK
  expect_no_error(smd_calc(x = x, y = y, tr = 0.2))

  # tr = 0.4 with n = 5: g = floor(0.4*5) = 2, h = 1, should fail
  expect_error(smd_calc(x = x, y = y, tr = 0.4),
               "too many observations")
})


# --- 11.12 Formula Interface with Trimming ---

test_that("formula interface passes tr correctly", {
  set.seed(505)
  df <- data.frame(
    value = c(rnorm(25), rnorm(25, mean = 0.5)),
    group = factor(rep(c("A", "B"), each = 25))
  )

  res_formula <- smd_calc(value ~ group, data = df, tr = 0.2)
  res_direct <- smd_calc(x = df$value[df$group == "A"],
                         y = df$value[df$group == "B"],
                         tr = 0.2)
  expect_equal(unname(res_formula$estimate), unname(res_direct$estimate))
  expect_equal(as.numeric(res_formula$conf.int), as.numeric(res_direct$conf.int))
})


# --- Additional: Internal helper functions ---

test_that("winvar returns correct values", {
  x <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)

  # With tr = 0, should equal var()
  expect_equal(winvar(x, tr = 0), var(x))

  # With tr = 0.2, g = 2, Winsorize: 3,3,3,4,5,6,7,8,8,8
  y <- c(3, 3, 3, 4, 5, 6, 7, 8, 8, 8)
  expect_equal(winvar(x, tr = 0.2), var(y))
})

test_that("trim_h returns correct effective sample size", {
  expect_equal(trim_h(30, 0), 30)
  expect_equal(trim_h(30, 0.2), 18)  # g = 6, h = 30 - 12 = 18
  expect_equal(trim_h(10, 0.1), 8)   # g = 1, h = 10 - 2 = 8
})


# --- Additional: SMD label includes trim note ---

test_that("SMD label includes trimming note when tr > 0", {
  set.seed(606)
  x <- rnorm(30)
  y <- rnorm(30, mean = 0.5)

  res_trim <- smd_calc(x = x, y = y, tr = 0.2, bias_correction = FALSE)
  expect_true(grepl("trimmed", res_trim$method))

  res_notrim <- smd_calc(x = x, y = y, tr = 0, bias_correction = FALSE)
  expect_false(grepl("trimmed", res_notrim$method))
})
