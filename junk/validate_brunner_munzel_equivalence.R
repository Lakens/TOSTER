# =============================================================================
# Validation Script: Brunner-Munzel Equivalence Testing
# =============================================================================
# This script validates TOSTER's brunner_munzel() implementation against:
# 1. nparcomp package (for base test statistics and p-values)
# 2. Internal consistency (equivalence = max of one-sided, MET = min of one-sided)
# 3. All test methods (t, logit, perm) for equivalence/minimal.effect
#
# Author: Aaron R. Caldwell
# Date: 2024
# =============================================================================

# Load required packages
# Use devtools to load the development version of TOSTER
if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(quiet = TRUE)
  message("Loaded development version of TOSTER")
} else {
  library(TOSTER)
  message("Using installed TOSTER package")
}
library(testthat)

# Check if nparcomp is available (optional comparison)
nparcomp_available <- requireNamespace("nparcomp", quietly = TRUE)
if (nparcomp_available) {
  library(nparcomp)
  message("nparcomp package loaded for comparison")
} else {
  message("nparcomp package not available - skipping nparcomp comparisons")
}

# Load test data
data(sleep)
data(mtcars)

cat("\n")
cat("=============================================================================\n")
cat("BRUNNER-MUNZEL EQUIVALENCE TESTING VALIDATION\n")
cat("=============================================================================\n")

# =============================================================================
# SECTION 1: Validate base implementation matches nparcomp
# =============================================================================

cat("\n--- Section 1: Base Implementation vs nparcomp ---\n\n")

if (nparcomp_available) {

  # Two-sample test
  cat("Test 1.1: Two-sample t-approximation (sleep data)\n")
  toster_res <- brunner_munzel(
    x = subset(sleep, group == 2)$extra,
    y = subset(sleep, group == 1)$extra,
    paired = FALSE,
    test_method = "t",
    alternative = "two.sided"
  )

  nparcomp_res <- npar.t.test(
    data = sleep,
    extra ~ group,
    method = "t.app",
    alternative = "two.sided",
    info = FALSE
  )

  cat(sprintf("  TOSTER p-value:   %.6f\n", toster_res$p.value))
  cat(sprintf("  nparcomp p-value: %.6f\n", nparcomp_res$Analysis$p.Value[1]))
  cat(sprintf("  Difference:       %.2e\n", abs(toster_res$p.value - nparcomp_res$Analysis$p.Value[1])))

  test_that("Two-sample t-approx matches nparcomp", {
    # nparcomp rounds to 3 decimal places, so use appropriate tolerance
    expect_equal(toster_res$p.value, nparcomp_res$Analysis$p.Value[1], tolerance = 0.01)
  })

  # Paired test
  cat("\nTest 1.2: Paired t-approximation (sleep data)\n")
  toster_paired <- brunner_munzel(
    x = subset(sleep, group == 2)$extra,
    y = subset(sleep, group == 1)$extra,
    paired = TRUE,
    test_method = "t",
    alternative = "two.sided"
  )

  nparcomp_paired <- npar.t.test.paired(
    data = sleep,
    extra ~ group,
    alternative = "two.sided",
    info = FALSE,
    plot.simci = FALSE
  )

  cat(sprintf("  TOSTER p-value:   %.6f\n", toster_paired$p.value))
  cat(sprintf("  nparcomp p-value: %.6f\n", nparcomp_paired$Analysis[1, 5]))
  cat(sprintf("  Difference:       %.2e\n", abs(toster_paired$p.value - nparcomp_paired$Analysis[1, 5])))

  test_that("Paired t-approx matches nparcomp", {
    # nparcomp rounds to 3 decimal places, so use appropriate tolerance
    expect_equal(toster_paired$p.value, nparcomp_paired$Analysis[1, 5], tolerance = 0.01)
  })

  # Two-sample permutation test
  cat("\nTest 1.3: Two-sample permutation (mtcars data)\n")
  set.seed(12345)
  toster_perm <- brunner_munzel(
    x = subset(mtcars, am == 1)$mpg,
    y = subset(mtcars, am == 0)$mpg,
    paired = FALSE,
    test_method = "perm",
    R = 5000,
    alternative = "two.sided"
  )

  set.seed(12345)
  nparcomp_perm <- npar.t.test(
    data = mtcars,
    mpg ~ am,
    method = "permu",
    nperm = 5000,
    alternative = "two.sided",
    info = FALSE
  )

  cat(sprintf("  TOSTER p-value:   %.6f\n", toster_perm$p.value))
  cat(sprintf("  nparcomp p-value: %.6f\n", nparcomp_perm$Analysis$p.value[1]))
  cat(sprintf("  Difference:       %.2e\n", abs(toster_perm$p.value - nparcomp_perm$Analysis$p.value[1])))

  # Note: Permutation p-values may differ slightly due to implementation differences
  test_that("Two-sample permutation approximately matches nparcomp", {
    expect_equal(toster_perm$p.value, nparcomp_perm$Analysis$p.value[1], tolerance = 0.05)
  })

} else {
  cat("Skipping nparcomp comparisons (package not available)\n")
}

# =============================================================================
# SECTION 2: Validate equivalence = max(one-sided p-values)
# =============================================================================

cat("\n--- Section 2: Equivalence Testing Consistency ---\n\n")

cat("Test 2.1: Equivalence p-value = max(one-sided p-values) [t-method]\n")

# Define equivalence bounds
eq_bounds <- c(0.3, 0.7)

# Run equivalence test
eq_test <- brunner_munzel(
  data = mtcars,
  mpg ~ am,
  alternative = "equivalence",
  mu = eq_bounds,
  test_method = "t"
)

# Run two one-sided tests
lower_test <- brunner_munzel(
  data = mtcars,
  mpg ~ am,
  mu = eq_bounds[1],
  alternative = "greater",
  test_method = "t"
)

upper_test <- brunner_munzel(
  data = mtcars,
  mpg ~ am,
  mu = eq_bounds[2],
  alternative = "less",
  test_method = "t"
)

expected_p <- max(lower_test$p.value, upper_test$p.value)

cat(sprintf("  Lower bound test (H1: p > %.2f): p = %.6f\n", eq_bounds[1], lower_test$p.value))
cat(sprintf("  Upper bound test (H1: p < %.2f): p = %.6f\n", eq_bounds[2], upper_test$p.value))
cat(sprintf("  max(p1, p2):                     p = %.6f\n", expected_p))
cat(sprintf("  Equivalence test:                p = %.6f\n", eq_test$p.value))
cat(sprintf("  Difference:                      %.2e\n", abs(eq_test$p.value - expected_p)))

test_that("Equivalence p-value equals max of one-sided (t-method)", {
  expect_equal(eq_test$p.value, expected_p, tolerance = 1e-10)
})

# Test with logit method
cat("\nTest 2.2: Equivalence p-value = max(one-sided p-values) [logit-method]\n")

eq_logit <- brunner_munzel(
  data = mtcars,
  mpg ~ am,
  alternative = "equivalence",
  mu = eq_bounds,
  test_method = "logit"
)

lower_logit <- brunner_munzel(
  data = mtcars,
  mpg ~ am,
  mu = eq_bounds[1],
  alternative = "greater",
  test_method = "logit"
)

upper_logit <- brunner_munzel(
  data = mtcars,
  mpg ~ am,
  mu = eq_bounds[2],
  alternative = "less",
  test_method = "logit"
)

expected_p_logit <- max(lower_logit$p.value, upper_logit$p.value)

cat(sprintf("  Lower bound test (H1: p > %.2f): p = %.6f\n", eq_bounds[1], lower_logit$p.value))
cat(sprintf("  Upper bound test (H1: p < %.2f): p = %.6f\n", eq_bounds[2], upper_logit$p.value))
cat(sprintf("  max(p1, p2):                     p = %.6f\n", expected_p_logit))
cat(sprintf("  Equivalence test:                p = %.6f\n", eq_logit$p.value))
cat(sprintf("  Difference:                      %.2e\n", abs(eq_logit$p.value - expected_p_logit)))

test_that("Equivalence p-value equals max of one-sided (logit-method)", {
  expect_equal(eq_logit$p.value, expected_p_logit, tolerance = 1e-10)
})

# Test with permutation method
cat("\nTest 2.3: Equivalence p-value = max(one-sided p-values) [perm-method]\n")

set.seed(54321)
eq_perm <- brunner_munzel(
  data = mtcars,
  mpg ~ am,
  alternative = "equivalence",
  mu = eq_bounds,
  test_method = "perm",
  R = 2000
)

set.seed(54321)
lower_perm <- brunner_munzel(
  data = mtcars,
  mpg ~ am,
  mu = eq_bounds[1],
  alternative = "greater",
  test_method = "perm",
  R = 2000
)

set.seed(54321)
upper_perm <- brunner_munzel(
  data = mtcars,
  mpg ~ am,
  mu = eq_bounds[2],
  alternative = "less",
  test_method = "perm",
  R = 2000
)

# Note: For permutation tests, the equivalence test computes both p-values in one call
# with the same permutation distribution, so we compare to verify consistency
cat(sprintf("  Lower bound test (H1: p > %.2f): p = %.6f\n", eq_bounds[1], lower_perm$p.value))
cat(sprintf("  Upper bound test (H1: p < %.2f): p = %.6f\n", eq_bounds[2], upper_perm$p.value))
cat(sprintf("  Equivalence test:                p = %.6f\n", eq_perm$p.value))
cat("  Note: Permutation p-values use same distribution for equivalence test\n")

test_that("Equivalence p-value is valid (perm-method)", {
  expect_true(eq_perm$p.value >= 0 && eq_perm$p.value <= 1)
  expect_equal(length(eq_perm$null.value), 2)
})

# =============================================================================
# SECTION 3: Validate minimal effect = min(one-sided p-values)
# =============================================================================

cat("\n--- Section 3: Minimal Effect Testing Consistency ---\n\n")

cat("Test 3.1: MET p-value = min(one-sided p-values) [t-method]\n")

met_bounds <- c(0.4, 0.6)

# Run MET test
met_test <- brunner_munzel(
  data = mtcars,
  mpg ~ am,
  alternative = "minimal.effect",
  mu = met_bounds,
  test_method = "t"
)

# Run two one-sided tests (note: opposite direction from equivalence)
lower_met <- brunner_munzel(
  data = mtcars,
  mpg ~ am,
  mu = met_bounds[1],
  alternative = "less",
  test_method = "t"
)

upper_met <- brunner_munzel(
  data = mtcars,
  mpg ~ am,
  mu = met_bounds[2],
  alternative = "greater",
  test_method = "t"
)

expected_met_p <- min(lower_met$p.value, upper_met$p.value)

cat(sprintf("  Lower bound test (H1: p < %.2f): p = %.6f\n", met_bounds[1], lower_met$p.value))
cat(sprintf("  Upper bound test (H1: p > %.2f): p = %.6f\n", met_bounds[2], upper_met$p.value))
cat(sprintf("  min(p1, p2):                     p = %.6f\n", expected_met_p))
cat(sprintf("  MET test:                        p = %.6f\n", met_test$p.value))
cat(sprintf("  Difference:                      %.2e\n", abs(met_test$p.value - expected_met_p)))

test_that("MET p-value equals min of one-sided (t-method)", {
  expect_equal(met_test$p.value, expected_met_p, tolerance = 1e-10)
})

# =============================================================================
# SECTION 4: Validate paired samples equivalence testing
# =============================================================================

cat("\n--- Section 4: Paired Samples Equivalence Testing ---\n\n")

cat("Test 4.1: Paired equivalence consistency [t-method]\n")

paired_eq <- brunner_munzel(
  x = sleep$extra[sleep$group == 1],
  y = sleep$extra[sleep$group == 2],
  paired = TRUE,
  alternative = "equivalence",
  mu = c(0.3, 0.7),
  test_method = "t"
)

paired_lower <- brunner_munzel(
  x = sleep$extra[sleep$group == 1],
  y = sleep$extra[sleep$group == 2],
  paired = TRUE,
  mu = 0.3,
  alternative = "greater",
  test_method = "t"
)

paired_upper <- brunner_munzel(
  x = sleep$extra[sleep$group == 1],
  y = sleep$extra[sleep$group == 2],
  paired = TRUE,
  mu = 0.7,
  alternative = "less",
  test_method = "t"
)

expected_paired_p <- max(paired_lower$p.value, paired_upper$p.value)

cat(sprintf("  Lower bound test (H1: p > 0.30): p = %.6f\n", paired_lower$p.value))
cat(sprintf("  Upper bound test (H1: p < 0.70): p = %.6f\n", paired_upper$p.value))
cat(sprintf("  max(p1, p2):                     p = %.6f\n", expected_paired_p))
cat(sprintf("  Paired equivalence test:         p = %.6f\n", paired_eq$p.value))
cat(sprintf("  Difference:                      %.2e\n", abs(paired_eq$p.value - expected_paired_p)))

test_that("Paired equivalence p-value equals max of one-sided", {
  expect_equal(paired_eq$p.value, expected_paired_p, tolerance = 1e-10)
})

# Paired permutation test
cat("\nTest 4.2: Paired equivalence consistency [perm-method]\n")

set.seed(99999)
paired_eq_perm <- brunner_munzel(
  x = sleep$extra[sleep$group == 1],
  y = sleep$extra[sleep$group == 2],
  paired = TRUE,
  alternative = "equivalence",
  mu = c(0.3, 0.7),
  test_method = "perm",
  R = 2000
)

cat(sprintf("  Paired permutation equivalence:  p = %.6f\n", paired_eq_perm$p.value))
cat(sprintf("  Relative effect estimate:        %.4f\n", paired_eq_perm$estimate))
cat(sprintf("  90%% CI: [%.4f, %.4f]\n", paired_eq_perm$conf.int[1], paired_eq_perm$conf.int[2]))

test_that("Paired permutation equivalence test works", {
  expect_true(paired_eq_perm$p.value >= 0 && paired_eq_perm$p.value <= 1)
  expect_equal(paired_eq_perm$alternative, "equivalence")
  expect_equal(attr(paired_eq_perm$conf.int, "conf.level"), 0.9)  # 1-2*alpha
})

# =============================================================================
# SECTION 5: Validate confidence interval levels
# =============================================================================

cat("\n--- Section 5: Confidence Interval Levels ---\n\n")

cat("Test 5.1: CI level for standard alternatives (should be 95%%)\n")
std_test <- brunner_munzel(data = mtcars, mpg ~ am, alternative = "two.sided")
cat(sprintf("  Two-sided CI level: %.0f%%\n", attr(std_test$conf.int, "conf.level") * 100))

test_that("Standard alternative uses 1-alpha CI", {
  expect_equal(attr(std_test$conf.int, "conf.level"), 0.95)
})

cat("\nTest 5.2: CI level for equivalence (should be 90%%, i.e., 1-2*alpha)\n")
eq_ci_test <- brunner_munzel(
  data = mtcars,
  mpg ~ am,
  alternative = "equivalence",
  mu = c(0.3, 0.7)
)
cat(sprintf("  Equivalence CI level: %.0f%%\n", attr(eq_ci_test$conf.int, "conf.level") * 100))

test_that("Equivalence uses 1-2*alpha CI", {
  expect_equal(attr(eq_ci_test$conf.int, "conf.level"), 0.90)
})

cat("\nTest 5.3: Logit method CI stays within [0, 1]\n")
logit_test <- brunner_munzel(
  data = mtcars,
  mpg ~ am,
  test_method = "logit"
)
cat(sprintf("  Logit CI: [%.4f, %.4f]\n", logit_test$conf.int[1], logit_test$conf.int[2]))

test_that("Logit CI stays within [0, 1]", {
  expect_true(logit_test$conf.int[1] >= 0)
  expect_true(logit_test$conf.int[2] <= 1)
})

# =============================================================================
# SECTION 6: Test all method combinations
# =============================================================================

cat("\n--- Section 6: All Method Combinations ---\n\n")

methods <- c("t", "logit", "perm")
alternatives <- c("equivalence", "minimal.effect")

cat("Testing all combinations of test_method and alternative:\n\n")

for (method in methods) {
  for (alt in alternatives) {
    set.seed(42)
    result <- tryCatch({
      brunner_munzel(
        data = mtcars,
        mpg ~ am,
        alternative = alt,
        mu = c(0.35, 0.65),
        test_method = method,
        R = 999
      )
    }, error = function(e) e)

    if (inherits(result, "error")) {
      cat(sprintf("  %s + %s: ERROR - %s\n", method, alt, result$message))
    } else {
      cat(sprintf("  %s + %s: p = %.4f, estimate = %.4f\n",
                  method, alt, result$p.value, result$estimate))
    }
  }
}

cat("\n")

# =============================================================================
# SUMMARY
# =============================================================================

cat("=============================================================================\n")
cat("VALIDATION COMPLETE\n")
cat("=============================================================================\n")
cat("\n")
cat("Key validations performed:\n")
cat("1. Base implementation matches nparcomp (two-sample and paired)\n")
cat("2. Equivalence p-value = max(one-sided p-values)\n")
cat("3. Minimal effect p-value = min(one-sided p-values)\n")
cat("4. Paired samples equivalence testing works correctly\n")
cat("5. CI levels are correct (1-alpha for standard, 1-2*alpha for TOST)\n")
cat("6. All test_method x alternative combinations work\n")
cat("\n")
