# Regression tests for brunner_munzel degenerate/edge cases
# These tests guard against the issues identified in diagnostic_report.md

hush = function(code){
  sink(nullfile())
  on.exit(sink())
  suppressMessages(suppressWarnings(code))
}

# --- Fix 1: Permutation CI index-0 bug (CRITICAL) ---

test_that("permutation CI does not contain NA for small exact tests", {
  # With n=3 per group, R_actual = choose(6,3) = 20.

  # floor(alpha/2 * 20) = floor(0.5) = 0, which caused index-0 -> NA
  res <- hush(brunner_munzel(
    x = c(0, 0, 0), y = c(0, 0, 0),
    test_method = "perm", alternative = "two.sided"
  ))
  expect_false(any(is.na(res$conf.int)))
  expect_true(all(is.finite(res$conf.int)))
})

test_that("permutation CI does not contain NA for complete separation", {
  res <- hush(brunner_munzel(
    x = c(1, 2, 3), y = c(10, 11, 12),
    test_method = "perm", alternative = "two.sided"
  ))
  expect_false(any(is.na(res$conf.int)))
  expect_true(all(is.finite(res$conf.int)))
})

# --- Fix 2: pd clamping restricted to logit path ---

test_that("estimate is not clamped for t-approx when pd is 0 or 1", {
  # Complete separation: all x < all y, so true pd = 0
  res <- hush(brunner_munzel(
    x = c(1, 1, 1, 1), y = c(2, 2, 2, 2),
    test_method = "t"
  ))
  expect_equal(unname(res$estimate), 0)
})

test_that("estimate is not clamped for perm when pd is 0 or 1", {
  res <- hush(brunner_munzel(
    x = c(1, 1, 1, 1), y = c(2, 2, 2, 2),
    test_method = "perm"
  ))
  expect_equal(unname(res$estimate), 0)
})

test_that("logit method warns when pd is exactly 0 or 1", {
  expect_warning(
    hush2 <- suppressMessages(
      brunner_munzel(
        x = c(1, 1, 1, 1), y = c(2, 2, 2, 2),
        test_method = "logit"
      )
    ),
    "Logit method is unreliable"
  )
})

# --- Fix 3: Paired CI clamped to [0, 1] ---

test_that("paired t-approx CI is within [0, 1]", {
  res <- hush(brunner_munzel(
    x = c(0, 0, 0), y = c(0, 0, 0),
    paired = TRUE, test_method = "t"
  ))
  expect_gte(res$conf.int[1], 0)
  expect_lte(res$conf.int[2], 1)
})

test_that("paired t-approx CI is within [0, 1] for varied identical data", {
  res <- hush(brunner_munzel(
    x = c(1, 2, 3, 4, 5), y = c(1, 2, 3, 4, 5),
    paired = TRUE, test_method = "t"
  ))
  expect_gte(res$conf.int[1], 0)
  expect_lte(res$conf.int[2], 1)
})

# --- Fix 4: Minimum sample size guard ---

test_that("n < 3 per group gives informative error", {
  expect_error(
    brunner_munzel(x = c(1), y = c(2), test_method = "t"),
    "at least 3 observations"
  )
  expect_error(
    brunner_munzel(x = c(1, 2), y = c(3, 4), test_method = "t"),
    "at least 3 observations"
  )
  expect_error(
    brunner_munzel(x = c(1), y = c(2), test_method = "perm"),
    "at least 3 observations"
  )
})

# --- Fix 5: Permutation p = 0 warning ---

test_that("exact enumeration never gives p = 0 for zero-inflated data", {
  skip_on_cran()
  # With full enumeration the observed arrangement is one of the permutations,
  # so b >= 1. Previously floating point differences between the observed and
  # permuted statistics made every permutation look less extreme (p = 0).
  # Here the single 1 lands in x in half of the permutations and in y in the
  # other half, so all |T*| equal |T| (two-sided p = 1) and half of T* >= T.
  x <- c(0, 0, 0, 0, 0, 1)
  y <- c(0, 0, 0, 0, 0, 0)
  expect_no_warning(
    res <- suppressMessages(brunner_munzel(x, y, test_method = "perm"))
  )
  expect_equal(res$p.value, 1)
  res_g <- suppressMessages(brunner_munzel(x, y, test_method = "perm",
                                           alternative = "greater"))
  expect_equal(res_g$p.value, 0.5)
})

test_that("exact p_method with sampled permutations can warn about p = 0", {
  skip_on_cran()
  # p = 0 is only possible when permutations are sampled and b/R is requested
  set.seed(1)
  expect_warning(
    suppressMessages(brunner_munzel(x = 1:10 + 20, y = 1:10, test_method = "perm",
                                    R = 200, p_method = "exact")),
    "Exact permutation p-value is 0"
  )
})

# --- Identical-samples sanity checks (the original brunnermunzel bug) ---

test_that("identical samples return p = 1 and t = 0 for all methods", {
  for (method in c("t", "logit", "perm")) {
    res <- hush(brunner_munzel(
      x = c(0, 0, 0), y = c(0, 0, 0),
      test_method = method
    ))
    expect_equal(unname(res$statistic), 0,
                 info = paste("statistic should be 0 for method =", method))
    expect_equal(res$p.value, 1,
                 info = paste("p.value should be 1 for method =", method))
    expect_equal(unname(res$estimate), 0.5,
                 info = paste("estimate should be 0.5 for method =", method))
  }
})

test_that("larger identical samples return p = 1 for all methods", {
  for (method in c("t", "logit", "perm")) {
    res <- hush(brunner_munzel(
      x = rep(5, 10), y = rep(5, 10),
      test_method = method
    ))
    expect_equal(res$p.value, 1,
                 info = paste("p.value should be 1 for method =", method))
    expect_equal(unname(res$estimate), 0.5,
                 info = paste("estimate should be 0.5 for method =", method))
  }
})

# --- CI clamping message ---

test_that("clamping message is emitted when CI exceeds [0, 1]", {
  expect_message(
    suppressWarnings(
      brunner_munzel(
        x = c(0, 0, 0), y = c(0, 0, 0),
        paired = TRUE, test_method = "t"
      )
    ),
    "clamped to the \\[0, 1\\] range"
  )
})
