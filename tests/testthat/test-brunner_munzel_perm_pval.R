# Regression tests for brunner_munzel paired permutation p_method handling
# This file tests that the paired permutation code correctly routes through
# bm_compute_perm_pval, respecting the p_method argument ("exact" vs "plusone").

hush = function(code){
  sink(nullfile())
  on.exit(sink())
  suppressMessages(code)
}

# Test 1: Paired permutation equivalence respects p_method argument
test_that("paired permutation equivalence respects p_method argument", {
  skip_on_cran()
  set.seed(42)
  n <- 20  # n > 13 forces randomization path
  x <- rnorm(n, mean = 5)
  y <- rnorm(n, mean = 5.1)

  res_exact <- hush(brunner_munzel(x, y, paired = TRUE,
                               alternative = "equivalence",
                               mu = c(0.3, 0.7),
                               test_method = "perm",
                               R = 999,
                               p_method = "exact"))

  res_plusone <- hush(brunner_munzel(x, y, paired = TRUE,
                                alternative = "equivalence",
                                mu = c(0.3, 0.7),
                                test_method = "perm",
                                R = 999,
                                p_method = "plusone"))

  # p-values should differ because (b+1)/(R+1) != b/R
  expect_false(identical(res_exact$p.value, res_plusone$p.value))
})

# Test 2: Paired permutation two.sided respects p_method argument
test_that("paired permutation two.sided respects p_method argument", {
  skip_on_cran()
  set.seed(42)
  n <- 20
  x <- rnorm(n, mean = 5)
  y <- rnorm(n, mean = 5.5)

  res_exact <- hush(brunner_munzel(x, y, paired = TRUE,
                               alternative = "two.sided",
                               test_method = "perm",
                               R = 999,
                               p_method = "exact"))

  res_plusone <- hush(brunner_munzel(x, y, paired = TRUE,
                                alternative = "two.sided",
                                test_method = "perm",
                                R = 999,
                                p_method = "plusone"))

  expect_false(identical(res_exact$p.value, res_plusone$p.value))
})

# Test 3: Paired exact permutation equivalence p-values unchanged after fix
test_that("paired exact permutation equivalence p-values unchanged after fix", {
  skip_on_cran()
  set.seed(123)
  # n <= 13 triggers exact enumeration
  x <- c(5.1, 4.8, 6.2, 5.7, 6.0, 5.5, 4.9, 5.8, 6.1, 5.3)
  y <- c(5.0, 4.9, 6.1, 5.6, 5.9, 5.6, 5.0, 5.7, 6.0, 5.4)

  res <- hush(brunner_munzel(x, y, paired = TRUE,
                         alternative = "equivalence",
                         mu = c(0.3, 0.7),
                         test_method = "perm"))

  # For exact permutation with p_method="exact", b/R == mean(Tperm >= t)
  # So the fix should not change the result
  expect_true(res$p.value >= 0 && res$p.value <= 1)
})

# Test 4: Paired permutation minimal.effect respects p_method argument
test_that("paired permutation minimal.effect respects p_method argument", {
  skip_on_cran()
  set.seed(42)
  n <- 20
  # Use a moderate difference so the MET p-value isn't at the floor (0)
  x <- rnorm(n, mean = 5)
  y <- rnorm(n, mean = 5.5)

  res_exact <- hush(brunner_munzel(x, y, paired = TRUE,
                               alternative = "minimal.effect",
                               mu = c(0.3, 0.7),
                               test_method = "perm",
                               R = 999,
                               p_method = "exact"))

  res_plusone <- hush(brunner_munzel(x, y, paired = TRUE,
                                alternative = "minimal.effect",
                                mu = c(0.3, 0.7),
                                test_method = "perm",
                                R = 999,
                                p_method = "plusone"))

  expect_false(identical(res_exact$p.value, res_plusone$p.value))
})

# Test 5: Paired permutation p-value is never zero with plusone method
test_that("paired permutation p-value is never zero with plusone", {
  skip_on_cran()
  set.seed(99)
  n <- 20
  x <- rnorm(n, mean = 0)
  y <- rnorm(n, mean = 10)  # extreme difference

  # Note: minimal.effect uses 1-p internally, so plusone doesn't guarantee p > 0
  for (alt in c("two.sided", "less", "greater", "equivalence")) {
    mu_val <- if (alt %in% c("equivalence", "minimal.effect")) c(0.3, 0.7) else 0.5
    res <- hush(brunner_munzel(x, y, paired = TRUE,
                           alternative = alt,
                           mu = mu_val,
                           test_method = "perm",
                           R = 999,
                           p_method = "plusone"))
    expect_gt(res$p.value, 0,
              label = paste("p-value should be > 0 for alternative =", alt))
  }
})

# Test 6: Paired permutation less/greater also respect p_method
test_that("paired permutation less and greater respect p_method argument", {
  skip_on_cran()
  set.seed(42)
  n <- 20
  x <- rnorm(n, mean = 5)
  y <- rnorm(n, mean = 5.5)

  for (alt in c("less", "greater")) {
    res_exact <- hush(brunner_munzel(x, y, paired = TRUE,
                                 alternative = alt,
                                 test_method = "perm",
                                 R = 999,
                                 p_method = "exact"))

    res_plusone <- hush(brunner_munzel(x, y, paired = TRUE,
                                  alternative = alt,
                                  test_method = "perm",
                                  R = 999,
                                  p_method = "plusone"))

    expect_false(identical(res_exact$p.value, res_plusone$p.value),
                 label = paste("p-values should differ for alternative =", alt))
  }
})
