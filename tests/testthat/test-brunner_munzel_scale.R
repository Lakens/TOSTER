# test-brunner_munzel_scale.R
# brunner_munzel scale argument --------

test_that("brunner_munzel scale='probability' is default and unchanged", {
  res_default <- brunner_munzel(mpg ~ am, data = mtcars)
  res_explicit <- brunner_munzel(mpg ~ am, data = mtcars, scale = "probability")
  expect_equal(res_default$estimate, res_explicit$estimate)
  expect_equal(res_default$conf.int, res_explicit$conf.int)
  expect_equal(res_default$stderr, res_explicit$stderr)
  expect_equal(res_default$p.value, res_explicit$p.value)
})

test_that("brunner_munzel scale='difference' transforms output correctly", {
  res_prob <- brunner_munzel(mpg ~ am, data = mtcars)
  res_diff <- brunner_munzel(mpg ~ am, data = mtcars, scale = "difference")

  expect_equal(as.numeric(res_diff$estimate), 2 * as.numeric(res_prob$estimate) - 1)
  expect_equal(as.numeric(res_diff$conf.int), 2 * as.numeric(res_prob$conf.int) - 1)
  expect_equal(as.numeric(res_diff$stderr), 2 * as.numeric(res_prob$stderr))

  # p-value and test statistic should be identical
  expect_equal(res_diff$p.value, res_prob$p.value)
  expect_equal(res_diff$statistic, res_prob$statistic)
  expect_equal(res_diff$parameter, res_prob$parameter)
})

test_that("brunner_munzel scale='logodds' transforms output correctly", {
  res_prob <- brunner_munzel(mpg ~ am, data = mtcars)
  res_lo <- brunner_munzel(mpg ~ am, data = mtcars, scale = "logodds")

  p <- as.numeric(res_prob$estimate)
  expect_equal(as.numeric(res_lo$estimate), qlogis(p), tolerance = 1e-10)
  expect_equal(as.numeric(res_lo$null.value), 0)  # logit(0.5) = 0

  # p-value unchanged
  expect_equal(res_lo$p.value, res_prob$p.value)
})

test_that("brunner_munzel scale='odds' transforms output correctly", {
  res_prob <- brunner_munzel(mpg ~ am, data = mtcars)
  res_odds <- brunner_munzel(mpg ~ am, data = mtcars, scale = "odds")

  p <- as.numeric(res_prob$estimate)
  expect_equal(as.numeric(res_odds$estimate), p / (1 - p), tolerance = 1e-10)
  expect_equal(as.numeric(res_odds$null.value), 1)  # 0.5/(1-0.5) = 1

  # p-value unchanged
  expect_equal(res_odds$p.value, res_prob$p.value)
})

test_that("brunner_munzel scale works with paired samples", {
  x <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
         11, 12, 13, 14, 15, 16, 17, 18, 19, 20)
  y <- c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
         12, 13, 14, 15, 16, 17, 18, 19, 20, 21)
  res <- brunner_munzel(x, y, paired = TRUE, scale = "difference")
  # Check paired label uses difference notation: 'x' - 'y'
  expect_true(grepl("'x' - 'y'", names(res$estimate), fixed = TRUE))
})

test_that("brunner_munzel scale works with equivalence tests", {
  res_prob <- brunner_munzel(mpg ~ am, data = mtcars,
                             alternative = "equivalence", mu = c(0.35, 0.65))
  res_diff <- brunner_munzel(mpg ~ am, data = mtcars,
                             alternative = "equivalence", mu = c(0.35, 0.65),
                             scale = "difference")
  # Null values should be transformed
  expect_equal(as.numeric(res_diff$null.value), c(-0.3, 0.3))
  # p-value unchanged
  expect_equal(res_diff$p.value, res_prob$p.value)
})

test_that("brunner_munzel scale works with permutation tests", {
  set.seed(42)
  res_prob <- brunner_munzel(mpg ~ am, data = mtcars, test_method = "perm",
                             scale = "probability")
  set.seed(42)
  res_diff <- brunner_munzel(mpg ~ am, data = mtcars, test_method = "perm",
                             scale = "difference")
  expect_equal(as.numeric(res_diff$estimate), 2 * as.numeric(res_prob$estimate) - 1)
})

test_that("brunner_munzel scale works with logit test_method", {
  res <- brunner_munzel(mpg ~ am, data = mtcars,
                        test_method = "logit", scale = "odds")
  expect_true(grepl("odds(", names(res$estimate), fixed = TRUE))
  # The test_method=logit affects inference; scale=odds affects reporting
  # Both should work independently
})

test_that("brunner_munzel formula method respects scale with group names", {
  res <- brunner_munzel(mpg ~ am, data = mtcars, scale = "difference")
  # Should contain quoted factor level names
  expect_true(grepl("'0'", names(res$estimate), fixed = TRUE))
  expect_true(grepl("'1'", names(res$estimate), fixed = TRUE))
  expect_true(grepl("P(", names(res$estimate), fixed = TRUE))
})

test_that("brunner_munzel default method uses quoted X/Y in labels", {
  x <- rnorm(20)
  y <- rnorm(20, mean = 1)
  res <- brunner_munzel(x, y)
  # Default method uses generic 'X' and 'Y' (quoted)
  expect_true(grepl("'X'", names(res$estimate), fixed = TRUE))
  expect_true(grepl("'Y'", names(res$estimate), fixed = TRUE))
})

test_that("brunner_munzel label quoting is consistent across scales", {
  res_prob <- brunner_munzel(mpg ~ am, data = mtcars, scale = "probability")
  res_diff <- brunner_munzel(mpg ~ am, data = mtcars, scale = "difference")
  res_lo <- brunner_munzel(mpg ~ am, data = mtcars, scale = "logodds")
  res_odds <- brunner_munzel(mpg ~ am, data = mtcars, scale = "odds")

  # All should have quoted group names
  for (res in list(res_prob, res_diff, res_lo, res_odds)) {
    expect_true(grepl("'0'", names(res$estimate), fixed = TRUE))
    expect_true(grepl("'1'", names(res$estimate), fixed = TRUE))
  }
})
