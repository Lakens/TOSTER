hush <- function(code) {
  suppressWarnings(suppressMessages(capture.output(code)))
}

# Test 1: Standard two-sided test at mu = 0 returns a plausible p-value --------
test_that("two-sided test at mu = 0 returns plausible result", {
  set.seed(4321)
  # Use large separation to ensure a clearly significant result
  x <- rnorm(40, mean = 0)
  y <- rnorm(40, mean = 2)

  res <- suppressWarnings(
    boot_ses_test(x = x, y = y, ses = "rb",
                  mu = 0, alternative = "two.sided", B = 999)
  )

  expect_s3_class(res, "htest")
  expect_true(res$p.value >= 0 && res$p.value <= 1)
  # Groups differ substantially, so p-value should be small
  expect_true(res$p.value < 0.05)
  expect_true(!is.null(res$estimate))
  expect_true(!is.null(res$model.param))
})

# Test 2: Equivalence test where obs_rb is clearly inside bounds --------
test_that("equivalence test with obs_rb inside bounds yields small p-value", {
  set.seed(1234)
  # Generate data with no real difference, very wide bounds
  x <- rnorm(50)
  y <- rnorm(50)

  res <- suppressWarnings(
    boot_ses_test(x = x, y = y, ses = "rb",
                  mu = c(-0.9, 0.9), alternative = "equivalence", B = 999)
  )

  expect_s3_class(res, "htest")
  expect_true(res$p.value < 0.05)
  expect_equal(res$alternative, "equivalence")
})

# Test 3: Equivalence test where obs_rb is clearly outside bounds --------
test_that("equivalence test with obs_rb outside bounds yields large p-value", {
  set.seed(5678)
  # Generate data with large difference, narrow bounds
  x <- rnorm(30, mean = 0)
  y <- rnorm(30, mean = 2)

  res <- suppressWarnings(
    boot_ses_test(x = x, y = y, ses = "rb",
                  mu = c(-0.1, 0.1), alternative = "equivalence", B = 599)
  )

  expect_s3_class(res, "htest")
  expect_true(res$p.value > 0.50)
})

# Test 4: Complete separation triggers a warning --------
test_that("complete separation triggers a warning", {
  # Use n >= 20 per group to avoid triggering the small-n warning
  x <- 1:25
  y <- 26:50

  expect_warning(
    boot_ses_test(x = x, y = y, ses = "rb",
                  mu = 0, alternative = "two.sided", B = 100),
    "Complete or near-complete separation"
  )
})

# Test 5: Small n triggers a warning --------
test_that("small n triggers a warning", {
  set.seed(99)
  x <- rnorm(10)
  y <- rnorm(10, mean = 0.5)

  expect_warning(
    boot_ses_test(x = x, y = y, ses = "rb",
                  mu = 0, alternative = "two.sided", B = 100),
    "Sample size is small"
  )
})

# Test 6: Invalid mu outside the valid ses range throws an error --------
test_that("invalid mu outside valid range throws error", {
  x <- rnorm(30)
  y <- rnorm(30)

  # rb must be in (-1, 1)
  expect_error(
    boot_ses_test(x = x, y = y, ses = "rb",
                  mu = 1.5, alternative = "two.sided", B = 100),
    "mu must be in the open interval"
  )

  # cstat must be in (0, 1)
  expect_error(
    boot_ses_test(x = x, y = y, ses = "cstat",
                  mu = -0.1, alternative = "two.sided", B = 100),
    "mu must be in the open interval"
  )

  # odds must be positive
  expect_error(
    boot_ses_test(x = x, y = y, ses = "odds",
                  mu = -1, alternative = "two.sided", B = 100),
    "mu must be positive"
  )
})

# Test 7: keep_boot = FALSE returns NULL for boot distribution fields --------
test_that("keep_boot = FALSE returns NULL for boot distributions", {
  set.seed(111)
  x <- rnorm(25)
  y <- rnorm(25, mean = 0.5)

  res <- suppressWarnings(
    boot_ses_test(x = x, y = y, ses = "rb",
                  mu = 0, alternative = "two.sided",
                  B = 100, keep_boot = FALSE)
  )

  expect_null(res$boot.dist)
  expect_null(res$boot.dist.low)
  expect_null(res$boot.dist.high)

  # Equivalence test with keep_boot = FALSE
  res_eq <- suppressWarnings(
    boot_ses_test(x = x, y = y, ses = "rb",
                  mu = c(-0.5, 0.5), alternative = "equivalence",
                  B = 100, keep_boot = FALSE)
  )

  expect_null(res_eq$boot.dist)
  expect_null(res_eq$boot.dist.low)
  expect_null(res_eq$boot.dist.high)
})

# Test 8: Formula interface produces same result as default interface --------
test_that("formula interface matches default interface", {
  set.seed(2222)
  dat <- data.frame(
    value = c(rnorm(25, 0), rnorm(25, 0.5)),
    group = factor(rep(c("a", "b"), each = 25))
  )

  x <- dat$value[dat$group == "a"]
  y <- dat$value[dat$group == "b"]

  res_default <- suppressWarnings(
    boot_ses_test(x = x, y = y, ses = "rb",
                  mu = 0, alternative = "two.sided", B = 599)
  )

  res_formula <- suppressWarnings(
    boot_ses_test(formula = value ~ group, data = dat,
                  ses = "rb", mu = 0,
                  alternative = "two.sided", B = 599)
  )

  # Estimates should be identical (same data, same method)
  expect_equal(res_default$estimate, res_formula$estimate, tolerance = 1e-10)
  # p-values may differ due to stochastic bootstrap if seed isn't perfectly aligned,
  # but the method and alternative should match
  expect_equal(res_default$method, res_formula$method)
  expect_equal(res_default$alternative, res_formula$alternative)
})

# Test 9: conf.int is NOT present in returned object --------
test_that("conf.int is intentionally absent from output", {
  set.seed(333)
  x <- rnorm(25)
  y <- rnorm(25, mean = 0.3)

  res <- suppressWarnings(
    boot_ses_test(x = x, y = y, ses = "rb",
                  mu = 0, alternative = "two.sided", B = 100)
  )

  expect_null(res$conf.int)
})

# Test 10: One-sample design is not supported --------
test_that("one-sample design errors", {
  expect_error(
    boot_ses_test(x = rnorm(20), ses = "rb",
                  mu = 0, alternative = "two.sided", B = 100),
    "requires two samples"
  )
})

# Test 11: B must be >= 100 --------
test_that("B must be at least 100", {
  x <- rnorm(20)
  y <- rnorm(20)
  expect_error(
    boot_ses_test(x = x, y = y, ses = "rb",
                  mu = 0, alternative = "two.sided", B = 50),
    "B must be an integer >= 100"
  )
})

# Test 12: Paired samples work --------
test_that("paired samples produce valid result", {
  set.seed(444)
  pre  <- rnorm(25, mean = 5)
  post <- pre + rnorm(25, mean = 0.3, sd = 0.5)

  res <- suppressWarnings(
    boot_ses_test(x = pre, y = post, paired = TRUE,
                  ses = "rb", mu = 0, alternative = "two.sided", B = 599)
  )

  expect_s3_class(res, "htest")
  expect_true(res$p.value >= 0 && res$p.value <= 1)
  expect_true(grepl("Paired", res$method))
})

# Test 13: Different ses scales produce valid results --------
test_that("all ses scales produce valid output", {
  set.seed(555)
  x <- rnorm(25)
  y <- rnorm(25, mean = 0.5)

  for (scale in c("rb", "cstat", "odds", "logodds")) {
    mu_val <- switch(scale,
                     "rb" = 0, "cstat" = 0.5, "odds" = 1, "logodds" = 0)
    res <- suppressWarnings(
      boot_ses_test(x = x, y = y, ses = scale,
                    mu = mu_val, alternative = "two.sided", B = 200)
    )
    expect_s3_class(res, "htest")
    expect_true(res$p.value >= 0 && res$p.value <= 1)
  }
})

# Test 14: Minimal effect test --------
test_that("minimal effect test returns valid result", {
  set.seed(666)
  # Very large group separation, narrow bounds: obs_rb should be well outside
  x <- rnorm(40, mean = 0)
  y <- rnorm(40, mean = 3)

  res <- suppressWarnings(
    boot_ses_test(x = x, y = y, ses = "rb",
                  mu = c(-0.2, 0.2), alternative = "minimal.effect", B = 999)
  )

  expect_s3_class(res, "htest")
  expect_equal(res$alternative, "minimal.effect")
  # Large effect well outside narrow bounds: should reject equivalence
  expect_true(res$p.value < 0.05)
})

# Test 15: keep_boot = TRUE returns bootstrap distributions --------
test_that("keep_boot = TRUE returns bootstrap distributions", {
  set.seed(777)
  x <- rnorm(25)
  y <- rnorm(25, mean = 0.3)

  # Standard alternative
  res <- suppressWarnings(
    boot_ses_test(x = x, y = y, ses = "rb",
                  mu = 0, alternative = "two.sided",
                  B = 200, keep_boot = TRUE)
  )
  expect_true(!is.null(res$boot.dist))
  expect_equal(length(res$boot.dist), 200)

  # Equivalence
  res_eq <- suppressWarnings(
    boot_ses_test(x = x, y = y, ses = "rb",
                  mu = c(-0.5, 0.5), alternative = "equivalence",
                  B = 200, keep_boot = TRUE)
  )
  expect_true(!is.null(res_eq$boot.dist.low))
  expect_true(!is.null(res_eq$boot.dist.high))
  expect_equal(length(res_eq$boot.dist.low), 200)
  expect_equal(length(res_eq$boot.dist.high), 200)
})
