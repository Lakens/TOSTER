# CRAN-tier smoke tests for resampling functions --------
# Fast checks (small R) so resampling code paths still run on CRAN; the
# thorough tests for these functions use skip_on_cran().

test_that("boot_cor_test smoke test", {
  set.seed(4321)
  x <- rnorm(30)
  y <- x + rnorm(30)

  res <- boot_cor_test(x, y, R = 99)
  expect_s3_class(res, "htest")
  expect_length(res$conf.int, 2)
  expect_true(res$p.value >= 0 && res$p.value <= 1)

  res_eq <- boot_cor_test(x, y, alternative = "equivalence", null = 0.3,
                          method = "spearman", R = 99)
  expect_s3_class(res_eq, "htest")
  expect_true(res_eq$p.value >= 0 && res_eq$p.value <= 1)
})

test_that("boot_compare_cor smoke test", {
  set.seed(4322)
  x1 <- rnorm(30)
  y1 <- x1 + rnorm(30)
  x2 <- rnorm(30)
  y2 <- rnorm(30)

  res <- boot_compare_cor(x1 = x1, y1 = y1, x2 = x2, y2 = y2, R = 99)
  expect_s3_class(res, "htest")
  expect_length(res$conf.int, 2)
  expect_true(res$p.value >= 0 && res$p.value <= 1)
})

test_that("boot_log_TOST smoke test", {
  set.seed(4323)
  dat <- data.frame(y = rlnorm(40, meanlog = 1),
                    group = rep(c("a", "b"), each = 20))

  res <- boot_log_TOST(y ~ group, data = dat, R = 99)
  expect_s3_class(res, "TOSTt")
  expect_true(all(res$TOST$p.value >= 0 & res$TOST$p.value <= 1))

  res_paired <- boot_log_TOST(x = dat$y[1:20], y = dat$y[21:40],
                              paired = TRUE, boot_ci = "perc", R = 99)
  expect_s3_class(res_paired, "TOSTt")
})
