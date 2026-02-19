# Tests for rank_diff (Kornbrot 1990 rank difference transformation)

# === Input validation ===

test_that("rank_diff errors when not numeric", {
  expect_error(rank_diff("a", "b"), "must be numeric")
  expect_error(rank_diff(1:5, letters[1:5]), "must be numeric")
})

test_that("rank_diff errors when lengths differ", {
  expect_error(rank_diff(1:5, 1:3), "same length")
})

test_that("rank_diff errors when names has wrong length", {
  expect_error(rank_diff(1:5, 6:10, names = "x"), "length 2")
})

test_that("rank_diff errors with no complete pairs", {
  expect_error(rank_diff(c(NA_real_, NA_real_), c(NA_real_, NA_real_)),
               "No complete pairs")
})

# === Basic functionality ===

test_that("rank_diff returns data frame with correct structure", {
  rd <- rank_diff(1:5, 6:10)
  expect_true(is.data.frame(rd))
  expect_equal(ncol(rd), 2)
  expect_equal(nrow(rd), 5)
  expect_equal(names(rd), c("x", "y"))
})

test_that("rank_diff respects custom names", {
  rd <- rank_diff(1:5, 6:10, names = c("pre", "post"))
  expect_equal(names(rd), c("pre", "post"))
})

test_that("rank_diff produces valid ranks", {
  x <- c(3, 1, 5, 2, 4)
  y <- c(8, 6, 10, 7, 9)
  rd <- rank_diff(x, y)

  # All values should be in [1, 2*n]
  all_ranks <- c(rd$x, rd$y)
  expect_true(all(all_ranks >= 1))
  expect_true(all(all_ranks <= 10))

  # Ranks should sum to n*(2n+1)/2 = 5*11/2 = 55
  expect_equal(sum(all_ranks), 55)
})

# === NA handling ===

test_that("rank_diff removes pairs with NAs and messages", {
  x <- c(1, NA, 3, 4, 5)
  y <- c(6, 7, NA, 9, 10)

  expect_message(rd <- rank_diff(x, y), "2 pairs with missing values removed")
  expect_equal(nrow(rd), 3)
  # Ranks should be based on 6 values: {1, 4, 5, 6, 9, 10}
  expect_equal(sum(c(rd$x, rd$y)), 21)  # 1+2+3+4+5+6 = 21
})

# === Monotone invariance (key property from Kornbrot 1990) ===

test_that("rank_diff gives same result for monotone transformations", {
  set.seed(42)
  x <- abs(rnorm(15, 5, 2))  # All positive for log/exp

  y <- abs(rnorm(15, 6, 2))

  rd_raw <- rank_diff(x, y)
  rd_log <- rank_diff(log(x), log(y))
  rd_exp <- rank_diff(exp(x), exp(y))
  rd_sq  <- rank_diff(x^2, y^2)

  # All monotone transformations should yield identical ranks

  expect_equal(rd_raw, rd_log)
  expect_equal(rd_raw, rd_exp)
  expect_equal(rd_raw, rd_sq)
})

# === Kornbrot (1990) Tables 1-2 consistency test ===

test_that("rank_diff: time and rate give same effect size (Kornbrot Tables 1-2)", {
  time_plac <- c(4.6, 4.3, 6.7, 5.8, 5.0, 4.2, 6.0,
                 2.0, 2.6, 10.0, 3.4, 7.1, 8.6)
  time_drug <- c(2.9, 2.8, 12.0, 3.8, 5.9, 6.5, 3.3,
                 2.3, 2.1, 14.3, 2.4, 14.0, 4.9)
  rate_plac <- 60 / time_plac
  rate_drug <- 60 / time_drug

  # Standard approach: time and rate give different absolute results
  res_time_std <- ses_calc(time_plac, time_drug, paired = TRUE, ses = "rb")
  res_rate_std <- ses_calc(rate_plac, rate_drug, paired = TRUE, ses = "rb")
  expect_false(isTRUE(all.equal(abs(unname(res_time_std$estimate)),
                                 abs(unname(res_rate_std$estimate)))))

  # Rank difference: time and rate give same absolute rb

  # (sign flips because 60/x is a decreasing monotone transform)
  rd_time <- rank_diff(time_plac, time_drug)
  rd_rate <- rank_diff(rate_plac, rate_drug)

  res_time_rd <- ses_calc(rd_time$x, rd_time$y, paired = TRUE, ses = "rb")
  res_rate_rd <- ses_calc(rd_rate$x, rd_rate$y, paired = TRUE, ses = "rb")
  expect_equal(abs(unname(res_time_rd$estimate)),
               abs(unname(res_rate_rd$estimate)))
})

test_that("rank_diff: increasing monotone transform gives identical results", {
  set.seed(42)
  x <- abs(rnorm(15, 5, 2))
  y <- abs(rnorm(15, 6, 2))

  # log is increasing on (0, Inf), so sign is preserved
  rd_raw <- rank_diff(x, y)
  rd_log <- rank_diff(log(x), log(y))

  res_raw <- ses_calc(rd_raw$x, rd_raw$y, paired = TRUE, ses = "rb")
  res_log <- ses_calc(rd_log$x, rd_log$y, paired = TRUE, ses = "rb")
  expect_equal(unname(res_raw$estimate), unname(res_log$estimate))
})

# === Works with ses_calc, perm_ses_test, boot_ses_calc ===

test_that("rank_diff output works with ses_calc", {
  set.seed(123)
  x <- rnorm(20)
  y <- rnorm(20, 0.5)

  rd <- rank_diff(x, y)
  res <- ses_calc(rd$x, rd$y, paired = TRUE, ses = "rb")
  expect_s3_class(res, "htest")
  expect_true(is.finite(unname(res$estimate)))
  expect_true(res$conf.int[1] < res$conf.int[2])
})

test_that("rank_diff output works with ses_calc hypothesis testing", {
  set.seed(123)
  x <- rnorm(20)
  y <- rnorm(20, 0.5)

  rd <- rank_diff(x, y)
  res <- ses_calc(rd$x, rd$y, paired = TRUE, ses = "cstat",
                  alternative = "two.sided", se_method = "score")
  expect_true(!is.null(res$p.value))
  expect_true(res$p.value >= 0 && res$p.value <= 1)
})

test_that("rank_diff output works with all se_methods", {
  x <- c(4.6, 4.3, 6.7, 5.8, 5.0, 4.2, 6.0, 2.0, 2.6, 10.0)
  y <- c(2.9, 2.8, 12.0, 3.8, 5.9, 6.5, 3.3, 2.3, 2.1, 14.3)

  rd <- rank_diff(x, y)
  for (method in c("score", "agresti", "fisher")) {
    res <- ses_calc(rd$x, rd$y, paired = TRUE, ses = "rb", se_method = method)
    expect_true(!is.na(unname(res$estimate)),
                info = paste("se_method =", method))
    expect_true(!is.na(res$stderr),
                info = paste("se_method =", method))
  }
})

test_that("rank_diff output works with perm_ses_test", {
  set.seed(456)
  x <- rnorm(10)
  y <- rnorm(10, 0.5)

  rd <- rank_diff(x, y)
  res <- perm_ses_test(rd$x, rd$y, paired = TRUE, ses = "rb", R = 199)
  expect_s3_class(res, "htest")
  expect_true(!is.null(res$p.value))
  expect_true(res$p.value >= 0 && res$p.value <= 1)
})

test_that("rank_diff output works with boot_ses_calc", {
  skip_on_cran()
  set.seed(789)
  x <- rnorm(15)
  y <- rnorm(15, 0.3)

  rd <- rank_diff(x, y)
  res <- suppressMessages(boot_ses_calc(rd$x, rd$y, paired = TRUE,
                                         ses = "rb", R = 99))
  expect_s3_class(res, "htest")
  expect_true(is.finite(unname(res$estimate)))
})

# === Ties handling ===

test_that("rank_diff handles ties correctly (midranks)", {
  x <- c(1, 2, 3)
  y <- c(2, 3, 4)

  rd <- rank_diff(x, y)
  # Pooled: 1, 2, 2, 3, 3, 4
  # Ranks:  1, 2.5, 2.5, 4.5, 4.5, 6
  expect_equal(rd$x, c(1, 2.5, 4.5))
  expect_equal(rd$y, c(2.5, 4.5, 6))
})
