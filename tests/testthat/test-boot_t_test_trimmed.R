# Tests for trimmed means support in boot_t_test

hush = function(code) {
  sink(nullfile())
  on.exit(sink())
  invisible(force(suppressMessages(code)))
}

# --- 10.1 Input validation ---
test_that("tr validation works", {
  x <- rnorm(20)

  # Invalid tr values
  expect_error(hush(boot_t_test(x, tr = -0.1)), "'tr' must be")
  expect_error(hush(boot_t_test(x, tr = 0.5)), "'tr' must be")
  expect_error(hush(boot_t_test(x, tr = 1.0)), "'tr' must be")
  expect_error(hush(boot_t_test(x, tr = c(0.1, 0.2))), "'tr' must be")
  expect_error(hush(boot_t_test(x, tr = NA)), "'tr' must be")

  # Sample too small for trimming
  expect_error(hush(boot_t_test(rnorm(3), tr = 0.4)), "Sample size too small")
})

# --- 10.2 Backward compatibility (tr = 0) ---
test_that("tr = 0 matches existing behavior", {
  data(sleep)

  set.seed(100)
  res_default <- hush(boot_t_test(extra ~ group, data = sleep, R = 599))

  set.seed(100)
  res_tr0 <- hush(boot_t_test(extra ~ group, data = sleep, tr = 0, R = 599))

  expect_equal(res_default$statistic, res_tr0$statistic)
  expect_equal(res_default$p.value, res_tr0$p.value)
  expect_equal(res_default$conf.int, res_tr0$conf.int)
  expect_equal(res_default$estimate, res_tr0$estimate)
})

# --- 10.3 One-sample trimmed bootstrap ---
test_that("one-sample trimmed bootstrap works", {
  set.seed(42)
  x <- c(rnorm(18), 10, -10)  # data with outliers

  set.seed(100)
  res <- hush(boot_t_test(x, tr = 0.1, R = 599))

  expect_s3_class(res, "htest")
  expect_true(grepl("Yuen", res$method))
  expect_true(grepl("trimmed mean", names(res$estimate)))
  expect_equal(unname(res$estimate), mean(x, trim = 0.1), tolerance = 1e-10)
})

# --- 10.4 Paired trimmed bootstrap ---
test_that("paired trimmed bootstrap works", {
  set.seed(42)
  x <- rnorm(20, mean = 5)
  y <- rnorm(20, mean = 4.5)

  set.seed(100)
  res <- hush(boot_t_test(x = x, y = y, paired = TRUE, tr = 0.2, R = 599))

  expect_s3_class(res, "htest")
  expect_true(grepl("Paired Yuen", res$method))
  expect_true(grepl("trimmed mean of the differences", names(res$estimate)))
})

# --- 10.5 Two-sample trimmed bootstrap (Welch and equal variance) ---
test_that("two-sample Welch trimmed bootstrap works", {
  set.seed(42)
  x <- c(rnorm(20, mean = 1), 15)   # outlier in x
  y <- c(rnorm(25, mean = 0), -12)   # outlier in y

  set.seed(100)
  res <- hush(boot_t_test(x = x, y = y, tr = 0.1, R = 599))

  expect_s3_class(res, "htest")
  expect_true(grepl("Welch Yuen", res$method))
  expect_length(res$estimate, 2)
  expect_true(all(grepl("trimmed mean", names(res$estimate))))
})

test_that("two-sample equal variance trimmed bootstrap works", {
  set.seed(42)
  x <- rnorm(20, mean = 1)
  y <- rnorm(20, mean = 0)

  set.seed(100)
  res <- hush(boot_t_test(x = x, y = y, var.equal = TRUE, tr = 0.1, R = 599))

  expect_s3_class(res, "htest")
  expect_true(grepl("Yuen", res$method))
  expect_false(grepl("Welch", res$method))
})

# --- 10.6 All CI methods with trimming ---
test_that("all boot_ci methods work with trimming", {
  set.seed(42)
  x <- rnorm(30)

  for (ci_method in c("stud", "basic", "perc", "bca")) {
    set.seed(100)
    res <- hush(boot_t_test(x, tr = 0.1, boot_ci = ci_method, R = 599))

    expect_s3_class(res, "htest")
    expect_true(!is.na(res$conf.int[1]) && !is.na(res$conf.int[2]),
                info = paste("CI method:", ci_method))
    expect_true(res$conf.int[1] < res$conf.int[2],
                info = paste("CI method:", ci_method))
  }
})

# --- 10.7 Equivalence and minimal effect with trimming ---
test_that("equivalence testing works with trimming", {
  set.seed(42)
  x <- rnorm(30, mean = 0.1)
  y <- rnorm(30, mean = 0)

  set.seed(100)
  res <- hush(boot_t_test(x = x, y = y, alternative = "equivalence",
                           mu = c(-1, 1), tr = 0.2, R = 599))

  expect_s3_class(res, "htest")
  expect_equal(res$alternative, "equivalence")
  expect_true(grepl("Yuen", res$method))
})

test_that("minimal.effect testing works with trimming", {
  set.seed(42)
  x <- rnorm(30, mean = 5)
  y <- rnorm(30, mean = 0)

  set.seed(100)
  res <- hush(boot_t_test(x = x, y = y, alternative = "minimal.effect",
                           mu = c(-1, 1), tr = 0.2, R = 599))

  expect_s3_class(res, "htest")
  expect_equal(res$alternative, "minimal.effect")
})

# --- 10.8 Formula interface passes tr through ---
test_that("formula interface passes tr correctly", {
  data(sleep)

  set.seed(100)
  res_formula <- hush(boot_t_test(extra ~ group, data = sleep, tr = 0.1, R = 599))

  set.seed(100)
  res_xy <- hush(boot_t_test(x = sleep$extra[sleep$group == 1],
                              y = sleep$extra[sleep$group == 2],
                              tr = 0.1, R = 599))

  expect_equal(res_formula$statistic, res_xy$statistic, tolerance = 1e-10)
  expect_equal(res_formula$p.value, res_xy$p.value, tolerance = 1e-10)
})

# --- 10.9 Trimming reduces outlier influence ---
test_that("trimming reduces outlier influence on CIs", {
  set.seed(42)
  x_clean <- rnorm(20)
  x_contaminated <- c(x_clean, 50, -50)  # extreme outliers

  set.seed(100)
  res_no_trim <- hush(boot_t_test(x_contaminated, R = 999))

  set.seed(100)
  res_trimmed <- hush(boot_t_test(x_contaminated, tr = 0.1, R = 999))

  # Trimmed CI should be narrower since outliers are downweighted
  width_no_trim <- diff(res_no_trim$conf.int)
  width_trimmed <- diff(res_trimmed$conf.int)
  expect_true(width_trimmed < width_no_trim)
})

# --- 10.10 Consistency with perm_t_test estimates ---
test_that("trimmed estimates match perm_t_test", {
  set.seed(42)
  x <- rnorm(25)
  y <- rnorm(25)

  set.seed(100)
  res_boot <- hush(boot_t_test(x = x, y = y, tr = 0.2, R = 599))

  set.seed(100)
  res_perm <- hush(perm_t_test(x = x, y = y, tr = 0.2, R = 599))

  # Point estimates should be identical (same trimmed means)
  expect_equal(unname(res_boot$estimate), unname(res_perm$estimate))

  # Observed t-statistic and df should be identical
  expect_equal(unname(res_boot$statistic), unname(res_perm$statistic),
               tolerance = 1e-10)
  expect_equal(unname(res_boot$parameter), unname(res_perm$parameter),
               tolerance = 1e-10)
})

# --- 10.11 Edge cases ---
test_that("trimming edge cases are handled", {
  # Very small tr (effectively no trimming for small n)
  set.seed(42)
  x <- rnorm(10)

  set.seed(100)
  res_tiny <- hush(boot_t_test(x, tr = 0.01, R = 599))
  expect_s3_class(res_tiny, "htest")

  # Larger trimming
  set.seed(100)
  res_large <- hush(boot_t_test(rnorm(30), tr = 0.4, R = 599))
  expect_s3_class(res_large, "htest")

  # NA handling with trimming
  x_na <- c(rnorm(20), NA, NA)
  set.seed(100)
  res_na <- hush(boot_t_test(x_na, tr = 0.1, R = 599))
  expect_s3_class(res_na, "htest")
})

# --- 10.12 Output structure completeness ---
test_that("trimmed boot_t_test returns complete htest object", {
  set.seed(42)
  x <- rnorm(20)

  set.seed(100)
  res <- hush(boot_t_test(x, tr = 0.2, R = 599))

  expect_true("statistic" %in% names(res))
  expect_true("parameter" %in% names(res))
  expect_true("p.value" %in% names(res))
  expect_true("stderr" %in% names(res))
  expect_true("conf.int" %in% names(res))
  expect_true("estimate" %in% names(res))
  expect_true("null.value" %in% names(res))
  expect_true("alternative" %in% names(res))
  expect_true("method" %in% names(res))
  expect_true("boot" %in% names(res))
  expect_true("data.name" %in% names(res))
  expect_true("call" %in% names(res))

  # p-value in valid range
  expect_true(res$p.value >= 0 && res$p.value <= 1)

  # boot vector has correct length
  expect_length(res$boot, 599)
})
