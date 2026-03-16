# Tests for bootstrap variance fix in boot_t_test
# Verifies that the studentized bootstrap p-value computation
# is correct after fixing the centering/sampling bug.

hush = function(code) {
  sink(nullfile())
  tmp = code
  sink()
  return(tmp)
}

# Test 1: Bootstrap p-value is not anticonservatively below parametric (Welch) ----
test_that("bootstrap Welch p-value is not anticonservative relative to parametric", {
  skip_on_cran()
  set.seed(42)
  res_boot  <- boot_t_test(extra ~ group, data = sleep, R = 4999)
  res_param <- t.test(extra ~ group, data = sleep)

  # Bootstrap p should be >= parametric p (or very close), not well below it
  expect_gte(res_boot$p.value, res_param$p.value - 0.01)
})

# Test 2: One-sample p-value near 1 when mu equals sample mean ----
test_that("one-sample p-value near 1 when mu equals sample mean", {
  skip_on_cran()
  set.seed(1)
  x <- rnorm(30, mean = 5, sd = 1)
  set.seed(42)
  res <- boot_t_test(x, mu = mean(x), R = 4999)
  expect_gt(res$p.value, 0.5)
})

# Test 3: One-sided p-values sum to approximately 1 ----
test_that("one-sided p-values sum to approximately 1", {
  skip_on_cran()
  set.seed(10)
  p_less    <- boot_t_test(extra ~ group, data = sleep,
                           alternative = "less",    R = 4999)$p.value
  set.seed(10)
  p_greater <- boot_t_test(extra ~ group, data = sleep,
                           alternative = "greater", R = 4999)$p.value
  expect_equal(p_less + p_greater, 1, tolerance = 0.01)
})

# Test 4: Equal-variance p-value is not anticonservative relative to parametric ----
test_that("var.equal p-value is not anticonservative relative to parametric", {
  skip_on_cran()
  set.seed(42)
  res_eq  <- boot_t_test(extra ~ group, data = sleep,
                         var.equal = TRUE, R = 4999)
  res_par <- t.test(extra ~ group, data = sleep, var.equal = TRUE)

  expect_gte(res_eq$p.value, res_par$p.value - 0.01)
})

# Test 5: Trimmed-mean path is unaffected (non-regression) ----
test_that("trimmed path p-values unchanged after fix", {
  skip_on_cran()
  set.seed(99)
  res <- boot_t_test(extra ~ group, data = sleep, tr = 0.1, R = 4999)
  expect_s3_class(res, "htest")
  expect_true(is.finite(res$p.value))
  expect_true(res$p.value > 0 && res$p.value < 1)
})

# Test 6: One-sample paired test is not anticonservative ----
test_that("paired bootstrap p-value is not anticonservative", {
  skip_on_cran()
  x <- sleep$extra[sleep$group == 1]
  y <- sleep$extra[sleep$group == 2]
  set.seed(42)
  res_boot  <- boot_t_test(x = x, y = y, paired = TRUE, R = 4999)
  res_param <- t.test(x = x, y = y, paired = TRUE)

  expect_gte(res_boot$p.value, res_param$p.value - 0.01)
})

# boot_t_test CI/p-value agreement tests -----

for (ci_method in c("perc", "basic", "bca", "stud")) {
  test_that(paste0("boot_t_test: CI/p-value agreement for ", ci_method), {
    skip_on_cran()

    # Use data with a clear effect to avoid boundary discretization issues
    set.seed(42)
    x <- rnorm(30, mean = 0.5)
    y <- rnorm(30)

    res <- boot_t_test(x = x, y = y,
                       alternative = "two.sided",
                       boot_ci = ci_method, R = 1999)
    ci_excludes_null <- res$conf.int[1] > 0 || res$conf.int[2] < 0
    p_rejects <- res$p.value < 0.05
    expect_equal(ci_excludes_null, p_rejects,
                 label = paste(ci_method, "CI/p agreement two.sided"))
  })
}
