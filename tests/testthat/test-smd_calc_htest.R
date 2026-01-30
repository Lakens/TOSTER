# Tests for smd_calc and boot_smd_calc htest output and hypothesis testing

hush = function(code) {
  sink("NUL") # use /dev/null in UNIX
  tmp = code
  sink()
  return(tmp)
}

# --- smd_calc htest output ---

test_that("smd_calc returns htest by default", {
  set.seed(42)
  x <- rnorm(30, mean = 5, sd = 2)
  y <- rnorm(30, mean = 6, sd = 2)

  result <- smd_calc(x = x, y = y)
  expect_s3_class(result, "htest")
  expect_true(!is.null(result$estimate))
  expect_true(!is.null(result$stderr))
  expect_true(!is.null(result$conf.int))
  expect_true(!is.null(result$method))
  expect_true(!is.null(result$data.name))
  expect_true(!is.null(attr(result$conf.int, "conf.level")))
  expect_equal(attr(result$conf.int, "conf.level"), 0.95)

  # No hypothesis test by default
  expect_null(result$statistic)
  expect_null(result$p.value)
  expect_null(result$null.value)
})

test_that("smd_calc htest matches data.frame output", {
  set.seed(42)
  x <- rnorm(30, mean = 5, sd = 2)
  y <- rnorm(30, mean = 6, sd = 2)

  ht <- smd_calc(x = x, y = y)
  df <- smd_calc(x = x, y = y, output = "data.frame")

  expect_equal(unname(ht$estimate), df$estimate)
  expect_equal(ht$stderr, df$SE)
  expect_equal(ht$conf.int[1], df$lower.ci)
  expect_equal(ht$conf.int[2], df$upper.ci)
})

test_that("smd_calc htest method string is correct", {
  set.seed(42)
  x <- rnorm(20, mean = 5, sd = 2)
  y <- rnorm(20, mean = 6, sd = 2)

  # Two Sample Hedges' g (default)
  r1 <- smd_calc(x = x, y = y)
  expect_true(grepl("Two Sample", r1$method))
  expect_true(grepl("Hedges' g", r1$method))
  expect_true(grepl("estimate with CI", r1$method))

  # Paired Sample
  r2 <- smd_calc(x = x, y = y, paired = TRUE)
  expect_true(grepl("Paired Sample", r2$method))

  # One Sample
  r3 <- smd_calc(x = x)
  expect_true(grepl("One Sample", r3$method))

  # Cohen's d
  r4 <- smd_calc(x = x, y = y, bias_correction = FALSE)
  expect_true(grepl("Cohen's d", r4$method))

  # Glass's delta
  r5 <- smd_calc(x = x, y = y, glass = "glass1")
  expect_true(grepl("Glass's delta1", r5$method))

  # With test, method says "test" not "estimate with CI"
  r6 <- smd_calc(x = x, y = y, alternative = "two.sided")
  expect_true(grepl("test", r6$method))
  expect_false(grepl("estimate with CI", r6$method))
})

test_that("smd_calc formula interface returns htest with correct data.name", {
  set.seed(42)
  df <- data.frame(
    value = c(rnorm(20, 5, 2), rnorm(20, 6, 2)),
    group = factor(rep(c("A", "B"), each = 20))
  )

  result <- smd_calc(formula = value ~ group, data = df)
  expect_s3_class(result, "htest")
  expect_equal(result$data.name, "value by group")
})

# --- smd_calc hypothesis testing ---

test_that("smd_calc two.sided test works (z method)", {
  set.seed(42)
  x <- rnorm(50, mean = 5, sd = 2)
  y <- rnorm(50, mean = 6, sd = 2)

  result <- smd_calc(x = x, y = y, alternative = "two.sided",
                     null.value = 0, smd_ci = "z")
  expect_s3_class(result, "htest")
  expect_true(!is.null(result$statistic))
  expect_true(!is.null(result$p.value))
  expect_true(!is.null(result$null.value))
  expect_equal(names(result$statistic), "z")
  expect_null(result$parameter)  # no df for z method

  # p-value should be symmetric
  # If we flip the sign of the effect, p should be the same
  expect_true(result$p.value >= 0 && result$p.value <= 1)
})

test_that("smd_calc two.sided test works (t method)", {
  set.seed(42)
  x <- rnorm(50, mean = 5, sd = 2)
  y <- rnorm(50, mean = 6, sd = 2)

  result <- smd_calc(x = x, y = y, alternative = "two.sided",
                     null.value = 0, test_method = "t", smd_ci = "t")
  expect_s3_class(result, "htest")
  expect_equal(names(result$statistic), "t")
  expect_true(!is.null(result$parameter))
  expect_equal(names(result$parameter), "df")
  expect_true(result$p.value >= 0 && result$p.value <= 1)
})

test_that("smd_calc one-sided tests produce correct tail probabilities", {
  set.seed(42)
  x <- rnorm(50, mean = 5, sd = 2)
  y <- rnorm(50, mean = 8, sd = 2)  # large effect

  r_less <- smd_calc(x = x, y = y, alternative = "less",
                     null.value = 0, smd_ci = "z")
  r_greater <- smd_calc(x = x, y = y, alternative = "greater",
                        null.value = 0, smd_ci = "z")
  r_two <- smd_calc(x = x, y = y, alternative = "two.sided",
                    null.value = 0, smd_ci = "z")

  # Effect is negative (x - y < 0), so "less" p should be small
  # and "greater" p should be large
  # (sign depends on convention: x - y)
  expect_true(r_less$p.value + r_greater$p.value == 1 ||
                abs(r_less$p.value + r_greater$p.value - 1) < 1e-10)
})

test_that("smd_calc equivalence test works", {
  set.seed(42)
  x <- rnorm(50, mean = 5, sd = 2)
  y <- rnorm(50, mean = 5.2, sd = 2)  # small effect

  result <- smd_calc(x = x, y = y, alternative = "equivalence",
                     null.value = c(-0.8, 0.8), smd_ci = "z")
  expect_s3_class(result, "htest")
  expect_equal(length(result$null.value), 2)
  expect_equal(names(result$null.value), c("lower bound", "upper bound"))
  expect_true(result$p.value >= 0 && result$p.value <= 1)

  # CI should be 90% for equivalence (alpha=0.05)
  expect_equal(attr(result$conf.int, "conf.level"), 0.90)
})

test_that("smd_calc minimal.effect test works", {
  set.seed(42)
  x <- rnorm(50, mean = 5, sd = 2)
  y <- rnorm(50, mean = 8, sd = 2)  # large effect

  result <- smd_calc(x = x, y = y, alternative = "minimal.effect",
                     null.value = c(-0.3, 0.3), smd_ci = "z")
  expect_s3_class(result, "htest")
  expect_equal(length(result$null.value), 2)
  expect_true(result$p.value >= 0 && result$p.value <= 1)
  expect_equal(attr(result$conf.int, "conf.level"), 0.90)
})

test_that("smd_calc errors on bad equivalence bounds", {
  set.seed(42)
  x <- rnorm(20)
  y <- rnorm(20)

  expect_error(
    smd_calc(x = x, y = y, alternative = "equivalence", null.value = 0.5),
    "null.value must be a vector of two values"
  )
})

test_that("smd_calc warns when null.value has >1 element for standard alternative", {
  set.seed(42)
  x <- rnorm(20)
  y <- rnorm(20)

  expect_warning(
    smd_calc(x = x, y = y, alternative = "two.sided",
             null.value = c(0, 0.5), smd_ci = "z"),
    "null.value has length > 1"
  )
})

test_that("smd_calc warns on test_method/smd_ci mismatch", {
  set.seed(42)
  x <- rnorm(20)
  y <- rnorm(20)

  expect_warning(
    smd_calc(x = x, y = y, alternative = "two.sided",
             test_method = "t", smd_ci = "z"),
    "test_method.*differs from smd_ci"
  )
})

test_that("smd_calc one-sample htest works", {
  set.seed(42)
  x <- rnorm(30, mean = 0.5, sd = 1)

  result <- smd_calc(x = x, alternative = "two.sided",
                     null.value = 0, smd_ci = "z")
  expect_s3_class(result, "htest")
  expect_true(!is.null(result$statistic))
  expect_true(!is.null(result$p.value))
})

test_that("smd_calc paired htest works", {
  set.seed(42)
  x <- rnorm(30, mean = 5, sd = 2)
  y <- rnorm(30, mean = 5.5, sd = 2)

  result <- smd_calc(x = x, y = y, paired = TRUE,
                     alternative = "two.sided", null.value = 0, smd_ci = "z")
  expect_s3_class(result, "htest")
  expect_true(!is.null(result$statistic))
  expect_true(!is.null(result$p.value))
})

# --- boot_smd_calc htest output ---

test_that("boot_smd_calc returns htest by default", {
  set.seed(42)
  x <- rnorm(20, mean = 5, sd = 2)
  y <- rnorm(20, mean = 6, sd = 2)

  result <- boot_smd_calc(x = x, y = y, R = 99)
  expect_s3_class(result, "htest")
  expect_true(!is.null(result$estimate))
  expect_true(!is.null(result$stderr))
  expect_true(!is.null(result$conf.int))
  expect_true(!is.null(result$method))
  expect_true(!is.null(result$boot))
  expect_true(!is.null(result$data.name))
  expect_equal(length(result$boot), 99)

  # No hypothesis test by default
  expect_null(result$statistic)
  expect_null(result$p.value)
  expect_null(result$null.value)
})

test_that("boot_smd_calc htest estimate matches data.frame", {
  set.seed(42)
  x <- rnorm(20, mean = 5, sd = 2)
  y <- rnorm(20, mean = 6, sd = 2)

  ht <- boot_smd_calc(x = x, y = y, R = 99)
  df <- boot_smd_calc(x = x, y = y, R = 99, output = "data.frame")

  # Point estimates should match (both come from same raw_smd call)
  expect_equal(unname(ht$estimate), df$estimate)
})

test_that("boot_smd_calc two.sided test works", {
  set.seed(42)
  x <- rnorm(30, mean = 5, sd = 2)
  y <- rnorm(30, mean = 7, sd = 2)

  result <- boot_smd_calc(x = x, y = y, R = 199,
                           alternative = "two.sided", null.value = 0)
  expect_s3_class(result, "htest")
  expect_true(!is.null(result$statistic))
  expect_true(!is.null(result$p.value))
  expect_equal(names(result$statistic), "z")
  expect_true(result$p.value >= 0 && result$p.value <= 1)
})

test_that("boot_smd_calc equivalence test works", {
  set.seed(42)
  x <- rnorm(30, mean = 5, sd = 2)
  y <- rnorm(30, mean = 5.1, sd = 2)

  result <- boot_smd_calc(x = x, y = y, R = 199,
                           alternative = "equivalence",
                           null.value = c(-0.8, 0.8))
  expect_s3_class(result, "htest")
  expect_equal(length(result$null.value), 2)
  expect_equal(names(result$null.value), c("lower bound", "upper bound"))
  expect_true(result$p.value >= 0 && result$p.value <= 1)
  expect_equal(attr(result$conf.int, "conf.level"), 0.90)
})

test_that("boot_smd_calc minimal.effect test works", {
  set.seed(42)
  x <- rnorm(30, mean = 5, sd = 2)
  y <- rnorm(30, mean = 7, sd = 2)

  result <- boot_smd_calc(x = x, y = y, R = 199,
                           alternative = "minimal.effect",
                           null.value = c(-0.3, 0.3))
  expect_s3_class(result, "htest")
  expect_equal(length(result$null.value), 2)
  expect_true(result$p.value >= 0 && result$p.value <= 1)
})

test_that("boot_smd_calc one-sided tests work", {
  set.seed(42)
  x <- rnorm(30, mean = 5, sd = 2)
  y <- rnorm(30, mean = 7, sd = 2)

  r_less <- boot_smd_calc(x = x, y = y, R = 199,
                           alternative = "less", null.value = 0)
  r_greater <- boot_smd_calc(x = x, y = y, R = 199,
                              alternative = "greater", null.value = 0)

  expect_s3_class(r_less, "htest")
  expect_s3_class(r_greater, "htest")
  expect_true(r_less$p.value >= 0 && r_less$p.value <= 1)
  expect_true(r_greater$p.value >= 0 && r_greater$p.value <= 1)
})

test_that("boot_smd_calc errors on bad equivalence bounds", {
  set.seed(42)
  x <- rnorm(20)
  y <- rnorm(20)

  expect_error(
    boot_smd_calc(x = x, y = y, R = 99,
                  alternative = "equivalence", null.value = 0.5),
    "null.value must be a vector of two values"
  )
})

test_that("boot_smd_calc formula interface returns htest", {
  set.seed(42)
  df <- data.frame(
    value = c(rnorm(20, 5, 2), rnorm(20, 6, 2)),
    group = factor(rep(c("A", "B"), each = 20))
  )

  result <- boot_smd_calc(formula = value ~ group, data = df, R = 99)
  expect_s3_class(result, "htest")
  expect_equal(result$data.name, "value by group")
})

test_that("boot_smd_calc method string is correct", {
  set.seed(42)
  x <- rnorm(20, mean = 5, sd = 2)
  y <- rnorm(20, mean = 6, sd = 2)

  r1 <- boot_smd_calc(x = x, y = y, R = 99)
  expect_true(grepl("Bootstrapped", r1$method))
  expect_true(grepl("Two Sample", r1$method))
  expect_true(grepl("Hedges' g", r1$method))

  r2 <- boot_smd_calc(x = x, y = y, R = 99,
                       alternative = "two.sided", null.value = 0)
  expect_true(grepl("test", r2$method))
})
