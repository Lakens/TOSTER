# Consistency of mu handling in t_TOST, tsum_TOST, and boot_t_TOST --------
# Raw estimate, CI, and raw bounds are on the original scale;
# the SMD and its bounds are relative to mu.

hush = function(code) {
  sink(nullfile())
  tmp = code
  sink()
  return(tmp)
}

set.seed(8421)
x_one <- rnorm(40, mean = 7.4, sd = 1.2)
x_two <- rnorm(30, mean = 5.6, sd = 1)
y_two <- rnorm(35, mean = 5, sd = 1.3)
x_pair <- rnorm(25, mean = 5.5, sd = 1)
y_pair <- x_pair - rnorm(25, mean = 0.4, sd = 0.6)

# One-sample --------

test_that("one-sample t_TOST with mu reports estimate, CI, and bounds on one scale", {
  m0 <- 7.5
  res <- t_TOST(x = x_one, mu = m0, eqb = c(5.5, 8.5),
                bias_correction = FALSE)

  expect_equal(res$mu, m0)
  expect_equal(res$effsize$estimate[1], mean(x_one))
  expect_equal(c(res$effsize$lower.ci[1], res$effsize$upper.ci[1]),
               as.numeric(t.test(x_one, conf.level = 0.9)$conf.int))
  expect_true(res$effsize$lower.ci[1] <= res$effsize$estimate[1] &&
                res$effsize$estimate[1] <= res$effsize$upper.ci[1])

  # tests
  expect_equal(res$TOST$p.value[1], t.test(x_one, mu = m0)$p.value)
  expect_equal(res$TOST$p.value[2],
               t.test(x_one, mu = 5.5, alternative = "greater")$p.value)
  expect_equal(res$TOST$p.value[3],
               t.test(x_one, mu = 8.5, alternative = "less")$p.value)
  expect_equal(res$eqb$low_eq[1], 5.5)
  expect_equal(res$eqb$high_eq[1], 8.5)

  # SMD and SMD bounds are relative to mu
  expect_equal(res$effsize$estimate[2], (mean(x_one) - m0) / sd(x_one))
  expect_equal(res$eqb$low_eq[2], (5.5 - m0) / sd(x_one))
  expect_equal(res$eqb$high_eq[2], (8.5 - m0) / sd(x_one))
  expect_true(res$effsize$lower.ci[2] <= res$effsize$estimate[2] &&
                res$effsize$estimate[2] <= res$effsize$upper.ci[2])
})

test_that("smd_calc one-sample subtracts mu", {
  m0 <- 7.5
  res <- smd_calc(x = x_one, mu = m0, bias_correction = FALSE)
  expect_equal(unname(res$estimate), (mean(x_one) - m0) / sd(x_one))

  res_t <- t_TOST(x = x_one, mu = m0, eqb = c(5.5, 8.5),
                  bias_correction = FALSE)
  expect_equal(unname(res$estimate), res_t$effsize$estimate[2])
})

# Two-sample and paired --------

test_that("two-sample t_TOST with mu uses x - y - mu for the SMD", {
  m0 <- 0.5
  res <- t_TOST(x = x_two, y = y_two, mu = m0, eqb = c(-0.5, 1.5),
                var.equal = TRUE, bias_correction = FALSE)
  diff <- mean(x_two) - mean(y_two)
  sp <- sqrt(((30 - 1) * var(x_two) + (35 - 1) * var(y_two)) / (30 + 35 - 2))

  expect_equal(res$effsize$estimate[1], diff)
  expect_equal(c(res$effsize$lower.ci[1], res$effsize$upper.ci[1]),
               as.numeric(t.test(x_two, y_two, var.equal = TRUE,
                                 conf.level = 0.9)$conf.int))
  expect_equal(res$effsize$estimate[2], (diff - m0) / sp)
  expect_equal(res$eqb$low_eq[2], (-0.5 - m0) / sp)

  smd <- smd_calc(x = x_two, y = y_two, mu = m0, var.equal = TRUE,
                  bias_correction = FALSE)
  expect_equal(unname(smd$estimate), res$effsize$estimate[2])

  # SMD is ~0 when mu equals the observed difference
  res0 <- t_TOST(x = x_two, y = y_two, mu = diff, eqb = c(-1, 2))
  expect_equal(res0$effsize$estimate[2], 0, tolerance = 1e-8)
})

test_that("paired t_TOST with mu uses x - y - mu for the SMD", {
  m0 <- 0.3
  res <- t_TOST(x = x_pair, y = y_pair, paired = TRUE, mu = m0,
                eqb = c(-0.5, 1), bias_correction = FALSE)
  diff <- mean(x_pair - y_pair)

  expect_equal(res$effsize$estimate[1], diff)
  expect_equal(c(res$effsize$lower.ci[1], res$effsize$upper.ci[1]),
               as.numeric(t.test(x_pair, y_pair, paired = TRUE,
                                 conf.level = 0.9)$conf.int))
  expect_true(res$effsize$lower.ci[1] <= res$effsize$estimate[1] &&
                res$effsize$estimate[1] <= res$effsize$upper.ci[1])

  smd <- smd_calc(x = x_pair, y = y_pair, paired = TRUE, mu = m0,
                  bias_correction = FALSE)
  expect_equal(unname(smd$estimate), res$effsize$estimate[2])

  res0 <- t_TOST(x = x_pair, y = y_pair, paired = TRUE, mu = diff,
                 eqb = c(-1, 2))
  expect_equal(res0$effsize$estimate[2], 0, tolerance = 1e-8)
})

# Location-shift invariance --------

test_that("shifting data, mu, and bounds together only moves the raw estimate and CI", {
  shift <- 3
  check_shift <- function(res, res_shift) {
    expect_equal(res_shift$TOST, res$TOST)
    expect_equal(res_shift$effsize[2, ], res$effsize[2, ])
    expect_equal(res_shift$eqb[2, ], res$eqb[2, ])
    expect_equal(res_shift$effsize$estimate[1], res$effsize$estimate[1] + shift)
    expect_equal(res_shift$effsize$lower.ci[1], res$effsize$lower.ci[1] + shift)
    expect_equal(res_shift$effsize$upper.ci[1], res$effsize$upper.ci[1] + shift)
    expect_equal(res_shift$eqb$low_eq[1], res$eqb$low_eq[1] + shift)
  }

  # one-sample
  check_shift(t_TOST(x = x_one, eqb = c(-1, 1)),
              t_TOST(x = x_one + shift, mu = shift,
                     eqb = c(-1, 1) + shift))
  # two-sample
  check_shift(t_TOST(x = x_two, y = y_two, eqb = c(-1, 1)),
              t_TOST(x = x_two + shift, y = y_two, mu = shift,
                     eqb = c(-1, 1) + shift))
  # paired
  check_shift(t_TOST(x = x_pair, y = y_pair, paired = TRUE, eqb = c(-1, 1)),
              t_TOST(x = x_pair + shift, y = y_pair, paired = TRUE,
                     mu = shift, eqb = c(-1, 1) + shift))
})

# SMD bounds --------

test_that("eqbound_type = 'SMD' bounds are standardized distances from mu", {
  m0 <- 7.5
  res <- suppressWarnings(suppressMessages(
    t_TOST(x = x_one, mu = m0, eqb = c(-0.5, 0.5), eqbound_type = "SMD",
           bias_correction = FALSE)
  ))
  expect_equal(res$eqb$low_eq[1], m0 - 0.5 * sd(x_one))
  expect_equal(res$eqb$high_eq[1], m0 + 0.5 * sd(x_one))
  expect_equal(res$eqb$low_eq[2], -0.5)
  expect_equal(res$TOST$p.value[2],
               t.test(x_one, mu = m0 - 0.5 * sd(x_one),
                      alternative = "greater")$p.value)
})

# tsum_TOST --------

test_that("tsum_TOST matches t_TOST when mu is not zero", {
  compare_res <- function(a, b) {
    expect_equal(b$TOST, a$TOST)
    expect_equal(b$effsize, a$effsize)
    expect_equal(b$eqb, a$eqb)
    expect_equal(b$mu, a$mu)
  }

  # one-sample
  compare_res(
    t_TOST(x = x_one, mu = 7.5, eqb = c(5.5, 8.5), smd_ci = "z"),
    tsum_TOST(m1 = mean(x_one), sd1 = sd(x_one), n1 = length(x_one),
              mu = 7.5, eqb = c(5.5, 8.5), smd_ci = "z")
  )
  # two-sample
  compare_res(
    t_TOST(x = x_two, y = y_two, mu = 0.5, eqb = c(-0.5, 1.5), smd_ci = "z"),
    tsum_TOST(m1 = mean(x_two), sd1 = sd(x_two), n1 = length(x_two),
              m2 = mean(y_two), sd2 = sd(y_two), n2 = length(y_two),
              mu = 0.5, eqb = c(-0.5, 1.5), smd_ci = "z")
  )
  # paired
  compare_res(
    t_TOST(x = x_pair, y = y_pair, paired = TRUE, mu = 0.3,
           eqb = c(-0.5, 1), smd_ci = "z"),
    tsum_TOST(m1 = mean(x_pair), sd1 = sd(x_pair), n1 = length(x_pair),
              m2 = mean(y_pair), sd2 = sd(y_pair), n2 = length(y_pair),
              r12 = cor(x_pair, y_pair), paired = TRUE,
              mu = 0.3, eqb = c(-0.5, 1), smd_ci = "z")
  )
})

# boot_t_TOST --------

test_that("boot_t_TOST results are invariant to shifting data, mu, and bounds", {
  m0 <- 7.5
  for (ci in c("stud", "perc")) {
    set.seed(1111)
    res <- boot_t_TOST(x = x_one, mu = m0, eqb = c(5.5, 8.5),
                       R = 199, boot_ci = ci)
    set.seed(1111)
    res0 <- boot_t_TOST(x = x_one - m0, eqb = c(5.5, 8.5) - m0,
                        R = 199, boot_ci = ci)

    expect_equal(res$mu, m0)
    expect_equal(res$TOST$p.value, res0$TOST$p.value)
    expect_equal(res$effsize$estimate[1], mean(x_one))
    expect_equal(res$effsize$lower.ci[1], res0$effsize$lower.ci[1] + m0)
    expect_equal(res$effsize$upper.ci[1], res0$effsize$upper.ci[1] + m0)
    expect_equal(res$effsize$estimate[2], res0$effsize$estimate[2])
    expect_true(res$effsize$lower.ci[1] <= res$effsize$estimate[1] &&
                  res$effsize$estimate[1] <= res$effsize$upper.ci[1])
  }

  set.seed(2222)
  res <- boot_t_TOST(x = x_two, y = y_two, mu = 0.5, eqb = c(-0.5, 1.5),
                     R = 199)
  set.seed(2222)
  res0 <- boot_t_TOST(x = x_two - 0.5, y = y_two, eqb = c(-1, 1), R = 199)
  expect_equal(res$TOST$p.value, res0$TOST$p.value)
  expect_equal(res$effsize$estimate[1], mean(x_two) - mean(y_two))
})

# Messages and S3 methods --------

test_that("bounds message checks whether the interval contains mu", {
  expect_no_message(t_TOST(x = x_one, mu = 7.5, eqb = c(5.5, 8.5)))
  expect_message(t_TOST(x = x_one, eqb = c(5.5, 8.5)),
                 "does not include zero")
  expect_message(t_TOST(x = x_one, mu = 10, eqb = c(5.5, 8.5)),
                 "does not include mu")
  expect_no_message(tsum_TOST(m1 = 7.4, sd1 = 1.2, n1 = 40,
                              mu = 7.5, eqb = c(5.5, 8.5)))
})

test_that("print and describe report mu", {
  res <- t_TOST(x = x_one, mu = 7.5, eqb = c(5.5, 8.5))
  expect_output(print(res), "relative to mu = 7.5")
  expect_output(print(res), "Equivalence Bounds: Raw \\[5.5, 8.5\\]")
  expect_true(grepl("7.5", describe(res)))

  res_sum <- tsum_TOST(m1 = mean(x_one), sd1 = sd(x_one), n1 = length(x_one),
                       mu = 7.5, eqb = c(5.5, 8.5))
  expect_true(grepl("7.5", describe(res_sum)))

  # no mu note when mu is zero; older objects without mu still work
  res_zero <- t_TOST(x = x_two, y = y_two, eqb = 1)
  expect_false(any(grepl("relative to mu", capture.output(print(res_zero)))))
  res_zero$mu <- NULL
  expect_type(hush(describe(res_zero)), "character")
})
