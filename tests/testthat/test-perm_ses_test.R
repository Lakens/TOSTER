# Tests for perm_ses_test

hush = function(code) {
  sink("NUL") # use /dev/null in UNIX
  tmp = code
  sink()
  return(tmp)
}

# === Basic two-sample tests ===

test_that("perm_ses_test: two-sample returns valid htest for all ses types", {
  set.seed(300)
  x <- rnorm(8)
  y <- rnorm(8, mean = 0.5)

  for (ses_type in c("rb", "cstat", "odds", "logodds")) {
    res <- hush(perm_ses_test(x = x, y = y, ses = ses_type, R = 199))
    expect_s3_class(res, "htest")
    expect_true(!is.null(res$p.value))
    expect_true(res$p.value >= 0 && res$p.value <= 1)
    expect_true(!is.null(res$estimate))
    expect_true(!is.null(res$conf.int))
    expect_true(!is.null(res$statistic))
    expect_true(!is.null(res$parameter))
    expect_equal(res$alternative, "two.sided")
  }
})

test_that("perm_ses_test: point estimates match ses_calc", {
  set.seed(301)
  x <- rnorm(10)
  y <- rnorm(10, mean = 0.5)

  for (ses_type in c("rb", "cstat", "odds", "logodds")) {
    res_perm <- hush(perm_ses_test(x = x, y = y, ses = ses_type, R = 99))
    res_ses <- ses_calc(x = x, y = y, ses = ses_type, output = "data.frame")

    expect_equal(unname(res_perm$estimate), res_ses$estimate,
                 info = paste("ses =", ses_type))
  }
})

# === Paired samples ===

test_that("perm_ses_test: paired samples", {
  set.seed(302)
  x <- rnorm(10)
  y <- x + rnorm(10, mean = 0.3, sd = 1)

  res <- hush(perm_ses_test(x = x, y = y, paired = TRUE, ses = "rb", R = 199))
  expect_s3_class(res, "htest")
  expect_true(is.numeric(unname(res$estimate)))
  expect_true(grepl("Paired", res$method))
})

test_that("perm_ses_test: paired point estimate matches ses_calc", {
  set.seed(303)
  x <- rnorm(8)
  y <- x + rnorm(8, mean = 0.5, sd = 1)

  res_perm <- hush(perm_ses_test(x = x, y = y, paired = TRUE, ses = "rb", R = 99))
  res_ses <- ses_calc(x = x, y = y, paired = TRUE, ses = "rb", output = "data.frame")

  expect_equal(unname(res_perm$estimate), res_ses$estimate)
})

# === One-sample ===

test_that("perm_ses_test: one-sample", {
  set.seed(304)
  x <- rnorm(10, mean = 1)

  res <- hush(perm_ses_test(x = x, ses = "rb", R = 199))
  expect_s3_class(res, "htest")
  expect_true(is.numeric(unname(res$estimate)))
  expect_true(grepl("One-Sample", res$method))
})

test_that("perm_ses_test: one-sample point estimate matches ses_calc", {
  set.seed(305)
  x <- rnorm(8, mean = 0.5)

  res_perm <- hush(perm_ses_test(x = x, ses = "cstat", R = 99))
  res_ses <- ses_calc(x = x, ses = "cstat", output = "data.frame")

  expect_equal(unname(res_perm$estimate), res_ses$estimate)
})

# === One-sided alternatives ===

test_that("perm_ses_test: one-sided alternatives", {
  set.seed(306)
  x <- rnorm(12, mean = 2)
  y <- rnorm(12)

  res_greater <- hush(perm_ses_test(x = x, y = y, ses = "rb",
                                     alternative = "greater", R = 499))
  res_less <- hush(perm_ses_test(x = x, y = y, ses = "rb",
                                  alternative = "less", R = 499))

  # With x >> y, rb should be positive
  expect_true(unname(res_greater$estimate) > 0)
  # p-value for greater should be smaller than for less
  expect_true(res_greater$p.value < res_less$p.value)
})

# === Equivalence testing ===

test_that("perm_ses_test: equivalence testing", {
  set.seed(307)
  x <- rnorm(15)
  y <- rnorm(15, mean = 0.1)

  res <- hush(perm_ses_test(x = x, y = y, ses = "rb",
                             alternative = "equivalence",
                             mu = c(-0.5, 0.5),
                             R = 499))
  expect_s3_class(res, "htest")
  expect_equal(res$alternative, "equivalence")
  expect_equal(length(res$null.value), 2)
  expect_true(res$p.value >= 0 && res$p.value <= 1)
})

test_that("perm_ses_test: equivalence requires two bounds", {
  set.seed(308)
  x <- rnorm(10)
  y <- rnorm(10)

  expect_error(
    perm_ses_test(x = x, y = y, alternative = "equivalence", mu = 0),
    "mu.*must specify two bounds"
  )
})

# === Minimal effect testing ===

test_that("perm_ses_test: minimal effect testing", {
  set.seed(309)
  x <- rnorm(15, mean = 2)
  y <- rnorm(15)

  res <- hush(perm_ses_test(x = x, y = y, ses = "rb",
                             alternative = "minimal.effect",
                             mu = c(-0.3, 0.3),
                             R = 499))
  expect_s3_class(res, "htest")
  expect_equal(res$alternative, "minimal.effect")
  expect_equal(length(res$null.value), 2)
  expect_true(res$p.value >= 0 && res$p.value <= 1)
})

# === Formula interface ===

test_that("perm_ses_test: formula interface", {
  set.seed(310)
  df <- data.frame(
    val = c(rnorm(10), rnorm(10, mean = 0.5)),
    grp = factor(rep(c("A", "B"), each = 10))
  )
  res <- hush(perm_ses_test(formula = val ~ grp, data = df,
                             ses = "rb", R = 199))
  expect_s3_class(res, "htest")
  expect_equal(res$data.name, "val by grp")
})

test_that("perm_ses_test: formula with paired warns", {
  set.seed(311)
  df <- data.frame(
    val = c(rnorm(10), rnorm(10, mean = 0.5)),
    grp = factor(rep(c("A", "B"), each = 10))
  )
  expect_message(
    hush(perm_ses_test(formula = val ~ grp, data = df,
                        paired = TRUE, ses = "rb", R = 199)),
    "paired = TRUE.*formula interface.*not recommended"
  )
})

# === Exact vs randomization ===

test_that("perm_ses_test: exact permutation for small samples", {
  x <- c(1, 2, 3, 4, 5)
  y <- c(3, 4, 5, 6, 7)

  res <- hush(perm_ses_test(x = x, y = y, ses = "rb"))
  expect_s3_class(res, "htest")
  expect_true(grepl("Exact", res$method))
  # R.used should be choose(10, 5) = 252
  expect_equal(res$R.used, choose(10, 5))
})

test_that("perm_ses_test: randomization when R specified", {
  set.seed(312)
  x <- rnorm(20)
  y <- rnorm(20, mean = 0.5)

  res <- hush(perm_ses_test(x = x, y = y, ses = "rb", R = 199))
  expect_s3_class(res, "htest")
  expect_true(grepl("Randomization", res$method))
})

# === Complete separation ===

test_that("perm_ses_test: complete separation on rb/cstat scale works", {
  x <- c(10, 11, 12)
  y <- c(1, 2, 3)

  res_rb <- hush(perm_ses_test(x = x, y = y, ses = "rb"))
  expect_s3_class(res_rb, "htest")
  expect_equal(unname(res_rb$estimate), 1)

  res_cs <- hush(perm_ses_test(x = x, y = y, ses = "cstat"))
  expect_equal(unname(res_cs$estimate), 1)
})

test_that("perm_ses_test: complete separation warns on odds scale", {
  x <- c(10, 11, 12)
  y <- c(1, 2, 3)

  expect_warning(
    hush(perm_ses_test(x = x, y = y, ses = "odds")),
    "Complete separation detected"
  )
})

# === Input validation ===

test_that("perm_ses_test: invalid alpha rejected", {
  x <- rnorm(10)
  expect_error(perm_ses_test(x = x, alpha = 0),
               "alpha must be a numeric value between 0 and 1")
  expect_error(perm_ses_test(x = x, alpha = 1),
               "alpha must be a numeric value between 0 and 1")
  expect_error(perm_ses_test(x = x, alpha = "a"),
               "alpha must be a numeric value between 0 and 1")
})

test_that("perm_ses_test: invalid R rejected", {
  x <- rnorm(10)
  expect_error(perm_ses_test(x = x, R = -1),
               "'R' must be NULL")
  expect_error(perm_ses_test(x = x, R = 0),
               "'R' must be NULL")
})

# === Permutation distribution stored ===

test_that("perm_ses_test: perm.dist stored when keep_perm = TRUE", {
  set.seed(313)
  x <- rnorm(8)
  y <- rnorm(8, mean = 0.5)

  res_keep <- hush(perm_ses_test(x = x, y = y, ses = "rb",
                                  R = 99, keep_perm = TRUE))
  expect_true(!is.null(res_keep$perm.dist))
  expect_equal(length(res_keep$perm.dist), res_keep$R.used)

  res_no <- hush(perm_ses_test(x = x, y = y, ses = "rb",
                                R = 99, keep_perm = FALSE))
  expect_null(res_no$perm.dist)
})

# === Confidence intervals ===

test_that("perm_ses_test: CI bounds are ordered", {
  set.seed(314)
  x <- rnorm(10)
  y <- rnorm(10, mean = 0.5)

  for (alt in c("two.sided", "equivalence", "minimal.effect")) {
    mu_val <- if (alt %in% c("equivalence", "minimal.effect")) c(-0.5, 0.5) else 0
    res <- hush(perm_ses_test(x = x, y = y, ses = "rb",
                               alternative = alt, mu = mu_val, R = 199))
    expect_true(res$conf.int[1] <= res$conf.int[2], info = paste("alt =", alt))
  }
})

test_that("perm_ses_test: CI within valid range for rb", {
  set.seed(315)
  x <- rnorm(10)
  y <- rnorm(10, mean = 0.5)

  res <- hush(perm_ses_test(x = x, y = y, ses = "rb", R = 199))
  expect_true(res$conf.int[1] >= -1)
  expect_true(res$conf.int[2] <= 1)
})

test_that("perm_ses_test: CI within valid range for cstat", {
  set.seed(316)
  x <- rnorm(10)
  y <- rnorm(10, mean = 0.5)

  res <- hush(perm_ses_test(x = x, y = y, ses = "cstat", R = 199))
  expect_true(res$conf.int[1] >= 0)
  expect_true(res$conf.int[2] <= 1)
})

# === Different ses scales produce consistent results ===

test_that("perm_ses_test: effect size transformations are consistent", {
  set.seed(317)
  x <- rnorm(8)
  y <- rnorm(8, mean = 0.5)

  res_rb <- hush(perm_ses_test(x = x, y = y, ses = "rb", R = 199))
  res_cs <- hush(perm_ses_test(x = x, y = y, ses = "cstat", R = 199))
  res_odds <- hush(perm_ses_test(x = x, y = y, ses = "odds", R = 199))
  res_lo <- hush(perm_ses_test(x = x, y = y, ses = "logodds", R = 199))

  # rb and cstat should be related: cstat = (rb + 1) / 2
  expect_equal(unname(res_cs$estimate), (unname(res_rb$estimate) + 1) / 2)

  # odds = cstat / (1 - cstat)
  expect_equal(unname(res_odds$estimate),
               unname(res_cs$estimate) / (1 - unname(res_cs$estimate)))

  # logodds = log(odds)
  expect_equal(unname(res_lo$estimate), log(unname(res_odds$estimate)))
})

# === Non-symmetric two-sided test ===

test_that("perm_ses_test: non-symmetric two-sided test works", {
  set.seed(318)
  x <- rnorm(10)
  y <- rnorm(10, mean = 0.5)

  res_sym <- hush(perm_ses_test(x = x, y = y, ses = "rb",
                                 symmetric = TRUE, R = 499))
  res_asym <- hush(perm_ses_test(x = x, y = y, ses = "rb",
                                  symmetric = FALSE, R = 499))

  # Both should return valid p-values (may differ slightly)
  expect_true(res_sym$p.value >= 0 && res_sym$p.value <= 1)
  expect_true(res_asym$p.value >= 0 && res_asym$p.value <= 1)
})

# === Exact paired permutation ===

test_that("perm_ses_test: exact paired permutation", {
  x <- c(1, 2, 3, 4, 5)
  y <- c(2, 3, 4, 5, 6)

  res <- hush(perm_ses_test(x = x, y = y, paired = TRUE, ses = "rb"))
  expect_s3_class(res, "htest")
  expect_true(grepl("Exact", res$method))
  expect_equal(res$R.used, 2^5)
})
