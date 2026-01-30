# Tests for ses_calc, boot_ses_calc, and helper functions

hush = function(code) {
  sink("NUL") # use /dev/null in UNIX
  tmp = code
  sink()
  return(tmp)
}

# === Helper function tests ===

test_that("rbs_calc: two-sample symmetry", {
  set.seed(8492)
  x <- rnorm(30)
  y <- rnorm(30, mean = 0.5)

  r1 <- TOSTER:::rbs_calc(x = x, y = y, mu = 0, paired = FALSE)
  r2 <- TOSTER:::rbs_calc(x = y, y = x, mu = 0, paired = FALSE)
  expect_equal(abs(r1), abs(r2))
  expect_equal(r1, -r2)
})

test_that("rbs_calc: paired symmetry", {
  set.seed(8493)
  x <- rnorm(30)
  y <- rnorm(30)

  r1 <- TOSTER:::rbs_calc(x = x, y = y, mu = 0, paired = TRUE)
  r2 <- TOSTER:::rbs_calc(x = y, y = x, mu = 0, paired = TRUE)
  expect_equal(abs(r1), abs(r2))
  expect_equal(r1, -r2)
})

test_that("rbs_calc: one-sample with mu shift", {
  set.seed(8494)
  x <- rnorm(30, mean = 2)

  # rbs_calc paired convention: z = (x - y) - mu
  # When x=zeros, y=positive_data: z is negative, giving negative rho
  # ses_calc swaps args (passes y,x) so end-user sees positive rb for positive mean
  # Here we test the raw internal convention
  r <- TOSTER:::rbs_calc(x = rep(0, length(x)), y = x, mu = 0, paired = TRUE)
  expect_true(r < 0)
})

test_that("Transformation round-trips", {
  # rb <-> cstat
  expect_equal(TOSTER:::cstat_to_rb(TOSTER:::rb_to_cstat(0.6)), 0.6)
  expect_equal(TOSTER:::cstat_to_rb(TOSTER:::rb_to_cstat(-0.3)), -0.3)
  expect_equal(TOSTER:::rb_to_cstat(0), 0.5)
  expect_equal(TOSTER:::rb_to_cstat(1), 1)
  expect_equal(TOSTER:::rb_to_cstat(-1), 0)

  # odds <-> pr
  expect_equal(TOSTER:::pr_to_odds(TOSTER:::odds_to_pr(2.5)), 2.5)
  expect_equal(TOSTER:::odds_to_pr(TOSTER:::pr_to_odds(0.7)), 0.7)

  # odds via log
  expect_equal(TOSTER:::pr_to_odds(0.5, log = TRUE), 0)  # log-odds of 0.5 = 0
  expect_equal(TOSTER:::odds_to_pr(0, log = TRUE), 0.5)  # p from log-odds 0 = 0.5

  # rho <-> z (Fisher)
  expect_equal(TOSTER:::z_to_rho(TOSTER:::rho_to_z(0.45)), 0.45)
})

test_that("to_logodds transformations are correct", {
  # rb = 0 -> cstat = 0.5 -> logodds = 0

  expect_equal(TOSTER:::to_logodds(0, "rb"), 0)

  # cstat = 0.5 -> logodds = 0
  expect_equal(TOSTER:::to_logodds(0.5, "cstat"), 0)

  # odds = 1 -> logodds = 0
  expect_equal(TOSTER:::to_logodds(1, "odds"), 0)

  # logodds pass-through
  expect_equal(TOSTER:::to_logodds(1.5, "logodds"), 1.5)

  # Consistency across scales for same effect
  rb_val <- 0.4
  cstat_val <- TOSTER:::rb_to_cstat(rb_val)
  odds_val <- TOSTER:::rb_to_odds(rb_val)
  logodds_val <- log(odds_val)

  expect_equal(TOSTER:::to_logodds(rb_val, "rb"),
               TOSTER:::to_logodds(cstat_val, "cstat"))
  expect_equal(TOSTER:::to_logodds(rb_val, "rb"),
               TOSTER:::to_logodds(odds_val, "odds"))
  expect_equal(TOSTER:::to_logodds(rb_val, "rb"),
               TOSTER:::to_logodds(logodds_val, "logodds"))
})

test_that("compute_placements returns expected structure", {
  x <- c(1, 3, 5, 7)
  y <- c(2, 4, 6, 8)
  pl <- TOSTER:::compute_placements(x, y)

  expect_true(is.list(pl))
  expect_equal(pl$n1, 4)
  expect_equal(pl$n2, 4)
  expect_equal(length(pl$V), 4)
  expect_equal(length(pl$W), 4)
  expect_true(pl$p_hat >= 0 && pl$p_hat <= 1)
})

test_that("compute_placements: no overlap gives p=0 or p=1", {
  x <- c(1, 2, 3)
  y <- c(10, 11, 12)

  pl <- TOSTER:::compute_placements(x, y)
  expect_equal(pl$p_hat, 0)

  pl2 <- TOSTER:::compute_placements(y, x)
  expect_equal(pl2$p_hat, 1)
})

test_that("var_concordance_twosample is non-negative", {
  set.seed(123)
  x <- rnorm(20)
  y <- rnorm(20, mean = 1)
  pl <- TOSTER:::compute_placements(x, y)
  v <- TOSTER:::var_concordance_twosample(pl)

  expect_true(v >= 0)
})

test_that("var_concordance_paired handles small samples", {
  d <- c(-2, -1, 1, 2)
  v <- TOSTER:::var_concordance_paired(d)
  expect_true(is.numeric(v) && v >= 0)

  # Too few non-zero differences
  v_small <- TOSTER:::var_concordance_paired(c(1, 0, 0, 0))
  expect_true(is.na(v_small))  # Only 1 non-zero
})

test_that("ses_compute_agresti returns expected structure", {
  set.seed(444)
  x <- rnorm(20)
  y <- rnorm(20, mean = 0.5)

  res <- TOSTER:::ses_compute_agresti(x, y, paired = FALSE)
  expect_true(is.list(res))
  expect_true(all(c("cstat", "rb", "odds", "logodds",
                     "se_cstat", "se_rb", "se_odds", "se_logodds",
                     "paired", "boundary_corrected") %in% names(res)))
  expect_false(res$boundary_corrected)
  expect_true(res$se_rb > 0)
})

test_that("ses_compute_agresti: boundary correction applied", {
  # Complete separation
  x <- c(10, 11, 12)
  y <- c(1, 2, 3)
  res <- TOSTER:::ses_compute_agresti(x, y, paired = FALSE)
  expect_true(res$boundary_corrected)
  # cstat should be near 1 but corrected: 1 - 0.5/(n1*n2) = 1 - 0.5/9 ≈ 0.944
  expect_true(res$cstat < 1)
  expect_true(res$cstat > 0.9)
})

test_that("ses_ci_logodds returns valid CIs", {
  set.seed(555)
  x <- rnorm(30)
  y <- rnorm(30, mean = 0.5)
  est <- TOSTER:::ses_compute_agresti(x, y, paired = FALSE)
  ci_res <- TOSTER:::ses_ci_logodds(est, conf.level = 0.95)

  expect_true(ci_res$ci_rb[1] < ci_res$ci_rb[2])
  expect_true(ci_res$ci_cstat[1] < ci_res$ci_cstat[2])
  expect_true(ci_res$ci_odds[1] < ci_res$ci_odds[2])
  expect_true(ci_res$ci_logodds[1] < ci_res$ci_logodds[2])

  # rb CI within [-1, 1]
  expect_true(ci_res$ci_rb[1] >= -1)
  expect_true(ci_res$ci_rb[2] <= 1)

  # cstat CI within [0, 1]
  expect_true(ci_res$ci_cstat[1] >= 0)
  expect_true(ci_res$ci_cstat[2] <= 1)

  # odds CI positive
  expect_true(ci_res$ci_odds[1] > 0)
})

test_that("ses_ci_fisher returns valid CIs", {
  set.seed(666)
  x <- rnorm(25)
  y <- rnorm(25, mean = 0.3)
  r_rbs <- TOSTER:::rbs_calc(x, y, mu = 0, paired = FALSE)

  ci_res <- TOSTER:::ses_ci_fisher(r_rbs, n1 = 25, n2 = 25,
                                    paired = FALSE, conf.level = 0.95)
  expect_true(ci_res$ci_rb[1] < ci_res$ci_rb[2])
  expect_true(ci_res$se_rb > 0)
  expect_true(ci_res$se_cstat > 0)
})

test_that("se_fisher_z: two-sample vs paired", {
  se_two <- TOSTER:::se_fisher_z(20, 20, paired = FALSE)
  se_paired <- TOSTER:::se_fisher_z(20, paired = TRUE)
  expect_true(is.numeric(se_two) && se_two > 0)
  expect_true(is.numeric(se_paired) && se_paired > 0)
})


# === ses_calc tests ===

test_that("ses_calc: default output is htest with no hypothesis test", {
  set.seed(100)
  x <- rnorm(20)
  y <- rnorm(20, mean = 0.5)

  res <- ses_calc(x = x, y = y)
  expect_s3_class(res, "htest")
  expect_equal(res$alternative, "none")
  # No test components when alternative = "none"
  expect_null(res$statistic)
  expect_null(res$p.value)
  expect_null(res$null.value)
  # Should have estimate and CI

  expect_true(!is.null(res$estimate))
  expect_true(!is.null(res$conf.int))
  expect_true(!is.null(res$stderr))
})

test_that("ses_calc: all effect size types return valid estimates", {
  set.seed(101)
  x <- rnorm(25)
  y <- rnorm(25, mean = 0.5)

  for (ses_type in c("rb", "cstat", "odds", "logodds")) {
    res <- ses_calc(x = x, y = y, ses = ses_type)
    expect_s3_class(res, "htest")
    expect_true(is.numeric(res$estimate))
  }

  # rb should be in [-1, 1]
  res_rb <- ses_calc(x = x, y = y, ses = "rb")
  expect_true(unname(res_rb$estimate) >= -1 && unname(res_rb$estimate) <= 1)

  # cstat should be in [0, 1]
  res_cs <- ses_calc(x = x, y = y, ses = "cstat")
  expect_true(unname(res_cs$estimate) >= 0 && unname(res_cs$estimate) <= 1)

  # odds should be positive
  res_odds <- ses_calc(x = x, y = y, ses = "odds")
  expect_true(unname(res_odds$estimate) > 0)
})

test_that("ses_calc: effect size transformations are consistent", {
  set.seed(102)
  x <- rnorm(30)
  y <- rnorm(30, mean = 0.3)

  res_rb <- ses_calc(x = x, y = y, ses = "rb", output = "data.frame")
  res_cs <- ses_calc(x = x, y = y, ses = "cstat", output = "data.frame")
  res_odds <- ses_calc(x = x, y = y, ses = "odds", output = "data.frame")
  res_lo <- ses_calc(x = x, y = y, ses = "logodds", output = "data.frame")

  # rb and cstat should be related: cstat = (rb + 1) / 2
  expect_equal(res_cs$estimate, (res_rb$estimate + 1) / 2)

  # odds = cstat / (1 - cstat)
  expect_equal(res_odds$estimate, res_cs$estimate / (1 - res_cs$estimate))

  # logodds = log(odds)
  expect_equal(res_lo$estimate, log(res_odds$estimate))
})

test_that("ses_calc: hypothesis test components present when alternative != 'none'", {
  set.seed(103)
  x <- rnorm(20)
  y <- rnorm(20, mean = 0.8)

  res <- ses_calc(x = x, y = y, alternative = "two.sided")
  expect_s3_class(res, "htest")
  expect_true(!is.null(res$statistic))
  expect_true(!is.null(res$p.value))
  expect_true(!is.null(res$null.value))
  expect_true(res$p.value >= 0 && res$p.value <= 1)
})

test_that("ses_calc: one-sided alternatives", {
  set.seed(104)
  x <- rnorm(25, mean = 2)  # Clearly larger
  y <- rnorm(25)

  res_greater <- ses_calc(x = x, y = y, ses = "rb", alternative = "greater")
  res_less <- ses_calc(x = x, y = y, ses = "rb", alternative = "less")

  # With x >> y, rb should be positive
  expect_true(unname(res_greater$estimate) > 0)
  # p-value for greater should be small, for less should be large
  expect_true(res_greater$p.value < 0.05)
  expect_true(res_less$p.value > 0.5)
})

test_that("ses_calc: equivalence and minimal effect testing", {
  set.seed(105)
  x <- rnorm(40)
  y <- rnorm(40, mean = 0.1)

  res_eq <- ses_calc(x = x, y = y, ses = "rb",
                     alternative = "equivalence",
                     null.value = c(-0.5, 0.5))
  expect_true(!is.null(res_eq$p.value))
  expect_equal(length(res_eq$null.value), 2)

  res_met <- ses_calc(x = x, y = y, ses = "rb",
                      alternative = "minimal.effect",
                      null.value = c(-0.5, 0.5))
  expect_true(!is.null(res_met$p.value))
  expect_equal(length(res_met$null.value), 2)
})

test_that("ses_calc: equivalence requires two bounds", {
  set.seed(106)
  x <- rnorm(20)
  y <- rnorm(20)

  expect_error(
    ses_calc(x = x, y = y, alternative = "equivalence", null.value = 0),
    "null.value must be a vector of two values"
  )
})

test_that("ses_calc: data.frame output format", {
  set.seed(107)
  x <- rnorm(20)
  y <- rnorm(20, mean = 0.5)

  res <- ses_calc(x = x, y = y, ses = "rb", output = "data.frame")
  expect_true(is.data.frame(res))
  expect_true(all(c("estimate", "SE", "lower.ci", "upper.ci",
                     "conf.level", "se_method") %in% names(res)))
  expect_true(res$lower.ci < res$upper.ci)
})

test_that("ses_calc: agresti vs fisher SE methods", {
  set.seed(108)
  x <- rnorm(30)
  y <- rnorm(30, mean = 0.3)

  res_ag <- ses_calc(x = x, y = y, ses = "rb", se_method = "agresti",
                     output = "data.frame")
  res_fi <- ses_calc(x = x, y = y, ses = "rb", se_method = "fisher",
                     output = "data.frame")

  # Point estimates should be the same
  expect_equal(res_ag$estimate, res_fi$estimate)
  # SEs may differ
  expect_true(res_ag$SE > 0)
  expect_true(res_fi$SE > 0)
})

test_that("ses_calc: paired samples", {
  data(sleep, envir = environment())
  res <- ses_calc(x = sleep$extra[sleep$group == 1],
                  y = sleep$extra[sleep$group == 2],
                  paired = TRUE, ses = "rb")
  expect_s3_class(res, "htest")
  expect_true(is.numeric(unname(res$estimate)))
})

test_that("ses_calc: one-sample", {
  set.seed(109)
  x <- rnorm(20, mean = 1)
  res <- ses_calc(x = x, ses = "rb")
  expect_s3_class(res, "htest")
  expect_true(is.numeric(unname(res$estimate)))
})

test_that("ses_calc: formula interface", {
  set.seed(110)
  df <- data.frame(
    val = c(rnorm(20), rnorm(20, mean = 0.5)),
    grp = factor(rep(c("A", "B"), each = 20))
  )
  res <- ses_calc(formula = val ~ grp, data = df, ses = "rb")
  expect_s3_class(res, "htest")
  expect_equal(res$data.name, "val by grp")
})

test_that("ses_calc: formula with paired warns", {
  set.seed(111)
  df <- data.frame(
    val = c(rnorm(20), rnorm(20, mean = 0.5)),
    grp = factor(rep(c("A", "B"), each = 20))
  )
  expect_message(
    ses_calc(formula = val ~ grp, data = df, paired = TRUE),
    "paired = TRUE.*formula interface.*not recommended"
  )
})

test_that("ses_calc: invalid alpha rejected", {
  x <- rnorm(10)
  expect_error(ses_calc(x = x, alpha = 0), "alpha must be a numeric value between 0 and 1")
  expect_error(ses_calc(x = x, alpha = 1), "alpha must be a numeric value between 0 and 1")
  expect_error(ses_calc(x = x, alpha = "a"), "alpha must be a numeric value between 0 and 1")
})

test_that("ses_calc: boundary correction messaging", {
  # Complete separation
  x <- c(10, 11, 12)
  y <- c(1, 2, 3)
  expect_message(
    ses_calc(x = x, y = y, ses = "rb", se_method = "agresti"),
    "Complete separation detected"
  )
})

test_that("ses_calc: hypothesis test on different SES scales", {
  set.seed(112)
  x <- rnorm(30)
  y <- rnorm(30, mean = 0.5)

  for (ses_type in c("rb", "cstat", "odds", "logodds")) {
    res <- ses_calc(x = x, y = y, ses = ses_type, alternative = "two.sided")
    expect_true(!is.null(res$p.value), info = paste("ses =", ses_type))
    expect_true(res$p.value >= 0 && res$p.value <= 1, info = paste("ses =", ses_type))
  }
})

test_that("ses_calc: fisher method hypothesis test", {
  set.seed(113)
  x <- rnorm(25, mean = 1)
  y <- rnorm(25)

  res <- ses_calc(x = x, y = y, ses = "rb",
                  se_method = "fisher", alternative = "two.sided")
  expect_true(!is.null(res$p.value))
  expect_true(!is.null(res$statistic))
})


# === boot_ses_calc tests ===

test_that("boot_ses_calc: default output is htest with no hypothesis test", {
  set.seed(200)
  x <- rnorm(15)
  y <- rnorm(15, mean = 0.5)

  res <- hush(boot_ses_calc(x = x, y = y, R = 99))
  expect_s3_class(res, "htest")
  expect_equal(res$alternative, "none")
  # No test components when alternative = "none"
  expect_null(res$statistic)
  expect_null(res$p.value)
  expect_null(res$null.value)
  # Should have estimate and CI
  expect_true(!is.null(res$estimate))
  expect_true(!is.null(res$conf.int))
  expect_true(!is.null(res$stderr))
  # Method should say "estimate with CI"
  expect_true(grepl("estimate with CI", res$method))
})

test_that("boot_ses_calc: hypothesis test when alternative specified", {
  set.seed(201)
  x <- rnorm(15)
  y <- rnorm(15, mean = 0.8)

  res <- hush(boot_ses_calc(x = x, y = y,
                            alternative = "two.sided",
                            R = 99))
  expect_s3_class(res, "htest")
  expect_equal(res$alternative, "two.sided")
  expect_true(!is.null(res$statistic))
  expect_true(!is.null(res$p.value))
  expect_true(!is.null(res$null.value))
  expect_true(res$p.value >= 0 && res$p.value <= 1)
  # Method should say "test"
  expect_true(grepl("test", res$method))
})

test_that("boot_ses_calc: point estimates match ses_calc", {
  set.seed(202)
  x <- rnorm(20)
  y <- rnorm(20, mean = 0.5)

  for (ses_type in c("rb", "cstat", "odds", "logodds")) {
    res_ses <- ses_calc(x = x, y = y, ses = ses_type, output = "data.frame")
    res_boot <- hush(boot_ses_calc(x = x, y = y, ses = ses_type, R = 99))
    expect_equal(unname(res_boot$estimate), res_ses$estimate,
                 info = paste("ses =", ses_type))
  }
})

test_that("boot_ses_calc: all bootstrap CI methods work", {
  set.seed(203)
  x <- rnorm(15)
  y <- rnorm(15, mean = 0.5)

  for (ci_type in c("basic", "stud", "perc")) {
    res <- hush(boot_ses_calc(x = x, y = y, boot_ci = ci_type, R = 99))
    expect_s3_class(res, "htest")
    expect_true(res$conf.int[1] < res$conf.int[2],
                info = paste("boot_ci =", ci_type))
  }
})

test_that("boot_ses_calc: data.frame output", {
  set.seed(204)
  x <- rnorm(15)
  y <- rnorm(15, mean = 0.5)

  res <- hush(boot_ses_calc(x = x, y = y, output = "data.frame", R = 99))
  expect_true(is.data.frame(res))
  expect_true(all(c("estimate", "SE", "lower.ci", "upper.ci",
                     "conf.level", "boot_ci") %in% names(res)))
})

test_that("boot_ses_calc: paired samples", {
  set.seed(205)
  x <- rnorm(15)
  y <- x + rnorm(15, mean = 0.3, sd = 1)

  res <- hush(boot_ses_calc(x = x, y = y, paired = TRUE, R = 99))
  expect_s3_class(res, "htest")
  expect_true(is.numeric(unname(res$estimate)))
})

test_that("boot_ses_calc: one-sample", {
  set.seed(206)
  x <- rnorm(15, mean = 1)

  res <- hush(boot_ses_calc(x = x, R = 99))
  expect_s3_class(res, "htest")
  expect_true(is.numeric(unname(res$estimate)))
})

test_that("boot_ses_calc: equivalence testing", {
  set.seed(207)
  x <- rnorm(20)
  y <- rnorm(20, mean = 0.1)

  res <- hush(boot_ses_calc(x = x, y = y, ses = "rb",
                            alternative = "equivalence",
                            null.value = c(-0.5, 0.5),
                            R = 99))
  expect_s3_class(res, "htest")
  expect_true(!is.null(res$p.value))
  expect_equal(length(res$null.value), 2)
  expect_equal(res$alternative, "equivalence")
})

test_that("boot_ses_calc: minimal effect testing", {
  set.seed(208)
  x <- rnorm(20, mean = 2)
  y <- rnorm(20)

  res <- hush(boot_ses_calc(x = x, y = y, ses = "rb",
                            alternative = "minimal.effect",
                            null.value = c(-0.3, 0.3),
                            R = 99))
  expect_s3_class(res, "htest")
  expect_true(!is.null(res$p.value))
  expect_equal(res$alternative, "minimal.effect")
})

test_that("boot_ses_calc: equivalence requires two bounds", {
  set.seed(209)
  x <- rnorm(15)
  y <- rnorm(15)

  expect_error(
    hush(boot_ses_calc(x = x, y = y, alternative = "equivalence",
                       null.value = 0, R = 99)),
    "null.value must be a vector of two values"
  )
})

test_that("boot_ses_calc: complete separation stops with error", {
  x <- c(10, 11, 12)
  y <- c(1, 2, 3)
  expect_error(
    hush(boot_ses_calc(x = x, y = y, R = 99)),
    "Complete separation detected"
  )
})

test_that("boot_ses_calc: formula interface", {
  set.seed(210)
  df <- data.frame(
    val = c(rnorm(15), rnorm(15, mean = 0.5)),
    grp = factor(rep(c("A", "B"), each = 15))
  )
  res <- hush(boot_ses_calc(formula = val ~ grp, data = df, R = 99))
  expect_s3_class(res, "htest")
  expect_equal(res$data.name, "val by grp")
})

test_that("boot_ses_calc: formula with paired warns", {
  set.seed(211)
  df <- data.frame(
    val = c(rnorm(15), rnorm(15, mean = 0.5)),
    grp = factor(rep(c("A", "B"), each = 15))
  )
  expect_message(
    hush(boot_ses_calc(formula = val ~ grp, data = df, paired = TRUE, R = 99)),
    "paired = TRUE.*formula interface.*not recommended"
  )
})

test_that("boot_ses_calc: one-sided alternatives", {
  set.seed(212)
  x <- rnorm(20, mean = 2)
  y <- rnorm(20)

  res_greater <- hush(boot_ses_calc(x = x, y = y, ses = "rb",
                                    alternative = "greater", R = 199))
  res_less <- hush(boot_ses_calc(x = x, y = y, ses = "rb",
                                 alternative = "less", R = 199))

  expect_true(res_greater$p.value < res_less$p.value)
})

test_that("boot_ses_calc: agresti vs fisher se_method", {
  set.seed(213)
  x <- rnorm(20)
  y <- rnorm(20, mean = 0.5)

  res_ag <- hush(boot_ses_calc(x = x, y = y, se_method = "agresti", R = 99))
  res_fi <- hush(boot_ses_calc(x = x, y = y, se_method = "fisher", R = 99))

  # Point estimates should match
  expect_equal(unname(res_ag$estimate), unname(res_fi$estimate))
})

test_that("boot_ses_calc: boot component included in htest output", {
  set.seed(214)
  x <- rnorm(15)
  y <- rnorm(15, mean = 0.5)

  res <- hush(boot_ses_calc(x = x, y = y, R = 99))
  expect_true(!is.null(res$boot))
  expect_equal(length(res$boot), 99)
})
