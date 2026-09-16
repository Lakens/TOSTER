#context("Run Examples for boot_t_TOST")

# need hush function to run print through examples

hush = function(code) {
  sink(nullfile())
  tmp = code
  sink()
  return(tmp)
}

test_that("Run examples for one sample", {
  skip_on_cran()
  set.seed(31653464)

  samp1 = rnorm(33)

  expect_error(boot_t_TOST())
  expect_error(boot_t_TOST(x = samp1,
                           eqb = "test",
                           R = 99))
  expect_error(boot_t_TOST(x = samp1,
                           low_eqbound = -.5,
                           high_eqbound = .5,
                           alpha = 1.22,
                           R = 99))
  expect_error(boot_t_TOST(Sepal.Width ~ Species, data = iris))

  expect_message({
    boot_t_TOST(x = samp1,
                hypothesis = "MET",
                eqb = .005,
                R = 99)
  })



  htest_alt1 = boot_t_test(x = samp1,
                           alternative = "t",
                           R = 99)
  htest_alt2 = boot_t_test(x = samp1,
                           alternative = "g",
                           R = 99)
  htest_alt3 = boot_t_test(x = samp1,
                           alternative = "l",
                           R = 99)

  # Normal one sample ----

  test1 = boot_t_TOST(x = samp1,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 R = 99)
  set.seed(2649)
  test1 = boot_t_TOST(x = samp1,
                      eqb = .5,
                      R = 99)
  set.seed(2649)
  htest1 = boot_t_test(x = samp1,
                       mu = .5,
                       R = 99,
                       alternative = "e")
  expect_equal(test1$TOST$p.value[3],
               htest1$p.value)
  expect_error( boot_log_TOST(x = samp1,
                      eqb = .5,
                      R = 99))



  set.seed(432020)
  test3 = boot_t_TOST(x = samp1,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 hypothesis = "MET",
                 R = 99)
  set.seed(432020)
  htest3 = boot_t_test(x = samp1,
                       mu = .5,
                       R = 99,
                       alternative = "m")
  expect_equal(test3$TOST$p.value[3],
               htest3$p.value)

  expect_equal(1-test1$TOST$p.value[2],
               test3$TOST$p.value[2],
               tolerance = .1)

  expect_equal(1-test1$TOST$p.value[3],
               test3$TOST$p.value[3],
               tolerance = .1)

  expect_equal(test1$effsize$estimate,
               test3$effsize$estimate,
               ignore_attr = TRUE)

  expect_equal(1-test1$TOST$p.value[3],
               test3$TOST$p.value[3],
               tolerance = .1)

  # Re-run with bias correction not run -----
  test1 = boot_t_TOST(x = samp1,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 bias_correction = FALSE,
                 R = 99)

  test3 = boot_t_TOST(x = samp1,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 hypothesis = "MET",
                 bias_correction = FALSE,
                 R = 99)


  expect_equal(1-test1$TOST$p.value[3],
               test3$TOST$p.value[3],
               tolerance = .1)

  expect_equal(test1$effsize$estimate,
               test3$effsize$estimate,
               ignore_attr = TRUE)

  expect_equal(1-test1$TOST$p.value[3],
               test3$TOST$p.value[3],
               tolerance = .1)


  prtest = hush(print(test1))
  prtest2 = hush(describe(test1))
  p1 = plot(test1)
  p2 = expect_warning(plot(test1,
            type = "c"))

})


test_that("Run examples for two sample", {
  skip_on_cran()
  set.seed(76584441)

  samp1 = rnorm(25)
  samp2 = rnorm(25)

  htest_alt1 = boot_t_test(x = samp1,
                           y= samp2,
                           alternative = "t",
                           R = 99)
  htest_alt2 = boot_t_test(x = samp1,
                           y=samp2,
                           alternative = "g",
                           R = 99)
  htest_alt3 = boot_t_test(x = samp1,
                           y = samp2,
                           alternative = "l",
                           R = 99)
  htest_alt4 = boot_t_test(extra ~ group,
                           data = sleep,
                           alternative = "t",
                           R = 99)

  df_samp = data.frame(y = c(samp1,samp2),
                       group = c(rep("g1",25),
                                 rep("g2",25)))

  expect_error(boot_t_TOST())

  test1 = boot_t_TOST(x = samp1,
                 y = samp2,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 R = 199)

  test3 = boot_t_TOST(x = samp1,
                 y = samp2,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 hypothesis = "MET",
                 R = 199)
  prtest = hush(print(test3))
  prtest2 = hush(describe(test3))
  expect_equal(1-test1$TOST$p.value[2],
               test3$TOST$p.value[2],
               tolerance = .1)

  expect_equal(1-test1$TOST$p.value[3],
               test3$TOST$p.value[3],
               tolerance = .1)

  expect_equal(test1$effsize$estimate,
               test3$effsize$estimate,
               ignore_attr = TRUE)

  expect_equal(1-test1$TOST$p.value[3],
               test3$TOST$p.value[3],
               tolerance = .1)



  # Re-run with bias correction not run and non-Welch ----
  set.seed(232642)
  test1 = boot_t_TOST(x = samp1,
                 y = samp2,
                 var.equal = TRUE,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 bias_correction = FALSE,
                 R = 199)

  set.seed(232642)
  test1 = boot_t_TOST(x = samp1,
                      y = samp2,
                      var.equal = TRUE,
                      eqb = .5,
                      bias_correction = FALSE,
                      R = 199)

  set.seed(232642)
  test3 = boot_t_TOST(x = samp1,
                 y = samp2,
                 var.equal = TRUE,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 hypothesis = "MET",
                 bias_correction = FALSE,
                 R = 199)

  prtest = hush(print(test3))
  prtest2 = hush(describe(test3))
  expect_equal(1-test1$TOST$p.value[2],
               test3$TOST$p.value[2],
               tolerance = .1)

  expect_equal(1-test1$TOST$p.value[3],
               test3$TOST$p.value[3],
               tolerance = .1)

  expect_equal(test1$effsize$estimate,
               test3$effsize$estimate,
               ignore_attr = TRUE)

  expect_equal(1-test1$TOST$p.value[3],
               test3$TOST$p.value[3],
               tolerance = .1)

})


test_that("Run examples for paired samples", {
  skip_on_cran()
  set.seed(921387)

  samp1 = rnorm(25)
  samp2 = rnorm(25)

  htest_alt1 = boot_t_test(x = samp1,
                           y= samp2,
                           paired = TRUE,
                           alternative = "t",
                           R = 99)
  htest_alt2 = boot_t_test(x = samp1,
                           y=samp2,
                           paired = TRUE,
                           alternative = "g",
                           R = 99)
  htest_alt3 = boot_t_test(x = samp1,
                           y = samp2,
                           paired = TRUE,
                           alternative = "l",
                           R = 99)

  cor12 = stats::cor(samp1,samp2)

  df_samp = data.frame(y = c(samp1,samp2),
                       group = c(rep("g1",25),
                                 rep("g2",25)))

  expect_error(boot_t_TOST())

  test1 = boot_t_TOST(x = samp1,
                 y = samp2,
                 paired = TRUE,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 R = 199)

  test3 = boot_t_TOST(x = samp1,
                 y = samp2,
                 paired = TRUE,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 hypothesis = "MET",
                 R = 199)

  expect_equal(1-test1$TOST$p.value[2],
               test3$TOST$p.value[2],
               tolerance = .1)

  expect_equal(1-test1$TOST$p.value[3],
               test3$TOST$p.value[3],
               tolerance = .1)

  expect_equal(test1$effsize$estimate,
               test3$effsize$estimate,
               ignore_attr = TRUE)

  expect_equal(1-test1$TOST$p.value[3],
               test3$TOST$p.value[3],
               tolerance = .1)

})

# boot_t_TOST CI/p-value agreement tests -----

for (ci_method in c("perc", "basic", "bca", "stud")) {
  test_that(paste0("boot_t_TOST: CI/p-value agreement for ", ci_method), {
    skip_on_cran()

    set.seed(42)
    x <- rnorm(30, mean = 0.3)
    y <- rnorm(30)

    hush = function(code) {
      sink(nullfile())
      tmp = code
      sink()
      return(tmp)
    }

    res <- hush(boot_t_TOST(x = x, y = y, eqb = 2,
                             boot_ci = ci_method, R = 1999))
    # Check two-sided on 90% CI (TOST uses 1-2*alpha = 90%)
    # 90% CI excludes null iff two-sided p < 2*alpha = 0.10
    ci <- c(res$effsize$lower.ci[1], res$effsize$upper.ci[1])
    ci_excludes_null <- ci[1] > 0 || ci[2] < 0
    p_rejects <- res$TOST$p.value[1] < 0.10
    expect_equal(ci_excludes_null, p_rejects,
                 label = paste(ci_method, "CI/p agreement two.sided"))
  })
}


# boot_t_TOST vs boot_t_test -----
# With the same seed both functions draw identical resamples, so p-values
# and raw confidence intervals should match exactly for every design and
# CI method.

test_that("boot_t_TOST matches boot_t_test for all designs and CI methods", {
  skip_on_cran()

  set.seed(8421)
  x1 <- rnorm(20, mean = 7.4, sd = 1.2)
  x2 <- rnorm(22, mean = 5.6, sd = 1)
  y2 <- rnorm(18, mean = 5, sd = 1.3)
  xp <- rnorm(15, mean = 5.5, sd = 1)
  yp <- xp - rnorm(15, mean = 0.4, sd = 0.6)

  designs <- list(
    one = list(args = list(x = x1), bounds = c(7, 8.5), mu = 7.5),
    welch = list(args = list(x = x2, y = y2), bounds = c(0, 1.5), mu = 0.5),
    pooled = list(args = list(x = x2, y = y2, var.equal = TRUE),
                  bounds = c(0, 1.5), mu = 0.5),
    paired = list(args = list(x = xp, y = yp, paired = TRUE),
                  bounds = c(0, 1), mu = 0.3)
  )

  for (d in names(designs)) {
    des <- designs[[d]]
    for (ci in c("stud", "basic", "perc", "bca")) {
      lab <- paste(d, ci)

      set.seed(99)
      res <- suppressMessages(do.call(boot_t_TOST, c(des$args, list(
        eqb = des$bounds, mu = des$mu, R = 199, boot_ci = ci))))
      set.seed(99)
      eq <- do.call(boot_t_test, c(des$args, list(
        mu = des$bounds, alternative = "equivalence", R = 199, boot_ci = ci)))
      set.seed(99)
      nhst <- do.call(boot_t_test, c(des$args, list(
        mu = des$mu, R = 199, boot_ci = ci)))
      set.seed(99)
      lower <- do.call(boot_t_test, c(des$args, list(
        mu = des$bounds[1], alternative = "greater", R = 199, boot_ci = ci)))
      set.seed(99)
      upper <- do.call(boot_t_test, c(des$args, list(
        mu = des$bounds[2], alternative = "less", R = 199, boot_ci = ci)))

      expect_equal(res$TOST$p.value[1], nhst$p.value, label = lab)
      expect_equal(res$TOST$p.value[2], lower$p.value, label = lab)
      expect_equal(res$TOST$p.value[3], upper$p.value, label = lab)
      expect_equal(max(res$TOST$p.value[2:3]), eq$p.value, label = lab)
      expect_equal(res$effsize$estimate[1], unname(nhst$estimate[length(nhst$estimate)]),
                   label = lab)
      expect_equal(c(res$effsize$lower.ci[1], res$effsize$upper.ci[1]),
                   as.numeric(eq$conf.int), label = lab)
    }
  }
})
