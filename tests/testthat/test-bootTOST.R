#context("Run Examples for boot_t_TOST")

# need hush function to run print through examples

hush = function(code) {
  sink("NUL") # use /dev/null in UNIX
  tmp = code
  sink()
  return(tmp)
}

test_that("Run examples for one sample", {

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

