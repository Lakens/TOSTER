#context("Run Examples for t_TOST")

# need hush function to run print through examples

hush = function(code) {
  sink("NUL") # use /dev/null in UNIX
  tmp = code
  sink()
  return(tmp)
}

test_that("Run examples for one sample", {

  hush = function(code) {
    sink("NUL") # use /dev/null in UNIX
    tmp = code
    sink()
    return(tmp)
  }


  set.seed(3164964)

  samp1 = rnorm(33)

  expect_error(t_TOST())
  expect_error(t_TOST(x = samp1,
                             eqb = c(-1,1,.5)))
  expect_error(t_TOST(x = samp1,
                     low_eqbound = -.5,
                     high_eqbound = .5,
                     alpha = 1.22))
  expect_error(t_TOST(Sepal.Width ~ Species, data = iris))

  expect_message({t_TOST(x = samp1,
                 eqb = .05,
                 hypothesis = "MET")})
  # Normal one sample ----

  test1 = t_TOST(x = samp1,
                 low_eqbound = -.5,
                 high_eqbound = .5)
  test1_smd = smd_calc(x = samp1,
                       alpha = .1)
  test1_smd_boot = boot_smd_calc(x = samp1,
                       alpha = .1,
                       R = 99)
  expect_equal(test1_smd$estimate, test1_smd_boot$estimate)
  expect_equal(test1_smd$estimate, test1$effsize$estimate[2])
  expect_equal(test1_smd$lower.ci, test1$effsize$lower.ci[2])
  expect_equal(test1_smd$upper.ci, test1$effsize$upper.ci[2])
  expect_equal(test1_smd$SE, test1$effsize$SE[2])

  test1 = t_TOST(x = samp1,
                 eqb = .5)

  test1 = t_TOST(x = samp1,
                 eqb = c(-.5,.5))

  test2 = suppressMessages(t_TOST(x = samp1,
                                #low_eqbound = -.5,
                                eqb = .5,
                                eqbound_type = "SMD"))

  test3 = t_TOST(x = samp1,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 hypothesis = "MET")

  test4 = suppressMessages( { t_TOST(x = samp1,
                                   #low_eqbound = -.5,
                                   eqb = .5,
                                   eqbound_type = "SMD",
                                   hypothesis = "MET")})

  expect_equal(1-test1$TOST$p.value[2],
               test3$TOST$p.value[2])

  expect_equal(1-test1$TOST$p.value[3],
               test3$TOST$p.value[3])

  expect_equal(1-test2$TOST$p.value[2],
               test4$TOST$p.value[2])

  expect_equal(1-test2$TOST$p.value[3],
               test4$TOST$p.value[3])

  # Compare to tsum --------

  tsum1 = tsum_TOST(m1 = mean(samp1),
                    sd1 = sd(samp1),
                    n1 = length(samp1),
                    low_eqbound = -.5,
                    high_eqbound = .5)
  tsum1 = tsum_TOST(m1 = mean(samp1),
                    sd1 = sd(samp1),
                    n1 = length(samp1),
                    eqb = .5)
  tsum1_tg = tsum_TOST(m1 = mean(samp1),
                    sd1 = sd(samp1),
                    n1 = length(samp1),
                    eqb = .5,
                    glass = "glass1")
  tsum1_tg = tsum_TOST(m1 = mean(samp1),
                       sd1 = sd(samp1),
                       n1 = length(samp1),
                       eqb = .5,
                       glass = "glass2")
  expect_error(tsum_TOST(m1 = mean(samp1),
                         sd1 = sd(samp1),
                         n1 = length(samp1),
                         eqb = c(.5,-.5,1)))
  expect_error(tsum_TOST(m1 = mean(samp1),
                         sd1 = sd(samp1),
                         n1 = length(samp1),
                         m12 = mean(samp1),
                         sd12 = sd(samp1),
                         n12 = length(samp1),
                         eqb = c(.5,-.5),
                         paired = TRUE))

  expect_warning(tsum_TOST(m1 = mean(samp1),
                         sd1 = sd(samp1),
                         n1 = length(samp1),
                         m2 = mean(samp1),
                         sd2 = sd(samp1),
                         r12 = .73,
                         n2 = 999,
                         eqb = c(.5,-.5),
                         paired = TRUE))
  expect_error(tsum_TOST(m1 = mean(samp1),
                    sd1 = sd(samp1),
                    n1 = length(samp1)))
  tsum2 = suppressMessages({ tsum_TOST(m1 = mean(samp1),
                                     sd1 = sd(samp1),
                                     n1 = length(samp1),
                                     low_eqbound = -.5,
                                     high_eqbound = .5,
                                     eqbound_type = "SMD") })

  tsum3 = tsum_TOST(m1 = mean(samp1),
                    sd1 = sd(samp1),
                    n1 = length(samp1),
                    low_eqbound = -.5,
                    high_eqbound = .5,
                    hypothesis = "MET")

  tsum4 = suppressMessages({tsum_TOST(m1 = mean(samp1),
                                    sd1 = sd(samp1),
                                    n1 = length(samp1),
                                    low_eqbound = -.5,
                                    high_eqbound = .5,
                                    eqbound_type = "SMD",
                                    hypothesis = "MET") })
  # Check internal consistency
  expect_equal(test1$TOST$p.value,
               tsum1$TOST$p.value,
               ignore_attr = TRUE)

  expect_equal(test2$TOST$p.value,
               tsum2$TOST$p.value,
               ignore_attr = TRUE)

  expect_equal(test3$TOST$p.value,
               tsum3$TOST$p.value,
               ignore_attr = TRUE)

  expect_equal(test4$TOST$p.value,
               tsum4$TOST$p.value,
               ignore_attr = TRUE)


  # Re-run with bias correction not run -----
  test1 = t_TOST(x = samp1,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 bias_correction = FALSE)

  test2 = suppressMessages({  t_TOST(x = samp1,
                                   low_eqbound = -.5,
                                   high_eqbound = .5,
                                   eqbound_type = "SMD",
                                   bias_correction = FALSE)
  })

  test3 = t_TOST(x = samp1,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 hypothesis = "MET",
                 bias_correction = FALSE)

  test4 = suppressMessages( { t_TOST(x = samp1,
                                   low_eqbound = -.5,
                                   high_eqbound = .5,
                                   eqbound_type = "SMD",
                                   hypothesis = "MET",
                                   bias_correction = FALSE,
                                   smd_ci = "goulet")
  })

  expect_equal(1-test1$TOST$p.value[2],
               test3$TOST$p.value[2])

  expect_equal(1-test1$TOST$p.value[3],
               test3$TOST$p.value[3])

  expect_equal(1-test2$TOST$p.value[2],
               test4$TOST$p.value[2])

  expect_equal(1-test2$TOST$p.value[3],
               test4$TOST$p.value[3])

  tsum1 = tsum_TOST(m1 = mean(samp1),
                    sd1 = sd(samp1),
                    n1 = length(samp1),
                    low_eqbound = -.5,
                    high_eqbound = .5,
                    bias_correction = FALSE)

  tsum2 = suppressMessages( tsum_TOST(m1 = mean(samp1),
                                    sd1 = sd(samp1),
                                    n1 = length(samp1),
                                    low_eqbound = -.5,
                                    high_eqbound = .5,
                                    eqbound_type = "SMD",
                                    bias_correction = FALSE) )

  tsum3 = tsum_TOST(m1 = mean(samp1),
                    sd1 = sd(samp1),
                    n1 = length(samp1),
                    low_eqbound = -.5,
                    high_eqbound = .5,
                    hypothesis = "MET",
                    bias_correction = FALSE)

  tsum4 = suppressMessages(tsum_TOST(m1 = mean(samp1),
                                   sd1 = sd(samp1),
                                   n1 = length(samp1),
                                   low_eqbound = -.5,
                                   high_eqbound = .5,
                                   eqbound_type = "SMD",
                                   hypothesis = "MET",
                                   bias_correction = FALSE))
  # Check internal consistency
  expect_equal(test1$TOST$p.value,
               tsum1$TOST$p.value,
               ignore_attr = TRUE)

  expect_equal(test2$TOST$p.value,
               tsum2$TOST$p.value,
               ignore_attr = TRUE)

  expect_equal(test3$TOST$p.value,
               tsum3$TOST$p.value,
               ignore_attr = TRUE)

  expect_equal(test4$TOST$p.value,
               tsum4$TOST$p.value,
               ignore_attr = TRUE)

  prtest = hush(print(test4))
  des = hush(describe(test4))
  p1 = plot(test4)
  p2 = plot(test4,
            type = "c")
  p3 = plot(test4,
            type = "tnull")
  p4 = plot(test4,
            type = "cd")

})


test_that("Run examples for two sample", {

  hush = function(code) {
    sink("NUL") # use /dev/null in UNIX
    tmp = code
    sink()
    return(tmp)
  }

  tt1 = t_TOST(extra ~ group,
               data = sleep,
               eqb = 2)
  des = describe(tt1)


  set.seed(651466441)

  samp1 = rnorm(25)
  samp2 = rnorm(25)

  df_samp = data.frame(y = c(samp1,samp2),
                       group = c(rep("g1",25),
                                 rep("g2",25)))

  expect_error(t_TOST())

  test1 = t_TOST(x = samp1,
                 y = samp2,
                 low_eqbound = -.5,
                 high_eqbound = .5)

  test1_smd = smd_calc(x = samp1,
                              y = samp2,
                       alpha = .1)

  test1_smd_boot_stud = boot_smd_calc(x = samp1,
                       y = samp2,
                       alpha = .1,
                       boot_ci = "stud",
                       R = 99)
  test1_smd_boot_basic = boot_smd_calc(x = samp1,
                                      y = samp2,
                                      alpha = .1,
                                      boot_ci = "b",
                                      R = 99)
  expect_error(boot_smd_calc(x = samp1,
                             y = samp2,
                             alpha = .1,
                             boot_ci = "n",
                             R = 99))
  test1_smd_boot_perc = boot_smd_calc(x = samp1,
                                      y = samp2,
                                      alpha = .1,
                                      boot_ci = "p",
                                      R = 99)
  expect_equal(rep(test1_smd$estimate,3),
               c(test1_smd_boot_stud$estimate,
                 test1_smd_boot_basic$estimate,
                 test1_smd_boot_perc$estimate))

  expect_error(smd_calc(x = samp1,
                        y = samp2,
                        alpha = -.1))

  expect_equal(test1_smd$estimate, test1$effsize$estimate[2])
  expect_equal(test1_smd$lower.ci, test1$effsize$lower.ci[2])
  expect_equal(test1_smd$upper.ci, test1$effsize$upper.ci[2])
  expect_equal(test1_smd$SE, test1$effsize$SE[2])

  test2 = suppressMessages( t_TOST(x = samp1,
                                 y = samp2,
                                 low_eqbound = -.5,
                                 high_eqbound = .5,
                                 eqbound_type = "SMD") )

  test3 = t_TOST(x = samp1,
                 y = samp2,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 hypothesis = "MET")

  test4 = suppressMessages( t_TOST(x = samp1,
                                 y = samp2,
                                 low_eqbound = -.5,
                                 high_eqbound = .5,
                                 eqbound_type = "SMD",
                                 hypothesis = "MET") )

  expect_equal(1-test1$TOST$p.value[2],
               test3$TOST$p.value[2])

  expect_equal(1-test1$TOST$p.value[3],
               test3$TOST$p.value[3])

  expect_equal(1-test2$TOST$p.value[2],
               test4$TOST$p.value[2])

  expect_equal(1-test2$TOST$p.value[3],
               test4$TOST$p.value[3])

  # t_sum comparison var.equal = FALSE ----

  tsum1 = tsum_TOST(m1 = mean(samp1),
                    sd1 = sd(samp1),
                    n1 = length(samp1),
                    m2 = mean(samp2),
                    sd2 = sd(samp2),
                    n2 = length(samp2),
                    low_eqbound = -.5,
                    high_eqbound = .5,
                    bias_correction = FALSE)

  tsum2 = suppressMessages( {  tsum_TOST(m1 = mean(samp1),
                                       sd1 = sd(samp1),
                                       n1 = length(samp1),
                                       m2 = mean(samp2),
                                       sd2 = sd(samp2),
                                       n2 = length(samp2),
                                       low_eqbound = -.5,
                                       high_eqbound = .5,
                                       eqbound_type = "SMD",
                                       bias_correction = FALSE) })

  tsum3 = tsum_TOST(m1 = mean(samp1),
                    sd1 = sd(samp1),
                    n1 = length(samp1),
                    m2 = mean(samp2),
                    sd2 = sd(samp2),
                    n2 = length(samp2),
                    low_eqbound = -.5,
                    high_eqbound = .5,
                    hypothesis = "MET",
                    bias_correction = FALSE)

  tsum4 = suppressMessages(  { tsum_TOST(m1 = mean(samp1),
                                       sd1 = sd(samp1),
                                       n1 = length(samp1),
                                       m2 = mean(samp2),
                                       sd2 = sd(samp2),
                                       n2 = length(samp2),
                                       low_eqbound = -.5,
                                       high_eqbound = .5,
                                       eqbound_type = "SMD",
                                       hypothesis = "MET",
                                       bias_correction = FALSE)})
  # Check internal consistency
  expect_equal(test1$TOST$p.value,
               tsum1$TOST$p.value,
               ignore_attr = TRUE)

  expect_equal(test2$TOST$p.value,
               tsum2$TOST$p.value,
               ignore_attr = TRUE)

  expect_equal(test3$TOST$p.value,
               tsum3$TOST$p.value,
               ignore_attr = TRUE)

  expect_equal(test4$TOST$p.value,
               tsum4$TOST$p.value,
               ignore_attr = TRUE)

  # Re-run with bias correction not run and non-Welch ----
  test1 = t_TOST(x = samp1,
                 y = samp2,
                 var.equal = TRUE,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 bias_correction = FALSE)

  test2 = suppressMessages( t_TOST(x = samp1,
                                 y = samp2,
                                 var.equal = TRUE,
                                 low_eqbound = -.5,
                                 high_eqbound = .5,
                                 eqbound_type = "SMD",
                                 bias_correction = FALSE) )

  test3 = t_TOST(x = samp1,
                 y = samp2,
                 var.equal = TRUE,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 hypothesis = "MET",
                 bias_correction = FALSE)

  test4 = suppressMessages(  t_TOST(x = samp1,
                                  y = samp2,
                                  var.equal = TRUE,
                                  low_eqbound = -.5,
                                  high_eqbound = .5,
                                  eqbound_type = "SMD",
                                  hypothesis = "MET",
                                  bias_correction = FALSE) )

  expect_equal(1-test1$TOST$p.value[2],
               test3$TOST$p.value[2])

  expect_equal(1-test1$TOST$p.value[3],
               test3$TOST$p.value[3])

  expect_equal(1-test2$TOST$p.value[2],
               test4$TOST$p.value[2])

  expect_equal(1-test2$TOST$p.value[3],
               test4$TOST$p.value[3])

  # t_sum comparison var.equal = TRUE ----

  tsum1 = tsum_TOST(m1 = mean(samp1),
                    sd1 = sd(samp1),
                    n1 = length(samp1),
                    m2 = mean(samp2),
                    sd2 = sd(samp2),
                    n2 = length(samp2),
                    var.equal = TRUE,
                    low_eqbound = -.5,
                    high_eqbound = .5,
                    bias_correction = FALSE)

  tsum2 = suppressMessages({  tsum_TOST(m1 = mean(samp1),
                                      sd1 = sd(samp1),
                                      n1 = length(samp1),
                                      m2 = mean(samp2),
                                      sd2 = sd(samp2),
                                      n2 = length(samp2),
                                      var.equal = TRUE,
                                      low_eqbound = -.5,
                                      high_eqbound = .5,
                                      eqbound_type = "SMD",
                                      bias_correction = FALSE) })

  tsum3 = tsum_TOST(m1 = mean(samp1),
                    sd1 = sd(samp1),
                    n1 = length(samp1),
                    m2 = mean(samp2),
                    sd2 = sd(samp2),
                    n2 = length(samp2),
                    var.equal = TRUE,
                    low_eqbound = -.5,
                    high_eqbound = .5,
                    hypothesis = "MET",
                    bias_correction = FALSE)

  tsum4 = suppressMessages(  {
    tsum_TOST(m1 = mean(samp1),
              sd1 = sd(samp1),
              n1 = length(samp1),
              m2 = mean(samp2),
              sd2 = sd(samp2),
              n2 = length(samp2),
              var.equal = TRUE,
              low_eqbound = -.5,
              high_eqbound = .5,
              eqbound_type = "SMD",
              hypothesis = "MET",
              bias_correction = FALSE)})
  # Check internal consistency
  expect_equal(test1$TOST$p.value,
               tsum1$TOST$p.value,
               ignore_attr = TRUE)

  expect_equal(test2$TOST$p.value,
               tsum2$TOST$p.value,
               ignore_attr = TRUE)

  expect_equal(test3$TOST$p.value,
               tsum3$TOST$p.value,
               ignore_attr = TRUE)

  expect_equal(test4$TOST$p.value,
               tsum4$TOST$p.value,
               ignore_attr = TRUE)

  # Run with formula
  test1 = t_TOST(formula = y ~ group,
                 data = df_samp,
                 var.equal = TRUE,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 bias_correction = FALSE)
  test1_smd = smd_calc(formula = y ~ group,
                       data = df_samp,
                       var.equal = TRUE,
                       bias_correction = FALSE)
  test1_smd_boot_stud = boot_smd_calc(x = samp1,
                                      y = samp2,
                                      alpha = .1,
                                      boot_ci = "stud",
                                      R = 99,
                                      var.equal = TRUE,
                                      bias_correction = FALSE)
  test1_smd_boot_basic = boot_smd_calc(x = samp1,
                                       y = samp2,
                                       alpha = .1,
                                       boot_ci = "b",
                                       R = 99,
                                       var.equal = TRUE,
                                       bias_correction = FALSE)
  expect_error(boot_smd_calc(x = samp1,
                             y = samp2,
                             alpha = .1,
                             boot_ci = "n",
                             R = 99,
                             var.equal = TRUE,
                             bias_correction = FALSE))
  test1_smd_boot_perc = boot_smd_calc(x = samp1,
                                      y = samp2,
                                      alpha = .1,
                                      boot_ci = "p",
                                      R = 99,
                                      var.equal = TRUE,
                                      bias_correction = FALSE)
  expect_equal(rep(test1_smd$estimate,3),
               c(test1_smd_boot_stud$estimate,
                 test1_smd_boot_basic$estimate,
                 test1_smd_boot_perc$estimate))
  # test htest
  ash = as_htest(test1)
  test2 = suppressMessages( t_TOST(formula = y ~ group,
                                 data = df_samp,
                                 var.equal = TRUE,
                                 low_eqbound = -.5,
                                 high_eqbound = .5,
                                 eqbound_type = "SMD",
                                 bias_correction = FALSE) )
  ash = as_htest(test2)
  test3 = t_TOST(formula = y ~ group,
                 data = df_samp,
                 var.equal = TRUE,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 hypothesis = "MET",
                 bias_correction = FALSE)

  test4 = suppressMessages( t_TOST(formula = y ~ group,
                                 data = df_samp,
                                 var.equal = TRUE,
                                 low_eqbound = -.5,
                                 high_eqbound = .5,
                                 eqbound_type = "SMD",
                                 hypothesis = "MET",
                                 bias_correction = FALSE,
                                 smd_ci = "g") )

  expect_equal(1-test1$TOST$p.value[2],
               test3$TOST$p.value[2])

  expect_equal(1-test1$TOST$p.value[3],
               test3$TOST$p.value[3])

  expect_equal(1-test2$TOST$p.value[2],
               test4$TOST$p.value[2])

  expect_equal(1-test2$TOST$p.value[3],
               test4$TOST$p.value[3])

  # Check internal consistency
  expect_equal(test1$TOST$p.value,
               tsum1$TOST$p.value,
               ignore_attr = TRUE)

  expect_equal(test2$TOST$p.value,
               tsum2$TOST$p.value,
               ignore_attr = TRUE)

  expect_equal(test3$TOST$p.value,
               tsum3$TOST$p.value,
               ignore_attr = TRUE)

  expect_equal(test4$TOST$p.value,
               tsum4$TOST$p.value,
               ignore_attr = TRUE)

  prtest = hush(print(test4))
  prtest2 = hush(describe(test4))
  des = hush(describe(test4))
  p1 = plot(test4)

})


test_that("Run examples for paired samples", {

  set.seed(789461245)

  samp1 = rnorm(25)
  samp2 = rnorm(25)

  cor12 = stats::cor(samp1,samp2)

  df_samp = data.frame(y = c(samp1,samp2),
                       group = c(rep("g1",25),
                                 rep("g2",25)))

  expect_error(t_TOST())

  test1 = t_TOST(x = samp1,
                 y = samp2,
                 paired = TRUE,
                 low_eqbound = -.5,
                 high_eqbound = .5)
  test1_smd = smd_calc(x = samp1,
                 y = samp2,
                 paired = TRUE,
                 alpha = .1)
  test1_smd_boot_stud = boot_smd_calc(x = samp1,
                                      y = samp2,
                                      paired = TRUE,
                                      alpha = .1,
                                      R = 99)
  test1_smd_boot_basic = boot_smd_calc(x = samp1,
                                       y = samp2,
                                       paired = TRUE,
                                       alpha = .1,
                                       boot_ci = "b",
                                       R = 99)
  expect_error(boot_smd_calc(x = samp1,
                             y = samp2,
                             alpha = .1,
                             boot_ci = "n",
                             paired = TRUE,
                             R = 99,
                             var.equal = TRUE,
                             bias_correction = FALSE))
  test1_smd_boot_perc = boot_smd_calc(x = samp1,
                                      y = samp2,
                                      paired = TRUE,
                                      alpha = .1,
                                      boot_ci = "p",
                                      R = 99)
  expect_equal(rep(test1_smd$estimate,3),
               c(test1_smd_boot_stud$estimate,
                 test1_smd_boot_basic$estimate,
                 test1_smd_boot_perc$estimate))

  expect_equal(test1_smd$estimate, test1$effsize$estimate[2])
  expect_equal(test1_smd$lower.ci, test1$effsize$lower.ci[2])
  expect_equal(test1_smd$upper.ci, test1$effsize$upper.ci[2])
  expect_equal(test1_smd$SE, test1$effsize$SE[2])
  ash = as_htest(test1)
  test2 = suppressMessages(  t_TOST(x = samp1,
                                  y = samp2,
                                  paired = TRUE,
                                  low_eqbound = -.5,
                                  high_eqbound = .5,
                                  eqbound_type = "SMD") )
  ash = as_htest(test2)
  test3 = t_TOST(x = samp1,
                 y = samp2,
                 paired = TRUE,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 hypothesis = "MET")

  test4 = suppressMessages( t_TOST(x = samp1,
                                 y = samp2,
                                 paired = TRUE,
                                 low_eqbound = -.5,
                                 high_eqbound = .5,
                                 eqbound_type = "SMD",
                                 hypothesis = "MET") )

  expect_equal(1-test1$TOST$p.value[2],
               test3$TOST$p.value[2])

  expect_equal(1-test1$TOST$p.value[3],
               test3$TOST$p.value[3])

  expect_equal(1-test2$TOST$p.value[2],
               test4$TOST$p.value[2])

  expect_equal(1-test2$TOST$p.value[3],
               test4$TOST$p.value[3])

  # t_sum paired ----

  tsum1 = tsum_TOST(m1 = mean(samp1),
                    sd1 = sd(samp1),
                    n1 = length(samp1),
                    m2 = mean(samp2),
                    sd2 = sd(samp2),
                    n2 = length(samp2),
                    r12 = cor12, paired = TRUE,
                    low_eqbound = -.5,
                    high_eqbound = .5,
                    bias_correction = FALSE)


  tsum2 = suppressMessages(  {  tsum_TOST(m1 = mean(samp1),
                                        sd1 = sd(samp1),
                                        n1 = length(samp1),
                                        m2 = mean(samp2),
                                        sd2 = sd(samp2),
                                        n2 = length(samp2),
                                        r12 = cor12, paired = TRUE,
                                        low_eqbound = -.5,
                                        high_eqbound = .5,
                                        eqbound_type = "SMD",
                                        bias_correction = FALSE) })

  tsum3 = tsum_TOST(m1 = mean(samp1),
                    sd1 = sd(samp1),
                    n1 = length(samp1),
                    m2 = mean(samp2),
                    sd2 = sd(samp2),
                    n2 = length(samp2),
                    r12 = cor12,
                    paired = TRUE,
                    low_eqbound = -.5,
                    high_eqbound = .5,
                    hypothesis = "MET",
                    bias_correction = FALSE)

  expect_equal(cor12,
               extract_r_paired(m1 = mean(samp1),
                                sd1 = sd(samp1),
                                m2 = mean(samp2),
                                sd2 = sd(samp2),
                                n = length(samp2),
                                tstat = tsum3$TOST$t[1]))
  expect_error(extract_r_paired(m1 = mean(samp1),
                                sd1 = sd(samp1),
                                m2 = mean(samp2),
                                sd2 = sd(samp2),
                                n = length(samp2),
                                tstat = NULL))
  test_ragain = extract_r_paired(m1 = mean(samp1),
                                 sd1 = sd(samp1),
                                 m2 = mean(samp2),
                                 #sd2 = sd(samp2),
                                 n = length(samp2),
                                 tstat = tsum3$TOST$t[1])

  something = extract_r_paired(m1 = mean(samp1),
                               sd1 = sd(samp1),
                               m2 = mean(samp2),
                               #sd2 = sd(samp2),
                               n = length(samp2),
                               tstat = tsum3$TOST$t[1])
  tsum4 = suppressMessages( { tsum_TOST(m1 = mean(samp1),
                                      sd1 = sd(samp1),
                                      n1 = length(samp1),
                                      m2 = mean(samp2),
                                      sd2 = sd(samp2),
                                      n2 = length(samp2),
                                      r12 = cor12, paired = TRUE,
                                      low_eqbound = -.5,
                                      high_eqbound = .5,
                                      eqbound_type = "SMD",
                                      hypothesis = "MET",
                                      bias_correction = FALSE)})

  expect_equal(test1$TOST$p.value,
               tsum1$TOST$p.value,
               ignore_attr = TRUE)

  expect_equal(test2$TOST$p.value,
               tsum2$TOST$p.value,
               ignore_attr = TRUE)

  expect_equal(test3$TOST$p.value,
               tsum3$TOST$p.value,
               ignore_attr = TRUE)

  expect_equal(test4$TOST$p.value,
               tsum4$TOST$p.value,
               ignore_attr = TRUE)

  # Re-run with bias correction not run and rm_correction
  # rm_correction = TRUE
  test1 = t_TOST(x = samp1,
                 y = samp2,
                 paired = TRUE,
                 rm_correction = TRUE,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 bias_correction = FALSE)

  test2 = suppressMessages( t_TOST(x = samp1,
                                 y = samp2,
                                 paired = TRUE,
                                 rm_correction = TRUE,
                                 low_eqbound = -.5,
                                 high_eqbound = .5,
                                 eqbound_type = "SMD",
                                 bias_correction = FALSE) )

  test3 = t_TOST(x = samp1,
                 y = samp2,
                 paired = TRUE,
                 rm_correction = TRUE,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 hypothesis = "MET",
                 bias_correction = FALSE)

  test4 = suppressMessages( t_TOST(x = samp1,
                                 y = samp2,
                                 paired = TRUE,
                                 rm_correction = TRUE,
                                 low_eqbound = -.5,
                                 high_eqbound = .5,
                                 eqbound_type = "SMD",
                                 hypothesis = "MET",
                                 bias_correction = FALSE) )

  expect_equal(1-test1$TOST$p.value[2],
               test3$TOST$p.value[2])

  expect_equal(1-test1$TOST$p.value[3],
               test3$TOST$p.value[3])

  expect_equal(1-test2$TOST$p.value[2],
               test4$TOST$p.value[2])

  expect_equal(1-test2$TOST$p.value[3],
               test4$TOST$p.value[3])

  # t_sum paired ----

  tsum1 = tsum_TOST(m1 = mean(samp1),
                    sd1 = sd(samp1),
                    n1 = length(samp1),
                    m2 = mean(samp2),
                    sd2 = sd(samp2),
                    n2 = length(samp2),
                    r12 = cor12, paired = TRUE,
                    low_eqbound = -.5,
                    high_eqbound = .5,
                    rm_correction = TRUE,
                    bias_correction = FALSE)


  tsum2 = suppressMessages( {  tsum_TOST(m1 = mean(samp1),
                                       sd1 = sd(samp1),
                                       n1 = length(samp1),
                                       m2 = mean(samp2),
                                       sd2 = sd(samp2),
                                       n2 = length(samp2),
                                       r12 = cor12, paired = TRUE,
                                       low_eqbound = -.5,
                                       high_eqbound = .5,
                                       eqbound_type = "SMD",
                                       rm_correction = TRUE,
                                       bias_correction = FALSE) })

  tsum3 = tsum_TOST(m1 = mean(samp1),
                    sd1 = sd(samp1),
                    n1 = length(samp1),
                    m2 = mean(samp2),
                    sd2 = sd(samp2),
                    n2 = length(samp2),
                    r12 = cor12,
                    paired = TRUE,
                    low_eqbound = -.5,
                    high_eqbound = .5,
                    hypothesis = "MET",
                    rm_correction = TRUE,
                    bias_correction = FALSE)

  tsum4 = suppressMessages( { tsum_TOST(m1 = mean(samp1),
                                      sd1 = sd(samp1),
                                      n1 = length(samp1),
                                      m2 = mean(samp2),
                                      sd2 = sd(samp2),
                                      n2 = length(samp2),
                                      r12 = cor12, paired = TRUE,
                                      low_eqbound = -.5,
                                      high_eqbound = .5,
                                      eqbound_type = "SMD",
                                      hypothesis = "MET",
                                      rm_correction = TRUE,
                                      bias_correction = FALSE)})

  expect_equal(test1$TOST$p.value,
               tsum1$TOST$p.value,
               ignore_attr = TRUE)

  expect_equal(test2$TOST$p.value,
               tsum2$TOST$p.value,
               ignore_attr = TRUE)

  expect_equal(test3$TOST$p.value,
               tsum3$TOST$p.value,
               ignore_attr = TRUE)

  expect_equal(test4$TOST$p.value,
               tsum4$TOST$p.value,
               ignore_attr = TRUE)

  # Use vector form for paired tests, not formula
  test1 = t_TOST(x = df_samp$y[df_samp$group == "g1"],
                 y = df_samp$y[df_samp$group == "g2"],
                 paired = TRUE,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 bias_correction = FALSE)

  test2 = suppressMessages( t_TOST(x = df_samp$y[df_samp$group == "g1"],
                                 y = df_samp$y[df_samp$group == "g2"],
                                 paired = TRUE,
                                 low_eqbound = -.5,
                                 high_eqbound = .5,
                                 eqbound_type = "SMD",
                                 bias_correction = FALSE) )

  test3 = t_TOST(x = df_samp$y[df_samp$group == "g1"],
                 y = df_samp$y[df_samp$group == "g2"],
                 paired = TRUE,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 hypothesis = "MET",
                 bias_correction = FALSE)

  test4 = suppressMessages(  t_TOST(x = df_samp$y[df_samp$group == "g1"],
                                  y = df_samp$y[df_samp$group == "g2"],
                                  paired = TRUE,
                                  low_eqbound = -.5,
                                  high_eqbound = .5,
                                  eqbound_type = "SMD",
                                  hypothesis = "MET",
                                  bias_correction = FALSE) )

  expect_equal(1-test1$TOST$p.value[2],
               test3$TOST$p.value[2])

  expect_equal(1-test1$TOST$p.value[3],
               test3$TOST$p.value[3])

  expect_equal(1-test2$TOST$p.value[2],
               test4$TOST$p.value[2])

  expect_equal(1-test2$TOST$p.value[3],
               test4$TOST$p.value[3])

  prtest = hush(print(test4))
  des = hush(describe(test4))
  p1 = plot(test4)

})

test_that("Run examples for plot_smd", {

  set.seed(1776)

  samp1 = rnorm(25)
  samp2 = rnorm(25)

  df_samp = data.frame(y = c(samp1,samp2),
                       group = c(rep("g1",25),
                                 rep("g2",25)))

  expect_error(t_TOST())

  test1 = t_TOST(x = samp1,
                 y = samp2,
                 paired = TRUE,
                 low_eqbound = -.5,
                 high_eqbound = .5)

  test2 = suppressMessages(  t_TOST(x = samp1,
                                  y = samp2,
                                  paired = TRUE,
                                  low_eqbound = -.5,
                                  high_eqbound = .5,
                                  eqbound_type = "SMD") )

  test3 = t_TOST(x = samp1,
                 y = samp2,
                 paired = TRUE,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 hypothesis = "MET")

  test4 = suppressMessages( t_TOST(x = samp1,
                                 y = samp2,
                                 paired = TRUE,
                                 low_eqbound = -.5,
                                 high_eqbound = .5,
                                 eqbound_type = "SMD",
                                 hypothesis = "MET") )



  p1 = plot_smd(lambda = c(test1$smd$d_lambda),
                df = c(test1$smd$d_df),
                d = c(test1$smd$d),
                type = "cd",
                smd_ci = "goulet")

  p1 = plot_smd(lambda = c(test1$smd$d_lambda),
                df = c(test1$smd$d_df),
                d = c(test1$smd$d),
                type = "c",
                smd_ci = "goulet")

  p2 = plot_smd(lambda = c(test2$smd$d_lambda),
                df = c(test2$smd$d_df),
                d = c(test2$smd$d),
                type = "cd",
                smd_ci = "goulet")

  p2 = plot_smd(lambda = c(test2$smd$d_lambda),
                df = c(test2$smd$d_df),
                d = c(test2$smd$d),
                type = "c",
                smd_ci = "goulet")


  expect_error(plot_smd(df = c(test1$smd$d_df),
                        SE = c(test1$smd$d_sigma)))



})

test_that("plot generic function",{
  set.seed(1812)

  samp1 = rnorm(25)
  samp2 = rnorm(25)

  df_samp = data.frame(y = c(samp1,samp2),
                       group = c(rep("g1",25),
                                 rep("g2",25)))

  test1 = t_TOST(x = samp1,
                 y = samp2,
                 paired = TRUE,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 smd_ci = "g")

  expect_error(plot(wilcox_TOST(x = samp1,
                                y = samp2,
                                paired = TRUE,
                                low_eqbound = -.5,
                                high_eqbound = .5)))

  p1 = plot(test1,
            type = "cd",
            estimates = "raw")
  p2 = plot(test1,
            type = "c",
            estimates = "raw")

  p3 = plot(test1,
            type = "cd",
            estimates = "SMD")
  p4 = plot(test1,
            type = "c",
            estimates = "SMD")

  p5 = suppressMessages(plot(test1,
                           type = "tnull",
                           estimates = "SMD"))
  p6 = suppressMessages(plot(test1,
                           type = "tnull"))
  p7 = plot(test1,
            type = "tnull",
            estimates = "raw")


})

test_that("Ensure paired output correct", {

  test1 = tsum_TOST(n1 = 23,
                   n2 = 23,
                   m2 = 14.2,
                   m1 = 13.8,
                   sd1 = 1.23,
                   sd2 = 1.78,
                   r12 = .41,
                   low_eqbound = -.5,
                   high_eqbound = .5,
                   paired = T,
                   bias_correction = FALSE,
                   eqbound_type = "raw",
                   rm_correction = T)

  expect_equal(sign(test1$effsize$estimate[1]),sign(test1$effsize$estimate[2]))
  test2 = tsum_TOST(n1 = 23,
                   n2 = 23,
                   m1 = 14.2,
                   m2 = 13.8,
                   sd2 = 1.23,
                   sd1 = 1.78,
                   r12 = .41,
                   low_eqbound = -.5,
                   high_eqbound = .5,
                   paired = T,
                   bias_correction = FALSE,
                   eqbound_type = "raw",
                   rm_correction = T)
  expect_equal(sign(test2$effsize$estimate[1]),sign(test2$effsize$estimate[2]))

  # Use vector form for paired tests, not formula
  test3 = t_TOST(x = sleep$extra[sleep$group == 1],
         y = sleep$extra[sleep$group == 2],
         low_eqbound = -.5,
         high_eqbound = .5,
         paired = T,
         bias_correction = FALSE,
         eqbound_type = "raw",
         rm_correction = T)

  expect_equal(sign(test3$effsize$estimate[1]),sign(test3$effsize$estimate[2]))

  set.seed(90183560)
  x1 = rnorm(30)
  y1 = rnorm(30)

  test4 = t_TOST(x=x1, y=y1,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 paired = T,
                 bias_correction = FALSE,
                 eqbound_type = "raw",
                 rm_correction = T)

  test5 = t_TOST(x=x1, y=y1,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 paired = T,
                 bias_correction = FALSE,
                 eqbound_type = "raw",
                 rm_correction = F)

  expect_equal(test4$effsize$estimate[2], .0952,
               tolerance = .001)

  expect_equal(test5$effsize$estimate[2], .0694,
               tolerance = .001)

  # mean(x1)
  # sd(x1)
  # mean(y1)
  # sd(y1)
  # x1: .14 (1.16)
  # x2: .04 (1.02)
  # r12 = .06

  # Use vector form for paired tests
  test4 = t_TOST(x = sleep$extra[sleep$group == 1],
                 y = sleep$extra[sleep$group == 2],
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 paired = T,
                 bias_correction = FALSE,
                 eqbound_type = "raw",
                 rm_correction = T)
  expect_equal(sign(test4$effsize$estimate[1]),sign(test4$effsize$estimate[2]))


})

test_that("Check NCT CIs for paired",{
  #effectsize::hedges_g(x=subset(sleep, group ==1)$extra,y=subset(sleep, group ==2)$extra,paired = TRUE, ci = .9)
  test1 = t_TOST(x=subset(sleep, group ==1)$extra,y=subset(sleep, group ==2)$extra,
               paired = TRUE,
               eqb = .5,
               smd_ci = "nct",
               bias_correction = T,
               glass = NULL)
  expect_equal(test1$effsize$estimate[2],
               -1.174, tolerance = .001)
  expect_equal(test1$effsize$lower.ci[2],
               -1.805, tolerance = .001)
  expect_equal(test1$effsize$upper.ci[2],
               -0.4977, tolerance = .001)
  #effectsize::cohens_d(x=subset(sleep, group ==1)$extra,y=subset(sleep, group ==2)$extra,paired = TRUE, ci = .9)

  test2= t_TOST(x=subset(sleep, group ==1)$extra,y=subset(sleep, group ==2)$extra,
               paired = TRUE,
               eqb = .5,
               smd_ci = "n",
               bias_correction = F,
               glass = NULL)
  test_plot = plot(test2)
  expect_equal(test2$effsize$estimate[2],
               -1.285, tolerance = .001)
  expect_equal(test2$effsize$lower.ci[2],
               -1.975, tolerance = .001)
  expect_equal(test2$effsize$upper.ci[2],
               -0.545, tolerance = .001)

  test3 = t_TOST(x=subset(sleep, group ==1)$extra,y=subset(sleep, group ==2)$extra,
                paired = TRUE,
                eqb = .5,
                smd_ci = "t",
                bias_correction = F,
                glass = NULL)
  p3 = plot(test3)
  test4 = t_TOST(x=subset(sleep, group ==1)$extra,y=subset(sleep, group ==2)$extra,
                 paired = TRUE,
                 eqb = .5,
                 smd_ci = "z",
                 bias_correction = F,
                 glass = NULL)
  p4 = plot(test4)
  test5 = t_TOST(x=subset(sleep, group ==1)$extra,y=subset(sleep, group ==2)$extra,
                 paired = TRUE,
                 eqb = .5,
                 smd_ci = "g",
                 bias_correction = F,
                 glass = "glass1")
  test5_smd = smd_calc(x=subset(sleep, group ==1)$extra,
                   y=subset(sleep, group ==2)$extra,
                 paired = TRUE,
                 smd_ci = "g",
                 bias_correction = F,
                 glass = "glass1",
                 alpha = .1)
  test5_smd = smd_calc(x=subset(sleep, group ==1)$extra,
                       y=subset(sleep, group ==2)$extra,
                       paired = TRUE,
                       smd_ci = "g",
                       bias_correction = F,
                       glass = "glass2",
                       alpha = .1)

  test6 = t_TOST(x=subset(sleep, group ==1)$extra,y=subset(sleep, group ==2)$extra,
                 paired = TRUE,
                 eqb = .5,
                 smd_ci = "z",
                 rm_correction = TRUE,
                 bias_correction = F)
  test7 = t_TOST(x=subset(sleep, group ==1)$extra,y=subset(sleep, group ==2)$extra,
                 paired = TRUE,
                 eqb = .5,
                 smd_ci = "n",
                 rm_correction = TRUE,
                 bias_correction = F)
  test8 = t_TOST(x=subset(sleep, group ==1)$extra,y=subset(sleep, group ==2)$extra,
                 paired = TRUE,
                 eqb = .5,
                 smd_ci = "t",
                 rm_correction = TRUE,
                 bias_correction = F)
})

test_that("Check NCT CIs for ind",{
  #effectsize::hedges_g(x=subset(sleep, group ==1)$extra,y=subset(sleep, group ==2)$extra,paired = FALSE, ci = .9, pooled_sd =TRUE)
  test1 = t_TOST(x=subset(sleep, group ==1)$extra,y=subset(sleep, group ==2)$extra,
                 paired = FALSE,
                 var.equal= TRUE,
                 eqb = .5,
                 smd_ci = "nct",
                 bias_correction = TRUE)
  expect_equal(test1$effsize$estimate[2],
               -0.7969352, tolerance = .01)
  expect_equal(test1$effsize$lower.ci[2],
               -1.523594, tolerance = .001)
  expect_equal(test1$effsize$upper.ci[2],
               -0.04942503, tolerance = .001)
  #effectsize::cohens_d(x=subset(sleep, group ==1)$extra,y=subset(sleep, group ==2)$extra,paired = FALSE, ci = .9)

  test2= t_TOST(x=subset(sleep, group ==1)$extra,y=subset(sleep, group ==2)$extra,
                paired = FALSE,
                var.equal=TRUE,
                eqb = .5,
                smd_ci = "n",
                bias_correction = FALSE,
                glass = NULL)
  test_plot = plot(test2)
  expect_equal(test2$effsize$estimate[2],
               -0.8321, tolerance = .001)
  expect_equal(test2$effsize$lower.ci[2],
               -1.5909, tolerance = .001)
  expect_equal(test2$effsize$upper.ci[2],
               -0.05161, tolerance = .001)

  test3 = t_TOST(x=subset(sleep, group ==1)$extra,y=subset(sleep, group ==2)$extra,
                 paired = FALSE,
                 eqb = .5,
                 smd_ci = "t",
                 bias_correction = F,
                 glass = NULL)
  p3 = plot(test3)
  test4 = t_TOST(x=subset(sleep, group ==1)$extra,y=subset(sleep, group ==2)$extra,
                 paired = FALSE,
                 eqb = .5,
                 smd_ci = "z",
                 bias_correction = F,
                 glass = NULL)
  p4 = plot(test4)
  test5 = t_TOST(x=subset(sleep, group ==1)$extra,y=subset(sleep, group ==2)$extra,
                 paired = FALSE,
                 eqb = .5,
                 smd_ci = "g",
                 bias_correction = F,
                 glass = "glass1")

  test6 = t_TOST(x=subset(sleep, group ==1)$extra,y=subset(sleep, group ==2)$extra,
                 paired = FALSE,
                 eqb = .5,
                 smd_ci = "z",
                 rm_correction = TRUE,
                 bias_correction = F)
  test7 = t_TOST(x=subset(sleep, group ==1)$extra,y=subset(sleep, group ==2)$extra,
                 paired = FALSE,
                 eqb = .5,
                 smd_ci = "n",
                 rm_correction = TRUE,
                 bias_correction = F)
  test8 = t_TOST(x=subset(sleep, group ==1)$extra,y=subset(sleep, group ==2)$extra,
                 paired = FALSE,
                 eqb = .5,
                 smd_ci = "t",
                 rm_correction = TRUE,
                 bias_correction = F)

  ## Check avg

  #t1=effectsize::hedges_g(x=subset(sleep, group ==1)$extra,y=subset(sleep, group ==2)$extra,paired = FALSE, ci = .9, pooled_sd =FALSE)
  test1 = t_TOST(x=subset(sleep, group ==1)$extra,y=subset(sleep, group ==2)$extra,
                 paired = FALSE,
                 var.equal= TRUE,
                 eqb = .5,
                 smd_ci = "nct",
                 bias_correction = TRUE)
  expect_equal(test1$effsize$estimate[2],
               -0.7969352, tolerance = .01)
  expect_equal(test1$effsize$lower.ci[2],
               -1.523594, tolerance = .001)
  expect_equal(test1$effsize$upper.ci[2],
               -0.04942503, tolerance = .001)
  #t2=effectsize::cohens_d(x=subset(sleep, group ==1)$extra,y=subset(sleep, group ==2)$extra,paired = FALSE, ci = .9, pooled_sd = FALSE)

  test2= t_TOST(x=subset(sleep, group ==1)$extra,y=subset(sleep, group ==2)$extra,
                paired = FALSE,
                var.equal=FALSE,
                eqb = .5,
                smd_ci = "n",
                bias_correction = FALSE,
                glass = NULL)
  test_plot = plot(test2)
  expect_equal(test2$effsize$estimate[2],
               -0.83218, tolerance = .001)
  expect_equal(test2$effsize$lower.ci[2],
               -1.59126, tolerance = .001)
  expect_equal(test2$effsize$upper.ci[2],
               -0.051069, tolerance = .001)
})

test_that("Check NCT CIs for one",{
  z1=subset(sleep, group ==1)$extra-subset(sleep, group ==2)$extra
  #t1=effectsize::hedges_g(x=z1,paired = FALSE, ci = .9, pooled_sd =TRUE)

  test1 = t_TOST(x=z1,
                 eqb = .5,
                 smd_ci = "nct",
                 bias_correction = TRUE)
  expect_equal(test1$effsize$estimate[2],
               -1.173925, tolerance = .01)
  expect_equal(test1$effsize$lower.ci[2],
               -1.804551, tolerance = .001)
  expect_equal(test1$effsize$upper.ci[2],
               -0.4977325, tolerance = .001)
  #t2 = effectsize::cohens_d(x=z1,paired = FALSE, ci = .9)

  test2= t_TOST(x=z1,
                eqb = .5,
                smd_ci = "nct",
                bias_correction = FALSE)
  test_plot = plot(test2)
  expect_equal(test2$effsize$estimate[2],
               -1.284558, tolerance = .001)
  expect_equal(test2$effsize$lower.ci[2],
               -1.974615, tolerance = .001)
  expect_equal(test2$effsize$upper.ci[2],
               -0.5446397, tolerance = .001)

  test3 = t_TOST(x=z1,
                 eqb = .5,
                 smd_ci = "t",
                 bias_correction = F,
                 glass = NULL)
  p3 = plot(test3)
  test4 = t_TOST(x=z1,
                 eqb = .5,
                 smd_ci = "z",
                 bias_correction = F,
                 glass = NULL)
  p4 = plot(test4)
  test5 = t_TOST(x=z1,
                 eqb = .5,
                 smd_ci = "g",
                 bias_correction = F,
                 glass = "glass1")

  test6 = t_TOST(x=z1,
                 eqb = .5,
                 smd_ci = "z",
                 rm_correction = TRUE,
                 bias_correction = F)
  test7 = t_TOST(x=z1,
                 eqb = .5,
                 smd_ci = "n",
                 rm_correction = TRUE,
                 bias_correction = F)
  test8 = t_TOST(x=z1,
                 eqb = .5,
                 smd_ci = "t",
                 rm_correction = TRUE,
                 bias_correction = F)

})

test_that("More tsum_test",{
  expect_error(TOSTER:::tsum_test(
    m1 = 12,
    sd1 = 1,
    n1 = 30,
    m2 = 11,
    sd2 = 1.5,
    n2 = 30,
    r12 = NULL,
    paired = FALSE,
    alternative = "two.sided",
    mu = c(0,1),
    var.equal = FALSE,
    conf.level = 0.95
  ))

  expect_error(TOSTER:::tsum_test(
    m1 = 12,
    sd1 = 1,
    n1 = 30,
    m2 = 11,
    sd2 = 1.5,
    n2 = 30,
    r12 = NULL,
    paired = FALSE,
    alternative = "two.sided",
    mu = 0,
    var.equal = FALSE,
    conf.level = 55
  ))

})

test_that("Formula methods reject paired = TRUE", {
  # Test that all formula methods properly reject paired = TRUE
  # to match base R behavior. paired = FALSE is allowed (redundant but harmless)
  data(sleep)
  
  # t_TOST.formula should reject paired = TRUE
  expect_error(
    t_TOST(extra ~ group, data = sleep, paired = TRUE, eqb = 1),
    "cannot use 'paired' in formula method"
  )
  
  # t_TOST.formula should allow paired = FALSE (redundant but harmless)
  expect_no_error(
    t_TOST(extra ~ group, data = sleep, paired = FALSE, eqb = 1)
  )
  
  # boot_t_TOST.formula should reject paired = TRUE
  expect_error(
    boot_t_TOST(extra ~ group, data = sleep, paired = TRUE, eqb = 1, R = 10),
    "cannot use 'paired' in formula method"
  )
  
  # boot_t_test.formula should reject paired = TRUE
  expect_error(
    boot_t_test(extra ~ group, data = sleep, paired = TRUE, R = 10),
    "cannot use 'paired' in formula method"
  )
  
  # wilcox_TOST.formula should reject paired = TRUE
  expect_error(
    wilcox_TOST(extra ~ group, data = sleep, paired = TRUE, eqb = 1),
    "cannot use 'paired' in formula method"
  )
  
  # simple_htest.formula should reject paired = TRUE
  expect_error(
    simple_htest(extra ~ group, data = sleep, paired = TRUE),
    "cannot use 'paired' in formula method"
  )
  
  # brunner_munzel.formula should reject paired = TRUE
  expect_error(
    brunner_munzel(extra ~ group, data = sleep, paired = TRUE),
    "cannot use 'paired' in formula method"
  )
  
  # Verify formula methods still work without paired parameter
  expect_no_error(
    t_TOST(extra ~ group, data = sleep, eqb = 1)
  )
  
  # Test calc functions also reject paired = TRUE with formula
  expect_error(
    smd_calc(extra ~ group, data = sleep, paired = TRUE),
    "cannot use 'paired' in formula method"
  )
  
  expect_error(
    ses_calc(extra ~ group, data = sleep, paired = TRUE),
    "cannot use 'paired' in formula method"
  )
  
  expect_error(
    boot_smd_calc(extra ~ group, data = sleep, paired = TRUE, R = 10),
    "cannot use 'paired' in formula method"
  )
  
  expect_error(
    boot_ses_calc(extra ~ group, data = sleep, paired = TRUE, R = 10),
    "cannot use 'paired' in formula method"
  )
  
})
