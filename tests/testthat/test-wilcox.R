hush = function(code) {
  sink("NUL") # use /dev/null in UNIX
  tmp = code
  sink()
  return(tmp)
}

# Normal one sample ----

test_that("Run examples for one sample", {

  set.seed(3164964)

  samp1 = rnorm(33)

  expect_error(wilcox_TOST())

  expect_message({
   wilcox_TOST(x = samp1,
               hypothesis = "MET",
               eqb = .05)
  })

  test1 = wilcox_TOST(x = samp1,
                 low_eqbound = -.5,
                 high_eqbound = .5)
  # Test data.frame output for backward compatibility
  test1_ses_df  = ses_calc(x = samp1,
                           alpha = .1,
                           se_method = "fisher",
                           output = "data.frame")
  expect_equal(test1$effsize$estimate[2],
               test1_ses_df$estimate)
  expect_equal(test1$effsize$lower.ci[2],
               test1_ses_df$lower.ci)
  expect_equal(test1$effsize$upper.ci[2],
               test1_ses_df$upper.ci)

  # Test htest output (new default)
  test1_ses_htest = ses_calc(x = samp1, alpha = .1)
  expect_s3_class(test1_ses_htest, "htest")
  expect_equal(unname(test1_ses_htest$estimate), test1$effsize$estimate[2])

  test1_ses  = boot_ses_calc(x = samp1,
                        alpha = .1,
                        boot_ci = "s")
  test1_ses  = boot_ses_calc(x = samp1,
                        alpha = .1,
                        boot_ci = "p")
  test1_ses  = boot_ses_calc(x = samp1,
                        alpha = .1)
  ash = as_htest(test1)

  test3 = wilcox_TOST(x = samp1,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 hypothesis = "MET")
  ash = as_htest(test3)

  expect_equal(1-test1$TOST$p.value[2],
               test3$TOST$p.value[2],
               tolerance = .001)

  expect_equal(1-test1$TOST$p.value[3],
               test3$TOST$p.value[3],
               tolerance = .001)

  prtest = hush(print(test3))

})

# Two sample ----
test_that("Run examples for two sample", {

  set.seed(651466441)

  samp1 = rnorm(25)
  samp2 = rnorm(25)

  df_samp = data.frame(y = c(samp1,samp2),
                       group = c(rep("g1",25),
                                 rep("g2",25)))

  expect_error(wilcox_TOST())

  test1 = wilcox_TOST(x = samp1,
                 y = samp2,
                 low_eqbound = -.5,
                 high_eqbound = .5)
  ash = as_htest(test1)

  test3 = wilcox_TOST(x = samp1,
                 y = samp2,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 hypothesis = "MET")
  ash = as_htest(test3)
  test1 = wilcox_TOST(formula = y ~ group,
                      data = df_samp,
                      low_eqbound = -.5,
                      high_eqbound = .5)
  # Test with data.frame output for comparison with boot_ses_calc
  test1_smd = ses_calc(formula = y ~ group,
                       data = df_samp,
                       output = "data.frame")
  test1_smd_cstat = ses_calc(formula = y ~ group,
                             data = df_samp,
                             ses = "cstat",
                             output = "data.frame")
  test1_smd_odds = ses_calc(formula = y ~ group,
                            data = df_samp,
                            ses = "odds",
                            output = "data.frame")
  test1_smd_logodds = ses_calc(formula = y ~ group,
                               data = df_samp,
                               ses = "logodds",
                               output = "data.frame")
  test1_smd_boot = boot_ses_calc(formula = y ~ group,
                       data = df_samp,
                       R = 99)
  test1_smd_boot_cstat = boot_ses_calc(formula = y ~ group,
                                 data = df_samp,
                                 ses = "cstat",
                                 R = 99)
  test1_smd_boot_odds = boot_ses_calc(formula = y ~ group,
                                 data = df_samp,
                                 ses = "odds",
                                 R = 99)
  test1_smd_boot_logodds = boot_ses_calc(formula = y ~ group,
                                      data = df_samp,
                                      ses = "logodds",
                                      R = 99)
  # Point estimates should match between ses_calc and boot_ses_calc
  # boot_ses_calc returns htest with named estimate, ses_calc data.frame has unnamed
  expect_equal(unname(test1_smd_boot$estimate),
               test1_smd$estimate)
  expect_equal(unname(test1_smd_boot_cstat$estimate),
               test1_smd_cstat$estimate)
  expect_equal(unname(test1_smd_boot_odds$estimate),
               test1_smd_odds$estimate)
  expect_equal(unname(test1_smd_boot_logodds$estimate),
               test1_smd_logodds$estimate)

  # Test htest output (new default)
  test1_smd_htest = ses_calc(formula = y ~ group, data = df_samp)
  expect_s3_class(test1_smd_htest, "htest")

  expect_error(ses_calc(formula = y ~ group,
                        data = df_samp,
                        alpha = 1.1))

  test3 = wilcox_TOST(formula = y ~ group,
                      data = df_samp,
                      low_eqbound = -.5,
                      high_eqbound = .5,
                      hypothesis = "MET")

  expect_equal(1-test1$TOST$p.value[2],
               test3$TOST$p.value[2],
               tolerance = .003)

  expect_equal(1-test1$TOST$p.value[3],
               test3$TOST$p.value[3],
               tolerance = .003)
  prtest = hush(print(test3))

})


test_that("Run examples for paired samples", {

  set.seed(789461245)

  samp1 = rnorm(25)
  expect_error(wilcox_TOST(x = samp1,
                      eqb = c(-1,1,.5)))
  expect_error(wilcox_TOST(x = samp1,
                           eqb = c(-1,1),
                           alpha= 1.22))
  expect_error(wilcox_TOST(x = samp1,
                           eqb = 1,
                           hypothesis = "DDDDD"))
  expect_error(wilcox_TOST(Sepal.Width ~ Species, data = iris))
  samp2 = rnorm(25)

  cor12 = stats::cor(samp1,samp2)

  df_samp = data.frame(y = c(samp1,samp2),
                       group = c(rep("g1",25),
                                 rep("g2",25)))

  expect_error(wilcox_TOST())

  test1 = wilcox_TOST(x = samp1,
                 y = samp2,
                 paired = TRUE,
                 low_eqbound = -.5,
                 high_eqbound = .5)
  test1 = wilcox_TOST(x = samp1,
                      y = samp2,
                      paired = TRUE,
                      eqb = c(-.5, .5))
  test1 = wilcox_TOST(x = samp1,
                      y = samp2,
                      paired = TRUE,
                      eqb =  .5)

  test1_ses = boot_ses_calc(x = samp1,
                      y = samp2,
                      paired = TRUE)

  test3 = wilcox_TOST(x = samp1,
                 y = samp2,
                 paired = TRUE,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 hypothesis = "MET")

  expect_equal(1-test1$TOST$p.value[2],
               test3$TOST$p.value[2],
               tolerance = .005)

  expect_equal(1-test1$TOST$p.value[3],
               test3$TOST$p.value[3],
               tolerance = .005)

  prtest = hush(print(test1))

})


test_that("Check rbs",{
  set.seed(1847501)
  z1 = rnorm(35)
  z2 = rnorm(35)

  # Two Sample ------
  rbs1 = TOSTER:::rbs_calc(x=z1, y=z2, mu=0, paired=FALSE)
  rbs2 = TOSTER:::rbs_calc(x=z2, y=z1, mu=0, paired=FALSE)
  expect_equal(abs(rbs1),abs(rbs2))

  # Paired Sample ------
  rbs1 = TOSTER:::rbs_calc(x=z1, y=z2, mu=0, paired=TRUE)
  rbs2 = TOSTER:::rbs_calc(x=z2, y=z1, mu=0, paired=TRUE)
  expect_equal(abs(rbs1),abs(rbs2))

  # Two Sample ------
  rbs1 = TOSTER:::rbs_calc(x=z1, y=z2, mu=.5, paired=FALSE)
  rbs2 = TOSTER:::rbs_calc(x=z2, y=z1, mu=.5, paired=FALSE)
  #expect_equal(abs(rbs1),abs(rbs2))

  # Paired Sample ------
  rbs1 = TOSTER:::rbs_calc(x=z1, y=z2, mu=.5, paired=TRUE)
  rbs2 = TOSTER:::rbs_calc(x=z2, y=z1, mu=.5, paired=TRUE)
  #expect_equal(abs(rbs1),abs(rbs2))
  x1 = 1.25
  expect_equal(TOSTER:::pr_to_odds(TOSTER:::odds_to_pr(x1)),
               1.25)

  rb1 = .75

  expect_equal(TOSTER:::cstat_to_rb(TOSTER:::rb_to_cstat(rb1)),
               rb1)

  z1 = .45

  expect_equal(TOSTER:::rho_to_z(TOSTER:::z_to_rho(z1)),
               z1)
  test1= TOSTER:::ranktransform(z1,sign=FALSE,method="first")
  test1= TOSTER:::ranktransform(z1,sign=TRUE,method="first")
  test2 = TOSTER:::ranktransform(1)
  test25 = TOSTER:::ranktransform(c(7,7,7,7,7))
  test3 = TOSTER:::ranktransform(c(NA,NA,NA))
  test3 = TOSTER:::ranktransform(c(TRUE,TRUE,FALSE))
})
