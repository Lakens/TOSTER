context("Run Examples for t_TOST")

# need hush function to run print through examples

hush = function(code) {
  sink("NUL") # use /dev/null in UNIX
  tmp = code
  sink()
  return(tmp)
}

test_that("Run examples for one sample", {

  set.seed(3164964)

  samp1 = rnorm(33)

  expect_error(t_TOST())

  test1 = t_TOST(x = samp1,
                 low_eqbound = -.5,
                 high_eqbound = .5)

  test2 = t_TOST(x = samp1,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 eqbound_type = "SMD")

  test3 = t_TOST(x = samp1,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 hypothesis = "MET")

  test4 = t_TOST(x = samp1,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 eqbound_type = "SMD",
                 hypothesis = "MET")

  expect_equal(1-test1$TOST$p.value[2],
               test3$TOST$p.value[2])

  expect_equal(1-test1$TOST$p.value[3],
               test3$TOST$p.value[3])

  expect_equal(1-test2$TOST$p.value[2],
               test4$TOST$p.value[2])

  expect_equal(1-test2$TOST$p.value[3],
               test4$TOST$p.value[3])


  # Re-run with bias correction not run
  test1 = t_TOST(x = samp1,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 bias_correction = FALSE)

  test2 = t_TOST(x = samp1,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 eqbound_type = "SMD",
                 bias_correction = FALSE)

  test3 = t_TOST(x = samp1,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 hypothesis = "MET",
                 bias_correction = FALSE)

  test4 = t_TOST(x = samp1,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 eqbound_type = "SMD",
                 hypothesis = "MET",
                 bias_correction = FALSE)

  expect_equal(1-test1$TOST$p.value[2],
               test3$TOST$p.value[2])

  expect_equal(1-test1$TOST$p.value[3],
               test3$TOST$p.value[3])

  expect_equal(1-test2$TOST$p.value[2],
               test4$TOST$p.value[2])

  expect_equal(1-test2$TOST$p.value[3],
               test4$TOST$p.value[3])

  prtest = hush(print(test4))
  p1 = plot(test4)

})


test_that("Run examples for two sample", {

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

  test2 = t_TOST(x = samp1,
                 y = samp2,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 eqbound_type = "SMD")

  test3 = t_TOST(x = samp1,
                 y = samp2,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 hypothesis = "MET")

  test4 = t_TOST(x = samp1,
                 y = samp2,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 eqbound_type = "SMD",
                 hypothesis = "MET")

  expect_equal(1-test1$TOST$p.value[2],
               test3$TOST$p.value[2])

  expect_equal(1-test1$TOST$p.value[3],
               test3$TOST$p.value[3])

  expect_equal(1-test2$TOST$p.value[2],
               test4$TOST$p.value[2])

  expect_equal(1-test2$TOST$p.value[3],
               test4$TOST$p.value[3])

  # Re-run with bias correction not run and non-Welch
  test1 = t_TOST(x = samp1,
                 y = samp2,
                 var.equal = TRUE,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 bias_correction = FALSE)

  test2 = t_TOST(x = samp1,
                 y = samp2,
                 var.equal = TRUE,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 eqbound_type = "SMD",
                 bias_correction = FALSE)

  test3 = t_TOST(x = samp1,
                 y = samp2,
                 var.equal = TRUE,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 hypothesis = "MET",
                 bias_correction = FALSE)

  test4 = t_TOST(x = samp1,
                 y = samp2,
                 var.equal = TRUE,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 eqbound_type = "SMD",
                 hypothesis = "MET",
                 bias_correction = FALSE)

  expect_equal(1-test1$TOST$p.value[2],
               test3$TOST$p.value[2])

  expect_equal(1-test1$TOST$p.value[3],
               test3$TOST$p.value[3])

  expect_equal(1-test2$TOST$p.value[2],
               test4$TOST$p.value[2])

  expect_equal(1-test2$TOST$p.value[3],
               test4$TOST$p.value[3])

  # Run with formula
  test1 = t_TOST(formula = y ~ group,
                 data = df_samp,
                 var.equal = TRUE,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 bias_correction = FALSE)

  test2 = t_TOST(formula = y ~ group,
                 data = df_samp,
                 var.equal = TRUE,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 eqbound_type = "SMD",
                 bias_correction = FALSE)

  test3 = t_TOST(formula = y ~ group,
                 data = df_samp,
                 var.equal = TRUE,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 hypothesis = "MET",
                 bias_correction = FALSE)

  test4 = t_TOST(formula = y ~ group,
                 data = df_samp,
                 var.equal = TRUE,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 eqbound_type = "SMD",
                 hypothesis = "MET",
                 bias_correction = FALSE)

  expect_equal(1-test1$TOST$p.value[2],
               test3$TOST$p.value[2])

  expect_equal(1-test1$TOST$p.value[3],
               test3$TOST$p.value[3])

  expect_equal(1-test2$TOST$p.value[2],
               test4$TOST$p.value[2])

  expect_equal(1-test2$TOST$p.value[3],
               test4$TOST$p.value[3])

  prtest = hush(print(test4))
  p1 = plot(test4)

})


test_that("Run examples for paired samples", {

  set.seed(789461245)

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

  test2 = t_TOST(x = samp1,
                 y = samp2,
                 paired = TRUE,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 eqbound_type = "SMD")

  test3 = t_TOST(x = samp1,
                 y = samp2,
                 paired = TRUE,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 hypothesis = "MET")

  test4 = t_TOST(x = samp1,
                 y = samp2,
                 paired = TRUE,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 eqbound_type = "SMD",
                 hypothesis = "MET")

  expect_equal(1-test1$TOST$p.value[2],
               test3$TOST$p.value[2])

  expect_equal(1-test1$TOST$p.value[3],
               test3$TOST$p.value[3])

  expect_equal(1-test2$TOST$p.value[2],
               test4$TOST$p.value[2])

  expect_equal(1-test2$TOST$p.value[3],
               test4$TOST$p.value[3])

  # Re-run with bias correction not run and rm_correction
  # rm_correction = TRUE
  test1 = t_TOST(x = samp1,
                 y = samp2,
                 paired = TRUE,
                 rm_correction = TRUE,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 bias_correction = FALSE)

  test2 = t_TOST(x = samp1,
                 y = samp2,
                 paired = TRUE,
                 rm_correction = TRUE,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 eqbound_type = "SMD",
                 bias_correction = FALSE)

  test3 = t_TOST(x = samp1,
                 y = samp2,
                 paired = TRUE,
                 rm_correction = TRUE,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 hypothesis = "MET",
                 bias_correction = FALSE)

  test4 = t_TOST(x = samp1,
                 y = samp2,
                 paired = TRUE,
                 rm_correction = TRUE,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 eqbound_type = "SMD",
                 hypothesis = "MET",
                 bias_correction = FALSE)

  expect_equal(1-test1$TOST$p.value[2],
               test3$TOST$p.value[2])

  expect_equal(1-test1$TOST$p.value[3],
               test3$TOST$p.value[3])

  expect_equal(1-test2$TOST$p.value[2],
               test4$TOST$p.value[2])

  expect_equal(1-test2$TOST$p.value[3],
               test4$TOST$p.value[3])

  # Run with formula
  test1 = t_TOST(formula = y ~ group,
                 data = df_samp,
                 paired = TRUE,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 bias_correction = FALSE)

  test2 = t_TOST(formula = y ~ group,
                 data = df_samp,
                 paired = TRUE,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 eqbound_type = "SMD",
                 bias_correction = FALSE)

  test3 = t_TOST(formula = y ~ group,
                 data = df_samp,
                 paired = TRUE,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 hypothesis = "MET",
                 bias_correction = FALSE)

  test4 = t_TOST(formula = y ~ group,
                 data = df_samp,
                 paired = TRUE,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 eqbound_type = "SMD",
                 hypothesis = "MET",
                 bias_correction = FALSE)

  expect_equal(1-test1$TOST$p.value[2],
               test3$TOST$p.value[2])

  expect_equal(1-test1$TOST$p.value[3],
               test3$TOST$p.value[3])

  expect_equal(1-test2$TOST$p.value[2],
               test4$TOST$p.value[2])

  expect_equal(1-test2$TOST$p.value[3],
               test4$TOST$p.value[3])

  prtest = hush(print(test4))
  p1 = plot(test4)

})

test_that("Run examples for plot_smd", {

  set.seed(777641421991)

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

  test2 = t_TOST(x = samp1,
                 y = samp2,
                 paired = TRUE,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 eqbound_type = "SMD")

  test3 = t_TOST(x = samp1,
                 y = samp2,
                 paired = TRUE,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 hypothesis = "MET")

  test4 = t_TOST(x = samp1,
                 y = samp2,
                 paired = TRUE,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 eqbound_type = "SMD",
                 hypothesis = "MET")

  p1 = plot_smd(lambda = c(test1$smd$d_lambda,
                           test2$smd$d_lambda,
                           test3$smd$d_lambda,
                           test4$smd$d_lambda),
                df = c(test1$smd$d_df,
                       test2$smd$d_df,
                       test3$smd$d_df,
                       test4$smd$d_df),
                SE = c(test1$smd$d_sigma,
                       test2$smd$d_sigma,
                       test3$smd$d_sigma,
                       test4$smd$d_sigma),
                smd_label = c("Study #1",
                              "Study #2",
                              "Study #3",
                              "Study #4"))

  p2 = plot_smd(lambda = c(test1$smd$d_lambda),
                df = c(test1$smd$d_df),
                SE = c(test1$smd$d_sigma))

  expect_error(plot_smd(df = c(test1$smd$d_df),
                        SE = c(test1$smd$d_sigma)))



})
