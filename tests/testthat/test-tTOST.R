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

  # Normal one sample ----

  test1 = t_TOST(x = samp1,
                 low_eqbound = -.5,
                 high_eqbound = .5)

  test2 = suppressMessages(t_TOST(x = samp1,
                                low_eqbound = -.5,
                                high_eqbound = .5,
                                eqbound_type = "SMD"))

  test3 = t_TOST(x = samp1,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 hypothesis = "MET")

  test4 = suppressMessages( { t_TOST(x = samp1,
                                   low_eqbound = -.5,
                                   high_eqbound = .5,
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
                                   bias_correction = FALSE)
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
  p1 = plot(test4)
  p2 = plot(test4,
            type = "c")

})


test_that("Run examples for two sample", {

  hush = function(code) {
    sink("NUL") # use /dev/null in UNIX
    tmp = code
    sink()
    return(tmp)
  }


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

  test2 = suppressMessages( t_TOST(formula = y ~ group,
                                 data = df_samp,
                                 var.equal = TRUE,
                                 low_eqbound = -.5,
                                 high_eqbound = .5,
                                 eqbound_type = "SMD",
                                 bias_correction = FALSE) )

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
                                 bias_correction = FALSE) )

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

  # Run with formula
  test1 = t_TOST(formula = y ~ group,
                 data = df_samp,
                 paired = TRUE,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 bias_correction = FALSE)

  test2 = suppressMessages( t_TOST(formula = y ~ group,
                                 data = df_samp,
                                 paired = TRUE,
                                 low_eqbound = -.5,
                                 high_eqbound = .5,
                                 eqbound_type = "SMD",
                                 bias_correction = FALSE) )

  test3 = t_TOST(formula = y ~ group,
                 data = df_samp,
                 paired = TRUE,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 hypothesis = "MET",
                 bias_correction = FALSE)

  test4 = suppressMessages(  t_TOST(formula = y ~ group,
                                  data = df_samp,
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
                type = "cd")

  p1 = plot_smd(lambda = c(test1$smd$d_lambda),
                df = c(test1$smd$d_df),
                d = c(test1$smd$d),
                type = "c")

  p2 = plot_smd(lambda = c(test2$smd$d_lambda),
                df = c(test2$smd$d_df),
                d = c(test2$smd$d),
                type = "cd")

  p2 = plot_smd(lambda = c(test2$smd$d_lambda),
                df = c(test2$smd$d_df),
                d = c(test2$smd$d),
                type = "c")


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
                 high_eqbound = .5)

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
