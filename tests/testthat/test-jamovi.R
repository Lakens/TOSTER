#library(testthat)
test_that("dataTOSTtwo: internal consistency",{
  data(sleep)
  t1 = t_TOST(formula = extra ~ group,
              data = sleep,
              low_eqbound = -2,
              high_eqbound = 2,
              smd_ci = "goulet")

  t2 = dataTOSTtwo(data = sleep,
                   deps = "extra",
                   group = "group",
                   low_eqbound = -2,
                   high_eqbound = 2)

  expect_equal(t1$TOST$p.value,
                    c(t2$tost$asDF$`p[0]`,
                    t2$tost$asDF$`p[1]`,
                    t2$tost$asDF$`p[2]`),
               ignore_attr = TRUE
                    )

  expect_equal(t1$effsize$estimate,
               c(t2$effsize$asDF$`est[raw]`,
                 t2$effsize$asDF$`est[cohen]`),
               ignore_attr = TRUE
  )

  # var.eqaul = TRUE
  t1 = t_TOST(formula = extra ~ group,
              data = sleep,
              var.equal = TRUE,
              low_eqbound = -2,
              high_eqbound = 2,
              smd_ci = "goulet")

  t2 = dataTOSTtwo(data = sleep,
                   deps = "extra",
                   group = "group",
                   var_equal = TRUE,
                   low_eqbound = -2,
                   high_eqbound = 2)

  expect_equal(t1$TOST$p.value,
               c(t2$tost$asDF$`p[0]`,
                 t2$tost$asDF$`p[1]`,
                 t2$tost$asDF$`p[2]`),
               ignore_attr = TRUE
  )

  expect_equal(t1$effsize$estimate,
               c(t2$effsize$asDF$`est[raw]`,
                 t2$effsize$asDF$`est[cohen]`),
               ignore_attr = TRUE
  )

  # bound type to SMD
  t1 = suppressMessages(t_TOST(formula = extra ~ group,
              data = sleep,
              var.equal = TRUE,
              #hypothesis = "MET",
              eqbound_type = "SMD",
              low_eqbound = -2,
              high_eqbound = 2,
              smd_ci = "goulet"))

  t2 = suppressMessages(dataTOSTtwo(data = sleep,
                   deps = "extra",
                   group = "group",
                   var_equal = TRUE,
                   eqbound_type = "SMD",
                   #hypothesis = "MET",
                   low_eqbound = -2,
                   high_eqbound = 2))

  expect_equal(t1$TOST$p.value,
               c(t2$tost$asDF$`p[0]`,
                 t2$tost$asDF$`p[1]`,
                 t2$tost$asDF$`p[2]`),
               ignore_attr = TRUE
  )

  expect_equal(t1$effsize$estimate,
               c(t2$effsize$asDF$`est[raw]`,
                 t2$effsize$asDF$`est[cohen]`),
               ignore_attr = TRUE
  )

  # hypothesis test to MET

  # bound type to SMD
  t1 = suppressMessages(t_TOST(formula = extra ~ group,
              data = sleep,
              var.equal = TRUE,
              hypothesis = "MET",
              eqbound_type = "SMD",
              low_eqbound = -2,
              high_eqbound = 2,
              smd_ci = "goulet"))

  t2 = suppressMessages(dataTOSTtwo(data = sleep,
                   deps = "extra",
                   group = "group",
                   var_equal = TRUE,
                   eqbound_type = "SMD",
                   hypothesis = "MET",
                   low_eqbound = -2,
                   high_eqbound = 2))

  expect_equal(t1$TOST$p.value,
               c(t2$tost$asDF$`p[0]`,
                 t2$tost$asDF$`p[1]`,
                 t2$tost$asDF$`p[2]`),
               ignore_attr = TRUE
  )

  expect_equal(t1$effsize$estimate,
               c(t2$effsize$asDF$`est[raw]`,
                 t2$effsize$asDF$`est[cohen]`),
               ignore_attr = TRUE
  )

  # bound type to SMD
  t1 = t_TOST(formula = extra ~ group,
              data = sleep,
              var.equal = TRUE,
              hypothesis = "MET",
              eqbound_type = "raw",
              bias_correction = FALSE,
              low_eqbound = -2,
              high_eqbound = 2,
              smd_ci = "goulet")

  t2 = dataTOSTtwo(data = sleep,
                   deps = "extra",
                   group = "group",
                   var_equal = TRUE,
                   eqbound_type = "raw",
                   hypothesis = "MET",
                   smd_type = "d",
                   low_eqbound = -2,
                   high_eqbound = 2)

  expect_equal(t1$TOST$p.value,
               c(t2$tost$asDF$`p[0]`,
                 t2$tost$asDF$`p[1]`,
                 t2$tost$asDF$`p[2]`),
               ignore_attr = TRUE
  )

  expect_equal(t1$effsize$estimate,
               c(t2$effsize$asDF$`est[raw]`,
                 t2$effsize$asDF$`est[cohen]`),
               ignore_attr = TRUE
  )
})


test_that("dataTOSTpaired: internal consistency",{
  data(sleep)

  pair1 = subset(sleep, group == 2)$extra
  pair2 = subset(sleep, group == 1)$extra
  sleep2 = data.frame(pair1 = pair1,
                      pair2 = pair2)

  t1 = t_TOST(formula = extra ~ group,
              data = sleep,
              paired = TRUE,
              low_eqbound = -2,
              high_eqbound = 2,
              smd_ci = "goulet")

  t2 = dataTOSTpaired(data = sleep2,
                      pair1 = "pair1",
                      pair2 = "pair2",
                      low_eqbound = -2,
                      high_eqbound = 2)

  expect_equal(t1$TOST$p.value,
               t2$tost$asDF$p,
               ignore_attr = TRUE
  )

  expect_equal(t1$effsize$estimate,
               t2$effsize$asDF$est,
               ignore_attr = TRUE
  )



  # bound type to SMD
  t1 = suppressMessages(t_TOST(formula = extra ~ group,
              data = sleep,
              paired = TRUE,
              #hypothesis = "MET",
              eqbound_type = "SMD",
              low_eqbound = -2,
              high_eqbound = 2,
              smd_ci = "goulet"))

  t2 = suppressMessages(dataTOSTpaired(data = sleep2,
                      pair1 = "pair1",
                      pair2 = "pair2",
                   eqbound_type = "SMD",
                   #hypothesis = "MET",
                   low_eqbound = -2,
                   high_eqbound = 2))

  expect_equal(t1$TOST$p.value,
               t2$tost$asDF$p,
               ignore_attr = TRUE
  )

  expect_equal(t1$effsize$estimate,
               t2$effsize$asDF$est,
               ignore_attr = TRUE
  )

  # Test MET

  # bound type to SMD
  t1 = t_TOST(formula = extra ~ group,
              data = sleep,
              paired = TRUE,
              hypothesis = "MET",
              #eqbound_type = "SMD",
              low_eqbound = -2,
              high_eqbound = 2,
              smd_ci = "goulet")

  t2 = dataTOSTpaired(data = sleep2,
                      pair1 = "pair1",
                      pair2 = "pair2",
                      #eqbound_type = "SMD",
                      hypothesis = "MET",
                      low_eqbound = -2,
                      high_eqbound = 2)

  expect_equal(t1$TOST$p.value,
               t2$tost$asDF$p,
               ignore_attr = TRUE
  )

  expect_equal(t1$effsize$estimate,
               t2$effsize$asDF$est,
               ignore_attr = TRUE
  )

  # bound type to SMD
  t1 = suppressMessages(t_TOST(formula = extra ~ group,
              data = sleep,
              paired = TRUE,
              #hypothesis = "MET",
              eqbound_type = "SMD",
              bias_correction = FALSE,
              low_eqbound = -2,
              high_eqbound = 2,
              smd_ci = "goulet"))

  t2 = suppressMessages(dataTOSTpaired(data = sleep2,
                      pair1 = "pair1",
                      pair2 = "pair2",
                      eqbound_type = "SMD",
                      smd_type = "d",
                      #hypothesis = "MET",
                      low_eqbound = -2,
                      high_eqbound = 2))

  expect_equal(t1$TOST$p.value,
               t2$tost$asDF$p,
               ignore_attr = TRUE
  )

  expect_equal(t1$effsize$estimate,
               t2$effsize$asDF$est,
               ignore_attr = TRUE
  )

  t2 = suppressMessages(dataTOSTpaired(data = sleep2,
                      pair1 = "pair1",
                      pair2 = "pair2",
                      eqbound_type = "SMD",
                      smd_type = "d",
                      #hypothesis = "MET",
                      low_eqbound = -2,
                      high_eqbound = 2,
                      desc = TRUE,
                      plots = TRUE,
                      indplot = TRUE,
                      diffplot = TRUE))
})

test_that("dataTOSTone: internal consistency",{
  data(sleep)

  t1 = t_TOST(x = sleep$extra,
              low_eqbound = -2,
              high_eqbound = 2,
              smd_ci = "goulet")

  t2 = dataTOSTone(data = sleep,
                   vars = "extra",
                   eqbound_type = "raw",
                   low_eqbound = -2,
                   high_eqbound = 2)

  expect_equal(t1$TOST$p.value,
               c(t2$tost$asDF$`p[0]`,
                 t2$tost$asDF$`p[1]`,
                 t2$tost$asDF$`p[2]`),
               ignore_attr = TRUE
  )

  expect_equal(t1$effsize$estimate,
               c(t2$effsize$asDF$`est[raw]`,
                 t2$effsize$asDF$`est[cohen]`),
               ignore_attr = TRUE
  )



  # bound type to SMD
  t1 = suppressMessages(t_TOST(x = sleep$extra,
              low_eqbound = -2,
              high_eqbound = 2,
              eqbound_type = "SMD",
              smd_ci = "goulet"))

  t2 = suppressMessages(dataTOSTone(data = sleep,
                   vars = "extra",
                   low_eqbound = -2,
                   high_eqbound = 2,
                   eqbound_type = "SMD"))

  expect_equal(t1$TOST$p.value,
               c(t2$tost$asDF$`p[0]`,
                 t2$tost$asDF$`p[1]`,
                 t2$tost$asDF$`p[2]`),
               ignore_attr = TRUE
  )

  expect_equal(t1$effsize$estimate,
               c(t2$effsize$asDF$`est[raw]`,
                 t2$effsize$asDF$`est[cohen]`),
               ignore_attr = TRUE
  )

  # hypothesis MET
  t1 = suppressMessages(t_TOST(x = sleep$extra,
              low_eqbound = -2,
              high_eqbound = 2,
              eqbound_type = "SMD",
              hypothesis = "MET",
              smd_ci = "goulet"))

  t2 = suppressMessages(dataTOSTone(data = sleep,
                      vars = extra,
                      eqbound_type = "SMD",
                      hypothesis = "MET",
                      low_eqbound = -2,
                      high_eqbound = 2))

  expect_equal(t1$TOST$p.value,
               c(t2$tost$asDF$`p[0]`,
                 t2$tost$asDF$`p[1]`,
                 t2$tost$asDF$`p[2]`),
               ignore_attr = TRUE
  )

  expect_equal(t1$effsize$estimate,
               c(t2$effsize$asDF$`est[raw]`,
                 t2$effsize$asDF$`est[cohen]`),
               ignore_attr = TRUE
  )

  # No bias correction

  t1 = suppressMessages(t_TOST(x = sleep$extra,
              low_eqbound = -2,
              high_eqbound = 2,
              eqbound_type = "SMD",
              hypothesis = "MET",
              bias_correction = FALSE,
              smd_ci = "goulet"))

  t2 = suppressMessages(dataTOSTone(data = sleep,
                   vars = extra,
                   eqbound_type = "SMD",
                   hypothesis = "MET",
                   low_eqbound = -2,
                   high_eqbound = 2,
                   smd_type = "d"))

  expect_equal(t1$TOST$p.value,
               c(t2$tost$asDF$`p[0]`,
                 t2$tost$asDF$`p[1]`,
                 t2$tost$asDF$`p[2]`),
               ignore_attr = TRUE
  )

  expect_equal(t1$effsize$estimate,
               c(t2$effsize$asDF$`est[raw]`,
                 t2$effsize$asDF$`est[cohen]`),
               ignore_attr = TRUE
  )
})

test_that("dataTOSTr",{
  data('iris')

  t1 = dataTOSTr(
    data = iris,
    pairs = list(
      list(
        i1="Sepal.Length",
        i2="Sepal.Width")),
    plot = TRUE
  )

  t2 = cor(iris$Sepal.Length,iris$Sepal.Width)

  expect_equal(t1$tost$asDF$r,
               t2)

  t1 = dataTOSTr(
    data = iris,
    pairs = list(
      list(
        i1="Sepal.Length",
        i2="Sepal.Width")),
    cor_type = "spearman",
    plot = TRUE
  )

  t2 = cor(iris$Sepal.Length,iris$Sepal.Width,
           method = "spearman")

  expect_equal(t1$tost$asDF$r,
               t2)

  t1 = dataTOSTr(
    data = iris,
    pairs = list(
      list(
        i1="Sepal.Length",
        i2="Sepal.Width")),
    cor_type = "kendall"
  )

  t2 = cor(iris$Sepal.Length,iris$Sepal.Width,
           method = "kendall")

  expect_equal(t1$tost$asDF$r,
               t2)

  t1_MET = dataTOSTr(
    data = iris,
    pairs = list(
      list(
        i1="Sepal.Length",
        i2="Sepal.Width")),
    hypothesis = "MET",
    cor_type = "kendall",
    plot = TRUE
  )

  t2_MET = cor(iris$Sepal.Length,iris$Sepal.Width,
           method = "kendall")

  expect_equal(t1_MET$tost$asDF$r,
               t2_MET)

  expect_equal(t1_MET$tost$asDF$p,
               t1$tost$asDF$p)


})

test_that("datatosttwoprop tests",{
  set.seed(8020)
  d1 = rbinom(10,1,.5)
  d2 = rbinom(10,1,.25)
  df1 = data.frame(
    d1 = d1,
    d2 = d2
  )

  t1 = datatosttwoprop(
    data = df1,
    var = d1,
    level = "0",
    group = d2)
  t2 = datatosttwoprop(
    data = df1,
    var = d1,
    level = "0",
    group = d2,
    plot = TRUE)

  expect_equal(t1$tost$asDF$p,
               t2$tost$asDF$p)

  t1 = datatosttwoprop(
    data = df1,
    var = d1,
    level = "0",
    hypothesis = "MET",
    group = d2)
  t2 = datatosttwoprop(
    data = df1,
    var = d1,
    hypothesis = "MET",
    level = "0",
    group = d2,
    plot = TRUE)


  expect_equal(t1$tost$asDF$p,
               t2$tost$asDF$p)

  # Run plot function

  # Both
  t1 = plot_cor(r = .3, n =22, method = "pearson")
  t2 = plot_cor(r = .43, n =32, method = "spearman")
  t3 = plot_cor(r = .21, n =100, method = "kendall")

  # density
  t1 = plot_cor(r = .3, n =22, method = "pearson",
                type = "cd")
  t2 = plot_cor(r = .43, n =32, method = "spearman",
                type = "cd")
  t3 = plot_cor(r = .21, n =100, method = "kendall",
                type = "cd")
  # consonance
  t1 = plot_cor(r = .3, n =22, method = "pearson",
                type = "c")
  t2 = plot_cor(r = .43, n =32, method = "spearman",
                type = "c")
  t3 = plot_cor(r = .21, n =100, method = "kendall",
                type = "c")
})

test_that("plot functions for jamovi work",{
  #skip_on_cran()

  set.seed(8020)
  d1 = rbinom(10,1,.5)
  d2 = rbinom(10,1,.25)
  df1 = data.frame(
    d1 = d1,
    d2 = d2
  )

  #proportions
  t1 = datatosttwoprop(
    data = df1,
    var = d1,
    level = "0",
    hypothesis = "MET",
    group = d2)

  t2 = datatosttwoprop(
    data = df1,
    var = d1,
    level = "0",
    group = d2,
    plot = TRUE)

  expect_equal(t1$tost$asDF$p,
               t2$tost$asDF$p)

  test = t2$plot
  test$plot

  #TOSTpaired

  data(sleep)

  pair1 = subset(sleep, group == 2)$extra
  pair2 = subset(sleep, group == 1)$extra
  sleep2 = data.frame(pair1 = pair1,
                      pair2 = pair2)

  t2 = dataTOSTpaired(data = sleep2,
                      pair1 = "pair1",
                      pair2 = "pair2",
                      low_eqbound = -2,
                      high_eqbound = 2,
                      desc = TRUE,
                      plots = TRUE,
                      indplot = TRUE,
                      diffplot = TRUE)

  t2$plots
  t2$indplot
  t2$diffplot
  t2$desc

  # Correlation plot

  data('iris')

  t1 = dataTOSTr(
    data = iris,
    pairs = list(
      list(
        i1="Sepal.Length",
        i2="Sepal.Width")),
    plot = TRUE
  )

  t1$plots

  # one sample

  t2 = dataTOSTone(data = sleep,
                   vars = "extra",
                   eqbound_type = "raw",
                   low_eqbound = -2,
                   high_eqbound = 2,
                   desc =TRUE,
                   plots = TRUE)

  t2$plots
  t2$desc

  # two sample

  t2 = dataTOSTtwo(data = sleep,
                   deps = "extra",
                   group = "group",
                   low_eqbound = -2,
                   high_eqbound = 2,
                   desc = TRUE,
                   plots = TRUE,
                   descplots = TRUE)

  t2$plots
  t2$descplots

})
