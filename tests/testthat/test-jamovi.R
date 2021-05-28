#library(testthat)
test_that("dataTOSTtwo: internal consistency",{
  data(sleep)
  t1 = t_TOST(formula = extra ~ group,
              data = sleep,
              low_eqbound = -2,
              high_eqbound = 2)

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
              high_eqbound = 2)

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
  t1 = t_TOST(formula = extra ~ group,
              data = sleep,
              var.equal = TRUE,
              #hypothesis = "MET",
              eqbound_type = "SMD",
              low_eqbound = -2,
              high_eqbound = 2)

  t2 = dataTOSTtwo(data = sleep,
                   deps = "extra",
                   group = "group",
                   var_equal = TRUE,
                   eqbound_type = "SMD",
                   #hypothesis = "MET",
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

  # hypothesis test to MET

  # bound type to SMD
  t1 = t_TOST(formula = extra ~ group,
              data = sleep,
              var.equal = TRUE,
              hypothesis = "MET",
              eqbound_type = "SMD",
              low_eqbound = -2,
              high_eqbound = 2)

  t2 = dataTOSTtwo(data = sleep,
                   deps = "extra",
                   group = "group",
                   var_equal = TRUE,
                   eqbound_type = "SMD",
                   hypothesis = "MET",
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
  t1 = t_TOST(formula = extra ~ group,
              data = sleep,
              var.equal = TRUE,
              hypothesis = "MET",
              eqbound_type = "raw",
              bias_correction = FALSE,
              low_eqbound = -2,
              high_eqbound = 2)

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
              high_eqbound = 2)

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
  t1 = t_TOST(formula = extra ~ group,
              data = sleep,
              paired = TRUE,
              #hypothesis = "MET",
              eqbound_type = "SMD",
              low_eqbound = -2,
              high_eqbound = 2)

  t2 = dataTOSTpaired(data = sleep2,
                      pair1 = "pair1",
                      pair2 = "pair2",
                   eqbound_type = "SMD",
                   #hypothesis = "MET",
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

  # Test MET

  # bound type to SMD
  t1 = t_TOST(formula = extra ~ group,
              data = sleep,
              paired = TRUE,
              hypothesis = "MET",
              #eqbound_type = "SMD",
              low_eqbound = -2,
              high_eqbound = 2)

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
  t1 = t_TOST(formula = extra ~ group,
              data = sleep,
              paired = TRUE,
              #hypothesis = "MET",
              eqbound_type = "SMD",
              bias_correction = FALSE,
              low_eqbound = -2,
              high_eqbound = 2)

  t2 = dataTOSTpaired(data = sleep2,
                      pair1 = "pair1",
                      pair2 = "pair2",
                      eqbound_type = "SMD",
                      smd_type = "d",
                      #hypothesis = "MET",
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
})

test_that("dataTOSTone: internal consistency",{
  data(sleep)

  t1 = t_TOST(x = sleep$extra,
              low_eqbound = -2,
              high_eqbound = 2)

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
  t1 = t_TOST(x = sleep$extra,
              low_eqbound = -2,
              high_eqbound = 2,
              eqbound_type = "SMD")

  t2 = dataTOSTone(data = sleep,
                   vars = "extra",
                   low_eqbound = -2,
                   high_eqbound = 2,
                   eqbound_type = "SMD")

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
  t1 = t_TOST(x = sleep$extra,
              low_eqbound = -2,
              high_eqbound = 2,
              eqbound_type = "SMD",
              hypothesis = "MET")

  t2 = dataTOSTone(data = sleep,
                      vars = extra,
                      eqbound_type = "SMD",
                      hypothesis = "MET",
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

  # No bias correction

  t1 = t_TOST(x = sleep$extra,
              low_eqbound = -2,
              high_eqbound = 2,
              eqbound_type = "SMD",
              hypothesis = "MET",
              bias_correction = FALSE)

  t2 = dataTOSTone(data = sleep,
                   vars = extra,
                   eqbound_type = "SMD",
                   hypothesis = "MET",
                   low_eqbound = -2,
                   high_eqbound = 2,
                   smd_type = "d")

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
