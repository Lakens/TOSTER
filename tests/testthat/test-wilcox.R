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

  # Normal one sample ----

  test1 = wilcox_TOST(x = samp1,
                 low_eqbound = -.5,
                 high_eqbound = .5)


  test3 = wilcox_TOST(x = samp1,
                 low_eqbound = -.5,
                 high_eqbound = .5,
                 hypothesis = "MET")


  expect_equal(1-test1$TOST$p.value[2],
               test3$TOST$p.value[2],
               tolerance = .001)

  expect_equal(1-test1$TOST$p.value[3],
               test3$TOST$p.value[3],
               tolerance = .001)

  prtest = hush(print(test3))

})

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


  test3 = wilcox_TOST(x = samp1,
                 y = samp2,
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
