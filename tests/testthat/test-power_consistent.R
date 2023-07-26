#context("Test if sensitivity, post-hoc, and a-priori power are internally consistent")

#library("TOSTER")

test_that("power functions are internally consistent", {
  hush = function(code) {
    sink("NUL") # use /dev/null in UNIX
    tmp = code
    sink()
    return(tmp)
  }
  #withr::local_options(lifecycle_verbosity = "quiet")

  suppressWarnings({ suppressMessages(hush({
  ## tests for one-sample t-test
  pow_n <-powerTOSTone(alpha=0.05, statistical_power=0.9,
                       low_eqbound_d=-0.3, high_eqbound_d=0.3)
  expect_equal(powerTOSTone(alpha=0.05, N=pow_n,
                            low_eqbound_d=-0.3, high_eqbound_d=0.3),
               0.9, tolerance = 0.001)
  expect_equal(powerTOSTone(alpha=0.05, N=pow_n,
                            statistical_power=0.9)[2],
               0.3, tolerance = 0.001)

  pow_power <- powerTOSTone(alpha=0.05, N=50,
                            low_eqbound_d=-0.3, high_eqbound_d=0.3)
  pow_power1 = powerTOSTone.raw(alpha=0.05, N=50,
                                low_eqbound=-0.3, high_eqbound=0.3,sd=1)
  #pow_power2 = powerTOSTone.raw(alpha=0.05, statistical_power = 0.3664239, low_eqbound=-0.3, high_eqbound=0.3,sd=1)
  #pow_power3 = powerTOSTone.raw(alpha=0.05, N=50, statistical_power = 0.3664239,sd=1)

  expect_equal(powerTOSTone(alpha=0.05, statistical_power = pow_power,
                            low_eqbound_d=-0.3, high_eqbound_d=0.3),
               50, tolerance = 0.001)
  expect_equal(powerTOSTone(alpha=0.05, statistical_power = pow_power,
                            N=50)[2], 0.30, tolerance = 0.001)

  pow_bound <- powerTOSTone(alpha=0.05, N=50, statistical_power = 0.9)[2]
  expect_equal(powerTOSTone(alpha=0.05, statistical_power = .9,
                            low_eqbound_d=-pow_bound, high_eqbound_d=pow_bound),
               50, tolerance = 0.001)
  expect_equal(powerTOSTone(alpha=0.05, N = 50,
                            low_eqbound_d=-pow_bound,
                            high_eqbound_d=pow_bound),
               .9, tolerance = 0.001)

  ## tests for two-sample independent t-test
  pow_n <- powerTOSTtwo(alpha=0.05, statistical_power=0.9,
                        low_eqbound_d=-0.3, high_eqbound_d=0.3)
  expect_equal(powerTOSTtwo(alpha=0.05, N=pow_n,
                            low_eqbound_d=-0.3, high_eqbound_d=0.3),
               0.9, tolerance = 0.001)
  expect_equal(powerTOSTtwo(alpha=0.05, N=pow_n,
                            statistical_power=0.9)[2],
               0.3, tolerance = 0.001)

  pow_power <- powerTOSTtwo(alpha=0.05, N=200,
                            low_eqbound_d=-0.3, high_eqbound_d=0.3)
  expect_equal(powerTOSTtwo(alpha=0.05, statistical_power = pow_power,
                            low_eqbound_d=-0.3, high_eqbound_d=0.3),
               200, tolerance = 0.001)
  expect_equal(powerTOSTtwo(alpha=0.05, statistical_power = pow_power,
                            N=200)[2], 0.30, tolerance = 0.001)

  pow_bound <- powerTOSTtwo(alpha=0.05, N=50, statistical_power = 0.9)[2]
  expect_equal(powerTOSTtwo(alpha=0.05, statistical_power = .9,
                            low_eqbound_d=-pow_bound,
                            high_eqbound_d=pow_bound),
               50, tolerance = 0.001)
  expect_equal(powerTOSTtwo(alpha=0.05, N = 50,
                            low_eqbound_d=-pow_bound,
                            high_eqbound_d=pow_bound),
               .9, tolerance = 0.001)

  ## tests for two-sample paired t-test
  pow_n <- powerTOSTpaired(alpha=0.05, statistical_power=0.9,
                           low_eqbound_dz=-0.3, high_eqbound_dz=0.3)
  test2 = powerTOSTpaired(alpha=0.05, N=pow_n,
                  low_eqbound_dz=-0.3,
                  high_eqbound_dz=0.3)
  expect_equal(test2, 0.9, tolerance = 0.001)

  test3 = powerTOSTpaired(alpha=0.05, N=pow_n,
                  statistical_power=0.9)[2]
  expect_equal(test3, 0.3, tolerance = 0.001)

  pow_power <- powerTOSTpaired(alpha=0.05, N=200,
                               low_eqbound_dz=-0.3,
                               high_eqbound_dz=0.3)

  test3 = powerTOSTpaired(alpha=0.05,
                          statistical_power = pow_power,
                          N=200)[2]
  expect_equal(test3, 0.30, tolerance = 0.001)

  pow_bound <- powerTOSTpaired(alpha=0.05, N=50,
                               statistical_power = 0.9)[2]
  test2 =powerTOSTpaired(alpha=0.05, statistical_power = .9,
                         low_eqbound_dz=-pow_bound,
                         high_eqbound_dz=pow_bound)
  expect_equal(test2, 50, tolerance = 0.001)
  test3 =powerTOSTpaired(alpha=0.05, N = 50,
                         low_eqbound_dz=-pow_bound,
                         high_eqbound_dz=pow_bound)
  expect_equal(test3, .9, tolerance = 0.001)

  ## tests for correlations
  pow_n <- powerTOSTr(alpha=0.05, statistical_power=0.9,
                      low_eqbound_r=-0.3, high_eqbound_r=0.3)
  expect_equal(powerTOSTr(alpha=0.05, N=pow_n,
                          low_eqbound_r=-0.3, high_eqbound_r=0.3),
               0.9, tolerance = 0.001)

  expect_equal(powerTOSTr(alpha=0.05, N=pow_n,
                          statistical_power=0.9)[2],
               0.3, tolerance = 0.001)

  pow_power <- powerTOSTr(alpha=0.05, N=100, low_eqbound_r=-0.3,
                          high_eqbound_r=0.3)
  expect_equal(powerTOSTr(alpha=0.05, statistical_power = pow_power,
                          low_eqbound_r=-0.3, high_eqbound_r=0.3),
               100, tolerance = 0.001)
  expect_equal(powerTOSTr(alpha=0.05, statistical_power = pow_power,
                          N=100)[2],
               0.30, tolerance = 0.001)

  pow_bound <- powerTOSTr(alpha=0.05, N=50, statistical_power = 0.9)[2]
  expect_equal(powerTOSTr(alpha=0.05, statistical_power = .9,
                          low_eqbound_r=-pow_bound,
                          high_eqbound_r=pow_bound),
               50, tolerance = 0.001)
  expect_equal(powerTOSTr(alpha=0.05, N = 50,
                          low_eqbound_r=-pow_bound, high_eqbound_r=pow_bound),
               .9, tolerance = 0.001)

  expect_equal(powerTOSTone.raw(alpha=0.05, statistical_power=0.9,
                                sd = 0.5, low_eqbound=-0.25, high_eqbound=0.25),
               powerTOSTone(alpha=0.05, statistical_power=0.9,
                            low_eqbound_d=-0.5, high_eqbound_d=0.5))
  expect_equal(powerTOSTtwo.raw(alpha=0.05, N = 20,
                                low_eqbound=-200, high_eqbound=200, sdpooled=200),
               powerTOSTtwo(alpha=0.05, N = 20,
                            low_eqbound_d = 1, high_eqbound_d = 1))

  expect_equal(powerTOSTpaired.raw(alpha = 0.05,statistical_power = 0.8,
                                   low_eqbound = -3, high_eqbound = 3,
                                   sdif = 10),
               powerTOSTpaired(alpha = 0.05, statistical_power = 0.8,
                               low_eqbound_dz = -0.3, high_eqbound_dz = 0.3))
  #
  }))
  })
})

