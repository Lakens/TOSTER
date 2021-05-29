#context("Identical results for raw and SMD power functions")

#library("TOSTER")

test_that("power values for raw and standardized functions are identical", {
  hush = function(code) {
    sink("NUL") # use /dev/null in UNIX
    tmp = code
    sink()
    return(tmp)
  }

  hush({
    expect_equal(powerTOSTone.raw(alpha=0.05, statistical_power=0.9, sd = 0.5, low_eqbound=-0.25, high_eqbound=0.25),
                 powerTOSTone(alpha=0.05, statistical_power=0.9, low_eqbound_d=-0.5, high_eqbound_d=0.5))

    expect_equal(powerTOSTtwo.raw(alpha=0.05, N = 20, low_eqbound=-200, high_eqbound=200, sdpooled=200),
                 powerTOSTtwo(alpha=0.05, N = 20, low_eqbound_d = 1, high_eqbound_d = 1))

    expect_equal(powerTOSTpaired.raw(alpha = 0.05,statistical_power = 0.8,low_eqbound = -3, high_eqbound = 3, sdif = 10),
                 powerTOSTpaired(alpha = 0.05, statistical_power = 0.8, low_eqbound_dz = -0.3, high_eqbound_dz = 0.3))
  })


})
