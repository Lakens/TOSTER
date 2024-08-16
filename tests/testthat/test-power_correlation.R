hush = function(code) {
  sink("NUL") # use /dev/null in UNIX
  tmp = code
  sink()
  return(tmp)
}

test_that("test messages and warnings", {

  test1 = power_z_cor(n = 90,
                      rho = .25,
                      power = NULL,
                      null = 0,
                      alpha = .05,
                      alternative = "t")

  test1 = power_z_cor(n = 90,
                      rho = .25,
                      power = NULL,
                      null = 0,
                      alpha = .05,
                      alternative = "l")

  test1 = power_z_cor(n = 90,
                      rho = .25,
                      power = NULL,
                      null = 0,
                      alpha = .05,
                      alternative = "g")

  test_e = power_z_cor(n = 90,
                       rho = 0,
                       power = NULL,
                       null = 0.25,
                       alpha = .05,
                       alternative = "e")

  expect_error(power_z_cor(n = NULL,
                           rho = .25,
                           power = NULL,
                           null = 0,
                           alpha = .05,
                           alternative = "t"))

  expect_error(power_z_cor(n = 90,
                           rho = .25,
                           power = NULL,
                           null = 0,
                           alpha = 1.05,
                           alternative = "t"))
  expect_error(power_z_cor(n = 90,
                           rho = .25,
                           power = 1.1,
                           null = 0,
                           alpha = NULL,
                           alternative = "t"))
  expect_error(power_z_cor(n = 3,
                           rho = .25,
                           power = NULL,
                           null = 0,
                           alpha = .05,
                           alternative = "t"))

  test1_n = power_z_cor(n = 90,
                      rho = NULL,
                      power = .8,
                      null = 0,
                      alpha = .05,
                      alternative = "t")
  test1_n = power_z_cor(n = 90,
                        rho = NULL,
                        power = .8,
                        null = 0,
                        alpha = .05,
                        alternative = "g")

  test1_alpha = power_z_cor(n = 90,
                        rho = .4,
                        power = .8,
                        null = 0,
                        alpha = NULL,
                        alternative = "t")

  test_e = power_z_cor(n = NULL,
                       rho = 0,
                       power = .8,
                       null = 0.25,
                       alpha = .05,
                       alternative = "e")

  test_e = power_z_cor(n = 90,
                       rho = 0,
                       power = .8,
                       null = 0.25,
                       alpha = NULL,
                       alternative = "e")

  expect_error(power_z_cor(
    n = 90,
    rho = 0,
    power = NULL,
    null = NULL,
    alpha = .05,
    alternative = "e"
  ))

  expect_error(power_z_cor(
    n = 90,
    rho = 0,
    power = NULL,
    null = NULL,
    alpha = NULL,
    alternative = "e"
  ))


  expect_error(power_z_cor(
    n = 3,
    rho = 0,
    power = NULL,
    null = .25,
    alpha = .05,
    alternative = "e"
  ))


})


test_that("Matching known results",{
  # From PASS

  test1 = power_z_cor(n = 20,
              rho = .3,
              power = NULL,
              null = 0,
              alpha = .01,
              alternative = "t")

  expect_equal(abs(test1$power-.094),0, tolerance = .01)

  test2 = power_z_cor(n = 10, rho = 0, power = NULL,
                      null = .001, alpha = .05,
                      alternative = "equivalence")

  expect_equal(test2$power,0)
})
