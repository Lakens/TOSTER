
#library(testthat)

test_that("PASS ex#1 pg 460-8", {
  # Solve for power
  t1 = power_t_TOST(
    n = 3,
    delta = -4,
    sd = 18,
    low_eqbound = -19.2,
    high_eqbound = 19.2,
    alpha = 0.05,
    type = "two.sample"
  )
  expect_equal(round(t1$power,5),0.0386)
  # solve n
  t1 = power_t_TOST(
    power = 0.0386,
    delta = -4,
    sd = 18,
    low_eqbound = -19.2,
    high_eqbound = 19.2,
    alpha = 0.05,
    type = "two.sample"
  )
  expect_equal(round(t1$n,0),3)

  # solve alpha

  t1 = power_t_TOST(
    power = 0.0386,
    n = 3,
    delta = -4,
    sd = 18,
    low_eqbound = -19.2,
    high_eqbound = 19.2,
    type = "two.sample"
  )
  expect_equal(round(t1$alpha,2),.05)


  t2 = power_t_TOST(
    n = 8,
    delta = -4,
    sd = 18,
    low_eqbound = -19.2,
    high_eqbound = 19.2,
    alpha = 0.05,
    type = "two.sample"
  )
  expect_equal(round(t2$power,5),0.28871)

  t3 = power_t_TOST(
    n = 15,
    delta = -4,
    sd = 18,
    low_eqbound = -19.2,
    high_eqbound = 19.2,
    alpha = 0.05,
    type = "two.sample"
  )
  expect_equal(round(t3$power,5),0.69339)

  t4 = power_t_TOST(
    n = 30,
    delta = -4,
    sd = 18,
    low_eqbound = -19.2,
    high_eqbound = 19.2,
    alpha = 0.05,
    type = "two.sample"
  )
  expect_equal(round(t4$power,5),0.94326)

  t5 = power_t_TOST(
    n = 50,
    delta = -4,
    sd = 18,
    low_eqbound = -19.2,
    high_eqbound = 19.2,
    alpha = 0.05,
    type = "two.sample"
  )
  expect_equal(round(t5$power,5),0.99458)

})

test_that("PASS ex#2 pg 460-10", {
  # Solve for N
  t1 = power_t_TOST(
    power = .8,
    delta = -4,
    sd = 18,
    low_eqbound = -19.2,
    high_eqbound = 19.2,
    alpha = 0.05,
    type = "two.sample"
  )
  expect_equal(round(t1$n,0),19)



})

test_that("PASS ex#3 pg 460-11", {
  # Solve for N
  t1 = power_t_TOST(
    power = .9,
    delta = 0,
    sd = 100,
    low_eqbound = -10,
    high_eqbound = 10,
    alpha = 0.025,
    type = "two.sample"
  )
  expect_equal(round(t1$n,0),2600)

  t2 = power_t_TOST(
    power = .9,
    delta = 2,
    sd = 100,
    low_eqbound = -10,
    high_eqbound = 10,
    alpha = 0.025,
    type = "two.sample"
  )
  expect_equal(round(t2$n,0),3305)

})

test_that("PASS ex#4 pg 460-12", {
  # Solve for N
  t1 = power_t_TOST(
    power = .8,
    delta = -2,
    sd = 8,
    low_eqbound = -5,
    high_eqbound = 5,
    alpha = 0.05,
    type = "two.sample"
  )
  expect_equal(round(t1$n,0),89)
})

test_that("PASS ex#5 pg 519-5", {
  # Solve for power
  t1 = power_t_TOST(
    n = 5,
    delta = -4,
    sd = 25.4558,
    low_eqbound = -19.2,
    high_eqbound = 19.2,
    alpha = 0.05,
    type = "one.sample"
  )

  expect_equal(round(t1$power,5),0.10599)

  t2 = power_t_TOST(
    n = 15,
    delta = -4,
    sd = 25.4558,
    low_eqbound = -19.2,
    high_eqbound = 19.2,
    alpha = 0.05,
    type = "paired"
  )
  expect_equal(round(t2$power,5),0.66629)

  t3 = power_t_TOST(
    n = 15,
    delta = -4,
    sd = 25.4558,
    low_eqbound = -19.2,
    high_eqbound = 19.2,
    alpha = 0.05,
    type = "paired"
  )
  expect_equal(round(t3$power,5),0.69339)

  t4 = power_t_TOST(
    n = 30,
    delta = -4,
    sd = 25.4558,
    low_eqbound = -19.2,
    high_eqbound = 19.2,
    alpha = 0.05,
    type = "paired"
  )
  expect_equal(round(t4$power,5),0.93855)

  t5 = power_t_TOST(
    n = 50,
    delta = -4,
    sd = 25.4558,
    low_eqbound = -19.2,
    high_eqbound = 19.2,
    alpha = 0.05,
    type = "paired"
  )
  expect_equal(round(t5$power,5),0.99410)

})

test_that("PASS ex#6 pg 519-5", {
  # Solve for power
  t1 = power_t_TOST(
    power = .8,
    delta = 0,
    sd = .1,
    low_eqbound = -.05,
    high_eqbound = .05,
    alpha = 0.05,
    type = "one.sample"
  )

  expect_equal(round(t1$n,0),36)
  })

test_that("PASS ex#7 pg 519-8", {
  # Solve for power
  t1 = power_t_TOST(
    power = .7,
    delta = 0,
    sd = 28.284,
    low_eqbound = -20,
    high_eqbound = 20,
    alpha = 0.05,
    type = "paired"
  )

  expect_equal(ceiling(t1$n),16)

  t2 = power_t_TOST(
    power = .7,
    delta = -5,
    sd = 28.284,
    low_eqbound = -20,
    high_eqbound = 20,
    alpha = 0.05,
    type = "paired"
  )

  expect_equal(ceiling(t2$n),20)

  t3 = power_t_TOST(
    power = .7,
    delta = -10,
    sd = 28.284,
    low_eqbound = -20,
    high_eqbound = 20,
    alpha = 0.05,
    type = "paired"
  )

  expect_equal(ceiling(t3$n),40)

  t4 = power_t_TOST(
    power = .7,
    delta = -15,
    sd = 28.284,
    low_eqbound = -20,
    high_eqbound = 20,
    alpha = 0.05,
    type = "paired"
  )

  expect_equal(ceiling(t4$n),152)

  t4 = power_t_TOST(
    n = 152,
    delta = -15,
    sd = 28.284,
    low_eqbound = -20,
    high_eqbound = 20,
    alpha = 0.05,
    type = "paired"
  )

  expect_equal(round(t4$power,1),.7)

  t4 = power_t_TOST(
    n = 152,
    power = .7,
    delta = -15,
    sd = 28.284,
    low_eqbound = -20,
    high_eqbound = 20,
    type = "paired"
  )

  expect_equal(round(t4$alpha,2),.05)
})

test_that("errors",{

  expect_error(  t4 = power_t_TOST(
    n = 152,
    power = .7,
    delta = -15,
    sd = 28.284,
    low_eqbound = -20,
    high_eqbound = 20,
    type = "paired",
    alpha = .05
  )
  )

  expect_error(  t4 = power_t_TOST(
    n = 152,
    delta = -15,
    sd = 28.284,
    low_eqbound = -20,
    high_eqbound = 20,
    type = "paired",
    alpha = 1
  )
  )

  expect_error(  t4 = power_t_TOST(
    n = 152,
    delta = -15,
    sd = 28.284,
    type = "paired",
    alpha = .05
  )
  )

})

