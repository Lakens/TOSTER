#context("Do dataTOST and TOST functions give identical results")
#library("TOSTER")

#data <- read.csv("https://raw.githubusercontent.com/jasp-stats/jasp-desktop/development/Resources/Data%20Sets/Big%205%20(Dolan%2C%20Oort%2C%20Stoel%20%26%20Wicherts%2C%202009).csv", sep="")
#data <- read.csv("C:/Users/Daniel/Downloads/Big 5 (Dolan, Oort, Stoel & Wicherts, 2009).csv", sep="")

test_that("p-values for TOSTr are identical using dataTOSTr and TOSTr", {

  hush = function(code) {
    sink("NUL") # use /dev/null in UNIX
    tmp = code
    sink()
    return(tmp)
  }

  data("cars")
  eq1_data <- suppressMessages(hush(dataTOSTr(
    data = cars,
    pairs = list(list(i1 = "speed", i2 = "dist")),
    low_eqbound_r = -0.15,
    high_eqbound_r = 0.15,
    desc = FALSE,
    plots = FALSE
  ) ))

  eq1_sum <- suppressMessages(hush(
    TOSTr(
      n = length(cars$dist),
      r = cor(cars$dist, cars$speed),
      # Correlation effect size
      low_eqbound_r = -0.15,
      # Value for the lower equivalence bound
      high_eqbound_r = 0.15,
      # Value for the higher equivalence bound
      alpha = 0.05,
      # Alpha level for TOST and NHST
      plot = FALSE
    )
  ))

  #Store dataTOST results as dataframe
  eq1_data <- as.data.frame(eq1_data$tost)

  #Check if p-values are equal
  expect_equal(as.numeric(eq1_data[4]), eq1_sum$NHST_p)
  expect_equal(as.numeric(eq1_data$cil),eq1_sum$LL_CI_TOST)
  expect_equal(as.numeric(eq1_data$ciu),eq1_sum$UL_CI_TOST)
  #expect_equal(as.numeric(eq1_data[11]), eq1_sum$TOST_p1)

  #flip dist
  cars$dist = cars$dist*-1
  eq1_data <- suppressMessages(hush(
    dataTOSTr(
    data = cars,
    pairs = list(list(i1 = "speed", i2 = "dist")),
    low_eqbound_r = -0.15,
    high_eqbound_r = 0.15,
    desc = FALSE,
    plots = FALSE
  )
  ))

  eq1_sum <- suppressMessages(hush(
    TOSTr(
      n = length(cars$dist),
      r = cor(cars$dist, cars$speed),
      # Correlation effect size
      low_eqbound_r = -0.15,
      # Value for the lower equivalence bound
      high_eqbound_r = 0.15,
      # Value for the higher equivalence bound
      alpha = 0.05,
      # Alpha level for TOST and NHST
      plot = FALSE
    )
  ))

  #Store dataTOST results as dataframe
  eq1_data <- as.data.frame(eq1_data$tost)

  #Check if p-values are equal
  expect_equal(as.numeric(eq1_data[4]), eq1_sum$NHST_p)
  expect_equal(as.numeric(eq1_data$cil),eq1_sum$LL_CI_TOST)
  expect_equal(as.numeric(eq1_data$ciu),eq1_sum$UL_CI_TOST)

  eq1_data <- suppressMessages(hush( dataTOSTr(
    data = cars,
    pairs = list(list(i1 = "speed", i2 = "dist")),
    low_eqbound_r = -0.15,
    high_eqbound_r = 0.15,
    desc = TRUE,
    plots = TRUE
  )
  ))
})
