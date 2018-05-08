context("Do dataTOST and TOST functions give identical results")
library("TOSTER")

data <- read.csv("https://raw.githubusercontent.com/jasp-stats/jasp-desktop/development/Resources/Data%20Sets/Big%205%20(Dolan%2C%20Oort%2C%20Stoel%20%26%20Wicherts%2C%202009).csv", sep="")
#data <- read.csv("C:/Users/Daniel/Downloads/Big 5 (Dolan, Oort, Stoel & Wicherts, 2009).csv", sep="")

test_that("p-values for TOSTr are identical using dataTOSTr and TOSTr", {
  eq1_data <- dataTOSTr(data = data,
            pairs = list( list( i1="Extraversion", i2="Agreeableness")),
            low_eqbound_r = -0.15,
            high_eqbound_r = 0.15,
            desc = FALSE,
            plots = FALSE)

  eq1_sum <- TOSTr(n = length(data$Extraversion),  r = cor(data$Extraversion,data$Agreeableness),  # Correlation effect size
        low_eqbound_r = -0.15,  # Value for the lower equivalence bound
        high_eqbound_r = 0.15,  # Value for the higher equivalence bound
        alpha = 0.05, # Alpha level for TOST and NHST
        plot = FALSE)

  #Store dataTOST results as dataframe
  eq1_data <- as.data.frame(eq1_data$tost)

  #Check if p-values are equal
  expect_equal(as.numeric(eq1_data[8]), eq1_sum$TOST_p2)
  expect_equal(as.numeric(eq1_data[11]), eq1_sum$TOST_p1)

  #flip Extraversion variable, so correlation becomes negative

  data$Extraversion <- -1*data$Extraversion

  eq1_data <- dataTOSTr(data = data,
                        pairs = list( list( i1="Extraversion", i2="Agreeableness")),
                        low_eqbound_r = -0.15,
                        high_eqbound_r = 0.15,
                        desc = FALSE,
                        plots = TRUE)

  eq1_sum <- TOSTr(n = length(data$Extraversion),  r = cor(data$Extraversion,data$Agreeableness),  # Correlation effect size
                   low_eqbound_r = -0.15,  # Value for the lower equivalence bound
                   high_eqbound_r = 0.15,  # Value for the higher equivalence bound
                   alpha = 0.05, # Alpha level for TOST and NHST
                   plot = FALSE)

  #Store dataTOST results as dataframe
  eq1_data <- as.data.frame(eq1_data$tost)

  #Check if p-values are equal
  expect_equal(as.numeric(eq1_data[8]), eq1_sum$TOST_p2)
  expect_equal(as.numeric(eq1_data[11]), eq1_sum$TOST_p1)
})
