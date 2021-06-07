
test_that("NCSS TOST 198-1 to 197-8",{


  yieldA = c(452, 874, 554, 447,
             356, 754, 558, 574,
             664, 682, 547, 435,
             245)

  yieldB =   c(546, 547, 774, 465,
               459, 665, 467, 365,
               589, 534, 456, 651,
               654,
               665, 546,
               537)


  test1 = t_TOST(x = yieldA,
                 y = yieldB,
                 low_eqbound = -110,
                 high_eqbound =110,
                 var.equal = TRUE,
                 alpha = .025)

  expect_equal(27, test1$TOST$df[1])

  expect_equal(51.114285, test1$TOST$SE[1])

  expect_equal(-8.1153846, test1$effsize$estimate[1])

  expect_equal(-112.99323,
               round(test1$effsize$lower.ci[1], 5))

  expect_equal(96.76246,
               round(test1$effsize$upper.ci[1], 5))

  test2 = t_TOST(x = yieldA,
                 y = yieldB,
                 low_eqbound = -110,
                 high_eqbound =110,
                 var.equal = FALSE,
                 alpha = .025)

  expect_equal(19.169, round(test2$TOST$df[1],3))

  expect_equal(53.618547, test2$TOST$SE[1])

  expect_equal(-8.1153846, test2$effsize$estimate[1])

  expect_equal(-120.27336,
               round(test2$effsize$lower.ci[1], 5))

  expect_equal(104.04259,
               round(test2$effsize$upper.ci[1], 5))

})
