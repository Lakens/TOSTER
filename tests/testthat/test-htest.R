# All htest extensions tested here

test_that("simple_htest: t-test & wilcox", {

  set.seed(3164964)

  expect_equal(t.test(1:10, y = c(7:20))$pvalue,
               simple_htest(1:10, y = c(7:20))$pvalue)
  expect_equal(t.test(1:10, y = c(7:20, 200))$pvalue,
               simple_htest(1:10, y = c(7:20, 200))$pvalue)

  testy1 = describe_htest(simple_htest(1:10, y = c(7:20, 200)))
  testy2 = describe_htest(t.test(1:10, y = c(7:20, 200)))
  expect_equivalent(testy1,testy2)

  testy1 = df_htest(simple_htest(1:10, y = c(7:20, 200)))
  testy2 = df_htest(t.test(1:10, y = c(7:20, 200)))
  expect_equivalent(testy1,testy2)

  # wilcox -- same data
  expect_equal(wilcox.test(1:10, y = c(7:20))$pvalue,
               simple_htest(1:10, y = c(7:20), test = "w")$pvalue)
  expect_equal(wilcox.test(1:10, y = c(7:20, 200))$pvalue,
               simple_htest(1:10, y = c(7:20, 200), test = "w")$pvalue)

  testy1 = describe_htest(simple_htest(1:10, y = c(7:20, 200), test = "w"))
  testy2 = describe_htest(wilcox.test(1:10, y = c(7:20, 200),
                                      conf.int = TRUE))
  expect_equivalent(testy1,testy2)

  testy1 = df_htest(simple_htest(1:10, y = c(7:20, 200), test = "w"))
  testy2 = df_htest(wilcox.test(1:10, y = c(7:20, 200),
                                      conf.int = TRUE))
  expect_equivalent(testy1,testy2)


  expect_error(simple_htest(1:10, y = c(7:20, 200), test = "w",
                            alternative = "e"))

  # Test equivalence -----
  testy1 = simple_htest(data = mtcars,
               mpg ~ am, test = "w",
               mu = 3,
               alternative = "e")
  testy2 = wilcox.test(data = mtcars,
                        mpg ~ am,
                        mu = -3,
                        alternative = "g")
  expect_equal(testy1$p.value,testy2$p.value)
  testy1 = simple_htest(data = mtcars,
                        mpg ~ am, test = "t",
                        mu = 3,
                        alternative = "e")
  testy2 = t.test(data = mtcars,
                       mpg ~ am,
                       mu = -3,
                       alternative = "g")
  expect_equal(testy1$p.value,testy2$p.value)


  # Test MET -----
  testy1 = simple_htest(data = mtcars,
                        mpg ~ am, test = "w",
                        mu = 3,
                        alternative = "m")
  testy2 = wilcox.test(data = mtcars,
                       mpg ~ am,
                       mu = -3,
                       alternative = "l")
  testydf = df_htest(testy1)
  expect_equal(testy1$p.value,testy2$p.value)
  expect_equal(testy1$p.value, testydf$p.value)
  testy1 = simple_htest(data = mtcars,
                        mpg ~ am, test = "t",
                        mu = 3,
                        alternative = "m")
  testy2 = t.test(data = mtcars,
                  mpg ~ am,
                  mu = -3,
                  alternative = "l")
  expect_equal(testy1$p.value,testy2$p.value)

  # Errors & Messages

  expect_error(df_htest(1),"htest must be of the class htest")
  expect_error(describe_htest(1),"htest must be of the class htest")
  testy2 = t.test(data = mtcars,
                  mpg ~ am,
                  mu = -3,
                  alternative = "l")
  testy2$p.value = NULL
  expect_error(describe_htest(testy2))
  ow = stats::oneway.test(extra ~ group, data = sleep)
  expect_message(describe_htest(ow))
  wil_noci =  wilcox.test(data = mtcars,
                                    mpg ~ am,
                                    mu = -3,
                                    alternative = "l")
  expect_message(describe_htest(wil_noci))

  expect_equal(df_htest(wil_noci, extract_names = FALSE)$statistic,
              73)

  expect_equal(df_htest(testy2, show_ci = FALSE)$conf.level,
               NULL)
  expect_equal(df_htest(testy2, test_statistics = FALSE)$t,
               NULL)
})

test_that("All other htests",{
  # Quade ---

  ## Conover (1999, p. 375f):
  ## Numbers of five brands of a new hand lotion sold in seven stores
  ## during one week.
  y <- matrix(c( 5,  4,  7, 10, 12,
                 1,  3,  1,  0,  2,
                 16, 12, 22, 22, 35,
                 5,  4,  3,  5,  4,
                 10,  9,  7, 13, 10,
                 19, 18, 28, 37, 58,
                 10,  7,  6,  8,  7),
              nrow = 7, byrow = TRUE,
              dimnames =
                list(Store = as.character(1:7),
                     Brand = LETTERS[1:5]))

  htest <- stats::quade.test(y)
  df1 = df_htest(htest = htest)
  pr = describe_htest(htest, alpha = .05)
  expect_equal(htest$p.value,df1$p.value)

  # McNemar ----

  ## Agresti (1990), p. 350.
  ## Presidential Approval Ratings.
  ##  Approval of the President's performance in office in two surveys,
  ##  one month apart, for a random sample of 1600 voting-age Americans.
  Performance <-
    matrix(c(794, 86, 150, 570),
           nrow = 2,
           dimnames = list("1st Survey" = c("Approve", "Disapprove"),
                           "2nd Survey" = c("Approve", "Disapprove")))

  htest = stats::mcnemar.test(Performance)
  df1 = df_htest(htest)
  pr = describe_htest(htest, alpha = .05)
  expect_equal(htest$p.value,df1$p.value)

  # Friedman ---
  ## Hollander & Wolfe (1973), p. 140ff.
  ## Comparison of three methods ("round out", "narrow angle", and
  ##  "wide angle") for rounding first base.  For each of 18 players
  ##  and the three method, the average time of two runs from a point on
  ##  the first base line 35ft from home plate to a point 15ft short of
  ##  second base is recorded.
  RoundingTimes <-
    matrix(c(5.40, 5.50, 5.55,
             5.85, 5.70, 5.75,
             5.20, 5.60, 5.50,
             5.55, 5.50, 5.40,
             5.90, 5.85, 5.70,
             5.45, 5.55, 5.60,
             5.40, 5.40, 5.35,
             5.45, 5.50, 5.35,
             5.25, 5.15, 5.00,
             5.85, 5.80, 5.70,
             5.25, 5.20, 5.10,
             5.65, 5.55, 5.45,
             5.60, 5.35, 5.45,
             5.05, 5.00, 4.95,
             5.50, 5.50, 5.40,
             5.45, 5.55, 5.50,
             5.55, 5.55, 5.35,
             5.45, 5.50, 5.55,
             5.50, 5.45, 5.25,
             5.65, 5.60, 5.40,
             5.70, 5.65, 5.55,
             6.30, 6.30, 6.25),
           nrow = 22,
           byrow = TRUE,
           dimnames = list(1 : 22,
                           c("Round Out", "Narrow Angle", "Wide Angle")))
  htest = friedman.test(RoundingTimes)
  df1 = df_htest(htest)
  pr = describe_htest(htest, alpha = .05)
  expect_equal(htest$p.value,df1$p.value)

  # Kruskal-Wallis -----

  ## Hollander & Wolfe (1973), 116.
  ## Mucociliary efficiency from the rate of removal of dust in normal
  ##  subjects, subjects with obstructive airway disease, and subjects
  ##  with asbestosis.
  x <- c(2.9, 3.0, 2.5, 2.6, 3.2) # normal subjects
  y <- c(3.8, 2.7, 4.0, 2.4)      # with obstructive airway disease
  z <- c(2.8, 3.4, 3.7, 2.2, 2.0) # with asbestosis
  htest = kruskal.test(list(x, y, z))
  df1 = df_htest(htest)
  pr = describe_htest(htest, alpha = .05)
  expect_equal(htest$p.value,df1$p.value)

  # Chi-squared -----

  ## Effect of simulating p-values
  x <- matrix(c(12, 5, 7, 7), ncol = 2)
  #chisq.test(x)$p.value           # 0.4233
  htest = chisq.test(x, simulate.p.value = TRUE, B = 10000)
  df1 = df_htest(htest)
  pr = describe_htest(htest, alpha = .05)
  expect_equal(htest$p.value,df1$p.value)

  # Fisher exact ----
  ## Agresti (1990, p. 61f; 2002, p. 91) Fisher's Tea Drinker
  ## A British woman claimed to be able to distinguish whether milk or
  ##  tea was added to the cup first.  To test, she was given 8 cups of
  ##  tea, in four of which milk was added first.  The null hypothesis
  ##  is that there is no association between the true order of pouring
  ##  and the woman's guess, the alternative that there is a positive
  ##  association (that the odds ratio is greater than 1).
  TeaTasting <-
    matrix(c(3, 1, 1, 3),
           nrow = 2,
           dimnames = list(Guess = c("Milk", "Tea"),
                           Truth = c("Milk", "Tea")))
  htest = fisher.test(TeaTasting, alternative = "greater")
  df1 = df_htest(htest)
  pr = describe_htest(htest, alpha = .05)
  expect_equal(htest$p.value,df1$p.value)

  # proportions test -----

  ## Data from Fleiss (1981), p. 139.
  ## H0: The null hypothesis is that the four populations from which
  ##     the patients were drawn have the same true proportion of smokers.
  ## A:  The alternative is that this proportion is different in at
  ##     least one of the populations.

  smokers  <- c( 83, 90, 129, 70 )
  patients <- c( 86, 93, 136, 82 )
  htest = prop.test(smokers, patients)
  df1 = df_htest(htest)
  pr = describe_htest(htest, alpha = .05)
  expect_equal(htest$p.value,df1$p.value)

  # binomial test -----

  ## Conover (1971), p. 97f.
  ## Under (the assumption of) simple Mendelian inheritance, a cross
  ##  between plants of two particular genotypes produces progeny 1/4 of
  ##  which are "dwarf" and 3/4 of which are "giant", respectively.
  ##  In an experiment to determine if this assumption is reasonable, a
  ##  cross results in progeny having 243 dwarf and 682 giant plants.
  ##  If "giant" is taken as success, the null hypothesis is that p =
  ##  3/4 and the alternative that p != 3/4.
  htest = binom.test(c(682, 243), p = 3/4)
  df1 = df_htest(htest)
  pr = describe_htest(htest, alpha = .05)
  expect_equal(htest$p.value,df1$p.value)

})
