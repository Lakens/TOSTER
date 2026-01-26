# All htest extensions tested here

test_that("simple_htest: t-test & wilcox", {

  set.seed(3164964)

  expect_equal(t.test(1:10, y = c(7:20))$p.value,
               simple_htest(1:10, y = c(7:20))$p.value)
  expect_equal(t.test(1:10, y = c(7:20, 200))$p.value,
               simple_htest(1:10, y = c(7:20, 200))$p.value)

  testy1 = describe_htest(simple_htest(1:10, y = c(7:20, 200)))
  testy2 = describe_htest(t.test(1:10, y = c(7:20, 200)))
  expect_equal(testy1,testy2)

  testy1 = df_htest(simple_htest(1:10, y = c(7:20, 200)))
  testy2 = df_htest(t.test(1:10, y = c(7:20, 200)))
  expect_equal(testy1,testy2)

  testy1 = describe_htest(as_htest(t_TOST(1:10, y = c(7:20, 200), eqb = 3)))
  testy1 = describe_htest(as_htest(t_TOST(1:10, y = c(7:20, 200), eqb = 3,
                                          hypothesis = "MET")))

  # wilcox -- same data
  expect_equal(wilcox.test(1:10, y = c(7:20))$p.value,
               simple_htest(1:10, y = c(7:20), test = "w")$p.value)
  expect_equal(wilcox.test(1:10, y = c(7:20, 200))$p.value,
               simple_htest(1:10, y = c(7:20, 200), test = "w")$p.value)

  testy1 = describe_htest(simple_htest(1:10, y = c(7:20, 200), test = "w"))
  testy2 = describe_htest(wilcox.test(1:10, y = c(7:20, 200),
                                      conf.int = TRUE))
  expect_equal(testy1,testy2)

  testy1 = df_htest(simple_htest(1:10, y = c(7:20, 200), test = "w"))
  testy2 = df_htest(wilcox.test(1:10, y = c(7:20, 200),
                                      conf.int = TRUE))
  expect_equal(testy1,testy2)

  testy1 = describe_htest(as_htest(wilcox_TOST(1:10, y = c(7:20, 200), eqb = 3)))
  testy1 = describe_htest(as_htest(wilcox_TOST(1:10, y = c(7:20, 200), eqb = 3,
                                               hypothesis = "MET")))

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

test_that("brunner_munzel",{

  set.seed(2808)

  expect_equal(brunner_munzel(1:10, y = c(7:20))$p.value,
               simple_htest(1:10, y = c(7:20),
                            test = "b",
                            mu = .5)$p.value)
  expect_equal(brunner_munzel(1:10, y = c(7:20, 200))$p.value,
               simple_htest(1:10, y = c(7:20, 200),
                            test = "b",
                            mu = .5)$p.value)

  testy1 = describe_htest(simple_htest(1:10, y = c(7:20, 200),
                                       test = "b",
                                       mu = .5))
  testy2 = describe_htest(brunner_munzel(1:10, y = c(7:20, 200)))
  expect_equal(testy1,testy2)

  testy1 = df_htest(simple_htest(1:10, y = c(7:20, 200),
                                 test = "b",
                                 mu = .5))
  testy2 = df_htest(brunner_munzel(1:10, y = c(7:20, 200)))
  expect_equal(testy1,testy2)

  testy1 = simple_htest(data = mtcars,
                        mpg ~ am,
                        test = "b",
                        mu = .75,
                        alternative = "e")
  testy2 = brunner_munzel(data = mtcars,
                       mpg ~ am,
                       mu = .25,
                       alternative = "g")
  expect_equal(testy1$p.value,testy2$p.value)
  set.seed(1944)
  testy1 = simple_htest(data = mtcars,
                        mpg ~ am,
                        test = "b",
                        mu = .25,
                        alternative = "g",
                        perm = TRUE)
  set.seed(1944)
  testy2 = brunner_munzel(data = mtcars,
                          mpg ~ am,
                          mu = .25,
                          alternative = "g",
                          perm = TRUE)
  expect_equal(testy1$p.value,testy2$p.value)

  set.seed(1945)
  testy1 = simple_htest(data = mtcars,
                        mpg ~ am,
                        test = "b",
                        mu = .25,
                        alternative = "l",
                        perm = TRUE)
  set.seed(1945)
  testy2 = brunner_munzel(data = mtcars,
                          mpg ~ am,
                          mu = .25,
                          alternative = "l",
                          perm = TRUE)
  expect_equal(testy1$p.value,testy2$p.value)

  set.seed(1946)
  testy1 = simple_htest(data = mtcars,
                        mpg ~ am,
                        test = "b",
                        mu = .5,
                        alternative = "t",
                        perm = TRUE)
  set.seed(1946)
  testy2 = brunner_munzel(data = mtcars,
                          mpg ~ am,
                          mu = .5,
                          alternative = "t",
                          perm = TRUE)
  expect_equal(testy1$p.value,testy2$p.value)

  # Test equivalence -----
  testy1 = simple_htest(data = mtcars,
                        mpg ~ am, test = "b",
                        mu = .4,
                        alternative = "e")
  testy2 = brunner_munzel(data = mtcars,
                       mpg ~ am,
                       mu = .4,
                       alternative = "g")
  expect_equal(testy1$p.value,testy2$p.value)
  testy1 = simple_htest(data = mtcars,
                        mpg ~ am, test = "b",
                        mu = .3,
                        alternative = "e")
  testy2 = brunner_munzel(data = mtcars,
                  mpg ~ am,
                  mu = .3,
                  alternative = "g")
  expect_equal(testy1$p.value,testy2$p.value)


  # Test MET -----
  testy1 = simple_htest(data = mtcars,
                        mpg ~ am, test = "b",
                        mu = .3,
                        alternative = "m")
  testy2 = brunner_munzel(data = mtcars,
                       mpg ~ am,
                       mu = .3,
                       alternative = "l")
  testydf = df_htest(testy1)
  expect_equal(testy1$p.value,testy2$p.value)
  expect_equal(testy1$p.value, testydf$p.value)

  # Paired ------
  data(sleep)

  set.seed(1944)
  # Use vector form for paired tests, not formula
  testy1 = simple_htest(x = sleep$extra[sleep$group == 1],
                        y = sleep$extra[sleep$group == 2],
                        test = "b",
                        mu = .25,
                        alternative = "g",
                        perm = TRUE,
                        paired = TRUE)
  set.seed(1944)
  testy2 = brunner_munzel(x = sleep$extra[sleep$group == 1],
                          y = sleep$extra[sleep$group == 2],
                          mu = .25,
                          alternative = "g",
                          perm = TRUE,
                          paired = TRUE)
  expect_equal(testy1$p.value,testy2$p.value)

  test_big = brunner_munzel(x = rnorm(100),
                            y = rnorm(100),
                          mu = .25,
                          alternative = "t",
                          perm = TRUE,
                          paired = TRUE)

  set.seed(1945)
  # Use vector form for paired tests, not formula
  testy1 = simple_htest(x = sleep$extra[sleep$group == 1],
                        y = sleep$extra[sleep$group == 2],
                        test = "b",
                        mu = .25,
                        alternative = "l",
                        perm = TRUE,
                        paired = TRUE)
  set.seed(1945)
  testy2 = brunner_munzel(x = sleep$extra[sleep$group == 1],
                          y = sleep$extra[sleep$group == 2],
                          mu = .25,
                          alternative = "l",
                          perm = TRUE,
                          paired = TRUE)
  expect_equal(testy1$p.value,testy2$p.value)

  set.seed(1946)
  # Use vector form for paired tests, not formula
  testy1 = simple_htest(x = sleep$extra[sleep$group == 1],
                        y = sleep$extra[sleep$group == 2],
                        test = "b",
                        mu = .5,
                        alternative = "t",
                        perm = TRUE,
                        paired = TRUE)
  set.seed(1946)
  testy2 = brunner_munzel(x = sleep$extra[sleep$group == 1],
                          y = sleep$extra[sleep$group == 2],
                          mu = .5,
                          alternative = "t",
                          perm = TRUE,
                          paired = TRUE)
  expect_equal(testy1$p.value,testy2$p.value)

  ## Equ ---
  # Use vector form for paired tests, not formula
  testy1 = simple_htest(x = sleep$extra[sleep$group == 1],
                        y = sleep$extra[sleep$group == 2],
                        paired = TRUE,
                        test = "b",
                        mu = .3,
                        alternative = "e")
  testy2 = brunner_munzel(x = sleep$extra[sleep$group == 1],
                          y = sleep$extra[sleep$group == 2],
                          paired = TRUE,
                          mu = .3,
                          alternative = "g")
  expect_equal(testy1$p.value,testy2$p.value)


  ## MET -----
  testy1 = simple_htest(x = sleep$extra[sleep$group == 1],
                        y = sleep$extra[sleep$group == 2],
                        paired = TRUE,
                        test = "b",
                        mu = .3,
                        alternative = "m")
  testy2 = brunner_munzel(x = sleep$extra[sleep$group == 1],
                          y = sleep$extra[sleep$group == 2],
                          paired = TRUE,
                          mu = .3,
                          alternative = "l")
  testydf = df_htest(testy1)
  expect_equal(testy1$p.value,testy2$p.value)
  expect_equal(testy1$p.value, testydf$p.value)

  # Direct equivalence/minimal.effect tests with mu as length-2 vector ----

  # Two-sample: equivalence with t-method
  eq_test <- brunner_munzel(data = mtcars,
                            mpg ~ am,
                            alternative = "equivalence",
                            mu = c(0.3, 0.7))
  expect_true(eq_test$alternative == "equivalence")
  expect_equal(length(eq_test$null.value), 2)
  expect_true(eq_test$p.value >= 0 && eq_test$p.value <= 1)
  # CI should be at 1-2*alpha level

  expect_equal(attr(eq_test$conf.int, "conf.level"), 0.9)

  # Two-sample: minimal.effect with t-method
  met_test <- brunner_munzel(data = mtcars,
                             mpg ~ am,
                             alternative = "minimal.effect",
                             mu = c(0.4, 0.6))
  expect_true(met_test$alternative == "minimal.effect")
  expect_equal(length(met_test$null.value), 2)
  expect_true(met_test$p.value >= 0 && met_test$p.value <= 1)

  # Two-sample: equivalence with logit method
  eq_logit <- brunner_munzel(data = mtcars,
                             mpg ~ am,
                             alternative = "equivalence",
                             mu = c(0.3, 0.7),
                             test_method = "logit")
  expect_true(eq_logit$alternative == "equivalence")
  # CI should stay within [0, 1]
  expect_true(eq_logit$conf.int[1] >= 0 && eq_logit$conf.int[2] <= 1)

  # Two-sample: equivalence with permutation method
  set.seed(12345)
  eq_perm <- brunner_munzel(data = mtcars,
                            mpg ~ am,
                            alternative = "equivalence",
                            mu = c(0.3, 0.7),
                            test_method = "perm",
                            R = 999)
  expect_true(eq_perm$alternative == "equivalence")
  expect_true(eq_perm$p.value >= 0 && eq_perm$p.value <= 1)

  # Two-sample: minimal.effect with permutation method
  set.seed(12346)
  met_perm <- brunner_munzel(data = mtcars,
                             mpg ~ am,
                             alternative = "minimal.effect",
                             mu = c(0.4, 0.6),
                             test_method = "perm",
                             R = 999)
  expect_true(met_perm$alternative == "minimal.effect")
  expect_true(met_perm$p.value >= 0 && met_perm$p.value <= 1)

  # Paired: equivalence with t-method
  eq_paired <- brunner_munzel(x = sleep$extra[sleep$group == 1],
                              y = sleep$extra[sleep$group == 2],
                              paired = TRUE,
                              alternative = "equivalence",
                              mu = c(0.3, 0.7))
  expect_true(eq_paired$alternative == "equivalence")
  expect_equal(length(eq_paired$null.value), 2)

  # Paired: minimal.effect with t-method
  met_paired <- brunner_munzel(x = sleep$extra[sleep$group == 1],
                               y = sleep$extra[sleep$group == 2],
                               paired = TRUE,
                               alternative = "minimal.effect",
                               mu = c(0.4, 0.6))
  expect_true(met_paired$alternative == "minimal.effect")

  # Paired: equivalence with logit method
  eq_paired_logit <- brunner_munzel(x = sleep$extra[sleep$group == 1],
                                    y = sleep$extra[sleep$group == 2],
                                    paired = TRUE,
                                    alternative = "equivalence",
                                    mu = c(0.3, 0.7),
                                    test_method = "logit")
  expect_true(eq_paired_logit$alternative == "equivalence")

  # Paired: equivalence with permutation method
  set.seed(12347)
  eq_paired_perm <- brunner_munzel(x = sleep$extra[sleep$group == 1],
                                   y = sleep$extra[sleep$group == 2],
                                   paired = TRUE,
                                   alternative = "equivalence",
                                   mu = c(0.3, 0.7),
                                   test_method = "perm",
                                   R = 999)
  expect_true(eq_paired_perm$alternative == "equivalence")

  # Consistency check: equivalence p-value should equal max of one-sided tests
  lo_test <- brunner_munzel(data = mtcars,
                            mpg ~ am,
                            mu = 0.3,
                            alternative = "greater")
  hi_test <- brunner_munzel(data = mtcars,
                            mpg ~ am,
                            mu = 0.7,
                            alternative = "less")
  eq_direct <- brunner_munzel(data = mtcars,
                              mpg ~ am,
                              alternative = "equivalence",
                              mu = c(0.3, 0.7))
  expect_equal(eq_direct$p.value, max(lo_test$p.value, hi_test$p.value),
               tolerance = 1e-10)

  # Consistency check: minimal.effect p-value should equal min of one-sided tests
  lo_test_met <- brunner_munzel(data = mtcars,
                                mpg ~ am,
                                mu = 0.4,
                                alternative = "less")
  hi_test_met <- brunner_munzel(data = mtcars,
                                mpg ~ am,
                                mu = 0.6,
                                alternative = "greater")
  met_direct <- brunner_munzel(data = mtcars,
                               mpg ~ am,
                               alternative = "minimal.effect",
                               mu = c(0.4, 0.6))
  expect_equal(met_direct$p.value, min(lo_test_met$p.value, hi_test_met$p.value),
               tolerance = 1e-10)

  # Error handling for equivalence/minimal.effect ----

  # mu must be length 2 for equivalence
  expect_error(brunner_munzel(data = mtcars, mpg ~ am,
                              alternative = "equivalence",
                              mu = 0.5))

  # mu[1] must be less than mu[2]
  expect_error(brunner_munzel(data = mtcars, mpg ~ am,
                              alternative = "equivalence",
                              mu = c(0.7, 0.3)))

  # mu values must be between 0 and 1
  expect_error(brunner_munzel(data = mtcars, mpg ~ am,
                              alternative = "equivalence",
                              mu = c(-0.1, 0.7)))

  expect_error(brunner_munzel(data = mtcars, mpg ~ am,
                              alternative = "equivalence",
                              mu = c(0.3, 1.1)))

  # Errors ----

  expect_error(brunner_munzel(x = rnorm(100),
                              #y = rnorm(100),
                              mu = .25,
                              alternative = "t",
                              perm = TRUE,
                              paired = TRUE))

  expect_error(brunner_munzel(x = "error",
                              #y = rnorm(100),
                              mu = .25,
                              alternative = "t",
                              perm = TRUE,
                              paired = TRUE))

  expect_error(brunner_munzel(x = rnorm(100),
                              y = rnorm(99),
                              mu = .25,
                              alternative = "t",
                              perm = TRUE,
                              paired = TRUE))

  expect_message(brunner_munzel(x = rnorm(300),
                              y = rnorm(300),
                              mu = .25,
                              alternative = "t",
                              perm = TRUE,
                              #paired = TRUE
                              ))



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

test_that("plot_htest_est works correctly", {

  # Standard two-sample t-test (should auto-convert to difference)
  t_two <- t.test(extra ~ group, data = sleep)
  p1 <- plot_htest_est(t_two)
  expect_s3_class(p1, "ggplot")

  # One-sample t-test
  t_one <- t.test(sleep$extra, mu = 0)
  p2 <- plot_htest_est(t_one)
  expect_s3_class(p2, "ggplot")

  # Correlation test
  cor_res <- cor.test(mtcars$mpg, mtcars$wt)
  p3 <- plot_htest_est(cor_res)
  expect_s3_class(p3, "ggplot")

  # TOST converted to htest (equivalence bounds)
  tost_res <- t_TOST(extra ~ group, data = sleep, eqb = 1)
  htest_tost <- as_htest(tost_res)
  p4 <- plot_htest_est(htest_tost)
  expect_s3_class(p4, "ggplot")

  # Error cases
  expect_error(plot_htest_est("not_htest"),
               "Input must be an object of class")

  # htest without estimate
  htest_no_est <- t_two
  htest_no_est$estimate <- NULL
  expect_error(plot_htest_est(htest_no_est),
               "htest object has no estimate")

  # htest without conf.int
  htest_no_ci <- t_two
  htest_no_ci$conf.int <- NULL
  expect_error(plot_htest_est(htest_no_ci),
               "htest object has no confidence interval")

  # Wilcoxon test with CI
  wilcox_res <- wilcox.test(extra ~ group, data = sleep, conf.int = TRUE)
  p5 <- plot_htest_est(wilcox_res)
  expect_s3_class(p5, "ggplot")

  # Test describe argument
  p_desc <- plot_htest_est(t_two, describe = TRUE)
  expect_s3_class(p_desc, "ggplot")
  expect_true(!is.null(p_desc$labels$subtitle))
  expect_true(grepl("t\\(", p_desc$labels$subtitle))  # Should contain t statistic

  p_no_desc <- plot_htest_est(t_two, describe = FALSE)
  expect_s3_class(p_no_desc, "ggplot")
  expect_true(is.null(p_no_desc$labels$subtitle))

  # Test with equivalence test (two null values)
  res_eq <- simple_htest(x = sleep$extra[sleep$group == 1],
                         y = sleep$extra[sleep$group == 2],
                         paired = TRUE, mu = 1, alternative = "e")
  p_eq <- plot_htest_est(res_eq)
  expect_s3_class(p_eq, "ggplot")
  expect_true(!is.null(p_eq$labels$subtitle))

})
