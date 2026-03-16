# need hush function to run print through examples

hush = function(code) {
  sink(nullfile())
  tmp = code
  sink()
  return(tmp)
}


test_that("Run examples for two sample", {
  skip_on_cran()

  set.seed(2525894)
  dat1 = mtcars
  expect_error(log_TOST(
    x = dat1$mpg,
    var.equal = TRUE,
    eqb = 1.25
  ))
  test1 = boot_log_TOST(data = dat1,
                        mpg ~ am,
                        var.equal = TRUE,
                        eqb = 1.25,
                        R = 199)

  expect_message({
    boot_log_TOST(data = dat1,
                  mpg ~ am,
                  var.equal = TRUE,
                  hypothesis = "MET",
                  eqb = 1.01,
                  R = 99)
  })

  expect_message({
   log_TOST(data = dat1,
                  mpg ~ am,
                  var.equal = TRUE,
                  hypothesis = "MET",
                  eqb = 1.01)
  })

  test1_t = log_TOST(data = dat1,
                        mpg ~ am,
                        var.equal = TRUE,
                        eqb = 1.25)

  expect_equal(test1$effsize$estimate,
               test1_t$effsize$estimate)

  test1 = boot_log_TOST(data = dat1,
                        mpg ~ am,
                        var.equal = TRUE,
                        eqb = 1.75,
                        R = 199)

  test1_t = log_TOST(data = dat1,
                     mpg ~ am,
                     var.equal = TRUE,
                     eqb = 1.75)

  expect_equal(test1$effsize$estimate,
               test1_t$effsize$estimate)

  test1 = boot_log_TOST(data = dat1,
                        mpg ~ am,
                        var.equal = TRUE,
                        hypothesis = "MET",
                        eqb = 1.25,
                        R = 199)

  test1_t = log_TOST(data = dat1,
                     mpg ~ am,
                     var.equal = TRUE,
                     hypothesis = "MET",
                     eqb = 1.25)

  expect_equal(test1$effsize$estimate,
               test1_t$effsize$estimate)

  test2 = boot_log_TOST(data = dat1,
                        mpg ~ am,
                        var.equal = FALSE,
                        eqb = 1.25,
                        R = 199)

  test2_t = log_TOST(data = dat1,
                     mpg ~ am,
                     var.equal = FALSE,
                     eqb = 1.25)

  expect_equal(test2$effsize$estimate,
               test2_t$effsize$estimate)

  test2 = boot_log_TOST(data = dat1,
                        mpg ~ am,
                        var.equal = FALSE,
                        eqb = 1.25,
                        R = 199,
                        hypothesis = "MET")

  test2_t = log_TOST(data = dat1,
                     mpg ~ am,
                     var.equal = FALSE,
                     eqb = 1.25,
                     hypothesis = "MET")

  expect_equal(test2$effsize$estimate,
               test2_t$effsize$estimate)

})


test_that("Run examples for paired samples", {
  skip_on_cran()

  set.seed(922287)
  sleep2 = sleep
  sleep2$sleep = sleep2$extra + 4

  # Test error for negative values in extra column (using vectors for paired test)
  expect_error(boot_log_TOST(
    x = sleep2$extra[sleep2$group == 1],
    y = sleep2$extra[sleep2$group == 2],
    paired = TRUE,
    eqb = 1.25,
    R = 199
  ))
  expect_error(log_TOST(
    x = sleep2$extra[sleep2$group == 1],
    y = sleep2$extra[sleep2$group == 2],
    paired = TRUE,
    eqb = 1.25
  ))

  # Test error for null = -1
  expect_error(boot_log_TOST(
    x = sleep2$sleep[sleep2$group == 1],
    y = sleep2$sleep[sleep2$group == 2],
    paired = TRUE,
    eqb = 1.25,
    null = -1,
    R = 199
  ))
  expect_error(log_TOST(
    x = sleep2$sleep[sleep2$group == 1],
    y = sleep2$sleep[sleep2$group == 2],
    paired = TRUE,
    eqb = 1.25,
    null = -1
  ))

  # Test error for alpha = -1
  expect_error(boot_log_TOST(
    x = sleep2$sleep[sleep2$group == 1],
    y = sleep2$sleep[sleep2$group == 2],
    paired = TRUE,
    eqb = 1.25,
    alpha = -1,
    R = 199
  ))
  expect_error(log_TOST(
    x = sleep2$sleep[sleep2$group == 1],
    y = sleep2$sleep[sleep2$group == 2],
    paired = TRUE,
    eqb = 1.25,
    alpha = -1
  ))

  # Valid tests using vectors for paired samples
  test1 = boot_log_TOST(
    x = sleep2$sleep[sleep2$group == 1],
    y = sleep2$sleep[sleep2$group == 2],
    paired = TRUE,
    eqb = 1.25,
    R = 199
  )

  test1_t = log_TOST(
    x = sleep2$sleep[sleep2$group == 1],
    y = sleep2$sleep[sleep2$group == 2],
    paired = TRUE,
    eqb = 1.25
  )

  expect_equal(test1$effsize$estimate,
               test1_t$effsize$estimate)

  test1 = boot_log_TOST(
    x = sleep2$sleep[sleep2$group == 1],
    y = sleep2$sleep[sleep2$group == 2],
    paired = TRUE,
    hypothesis = "MET",
    eqb = 1.25,
    R = 199
  )

  test1_t = log_TOST(
    x = sleep2$sleep[sleep2$group == 1],
    y = sleep2$sleep[sleep2$group == 2],
    paired = TRUE,
    hypothesis = "MET",
    eqb = 1.25
  )

  expect_equal(test1$effsize$estimate,
               test1_t$effsize$estimate)

  test1 = boot_log_TOST(
    x = sleep2$sleep[sleep2$group == 1],
    y = sleep2$sleep[sleep2$group == 2],
    paired = TRUE,
    hypothesis = "MET",
    eqb = 1.05,
    R = 199
  )

  test1_t = log_TOST(
    x = sleep2$sleep[sleep2$group == 1],
    y = sleep2$sleep[sleep2$group == 2],
    paired = TRUE,
    hypothesis = "MET",
    eqb = 1.05
  )

  expect_equal(test1$effsize$estimate,
               test1_t$effsize$estimate)


  test1 = boot_log_TOST(
    x = sleep2$sleep[sleep2$group == 1],
    y = sleep2$sleep[sleep2$group == 2],
    paired = TRUE,
    eqb = 1.9,
    R = 199
  )

  test1_t = log_TOST(
    x = sleep2$sleep[sleep2$group == 1],
    y = sleep2$sleep[sleep2$group == 2],
    paired = TRUE,
    eqb = 1.9
  )

  expect_equal(test1$effsize$estimate,
               test1_t$effsize$estimate)
  
  # Test that paired = TRUE shows a message with formula method
  expect_message(
    log_TOST(sleep ~ group, data = sleep2, paired = TRUE, eqb = 1.25),
    "Using 'paired = TRUE' with the formula interface is not recommended"
  )
  
  expect_message(
    boot_log_TOST(sleep ~ group, data = sleep2, paired = TRUE, eqb = 1.25, R = 199),
    "Using 'paired = TRUE' with the formula interface is not recommended"
  )

})

# boot_log_TOST CI/p-value agreement tests -----

for (ci_method in c("perc", "basic", "bca", "stud")) {
  test_that(paste0("boot_log_TOST: CI/p-value agreement for ", ci_method), {
    skip_on_cran()

    set.seed(42)
    x <- rlnorm(30, meanlog = 3.5, sdlog = 0.4)
    y <- rlnorm(30, meanlog = 3.6, sdlog = 0.4)

    hush = function(code) {
      sink(nullfile())
      tmp = code
      sink()
      return(tmp)
    }

    res <- hush(boot_log_TOST(x = x, y = y, eqb = 1.25,
                               boot_ci = ci_method, R = 1999))
    # Check two-sided on log scale with 90% CI (TOST uses 1-2*alpha = 90%)
    # 90% CI excludes null iff two-sided p < 2*alpha = 0.10
    ci <- c(res$effsize$lower.ci[1], res$effsize$upper.ci[1])
    log_null <- log(1)  # = 0
    ci_excludes_null <- ci[1] > log_null || ci[2] < log_null
    p_rejects <- res$TOST$p.value[1] < 0.10
    expect_equal(ci_excludes_null, p_rejects,
                 label = paste(ci_method, "CI/p agreement two.sided"))
  })
}
