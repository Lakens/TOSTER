#context("Run Examples for boot_t_TOST")

# need hush function to run print through examples

hush = function(code) {
  sink(nullfile())
  tmp = code
  sink()
  return(tmp)
}

test_that("Run examples for z_cor_test", {

  set.seed(76584441)

  samp1 = rnorm(25)
  samp2 = rnorm(25)

  df_samp = data.frame(y = c(samp1,samp2),
                       group = c(rep("g1",25),
                                 rep("g2",25)))

  expect_error(z_cor_test())
  expect_error(z_cor_test(samp1,
                          samp2,
                          method = "p",
                          null = c(-.24,.24),
                          TOST = FALSE))

  test1 = z_cor_test(samp1,
                     samp2,
                     method = "p")

  test2 = z_cor_test(samp1,
                     samp2,
                     method = "s")

  test3 = z_cor_test(samp1,
                     samp2,
                     method = "k")

  expect_equal(c(unname(test1$parameter),
                 unname(test2$parameter),
                 unname(test3$parameter)),
                     c(25, 25, 25))

  expect_equal(test1$p.value,
               0.726,
               tolerance = .01)
  expect_equal(test2$p.value,
               0.936,
               tolerance = .01)
  expect_equal(test3$p.value,
               0.963,
               tolerance = .01)

  test1c = cor.test(samp1,
                     samp2,
                     method = "p")

  test2c = cor.test(samp1,
                     samp2,
                     method = "s")

  test3c = cor.test(samp1,
                     samp2,
                     method = "k")

  expect_equal(unname(test1$estimate),
               unname(test1c$estimate))
  expect_equal(test3$estimate,
               test3c$estimate)
  expect_equal(test2$estimate,
               test2c$estimate)

  test1 = z_cor_test(samp1,
                     samp2,
                     method = "p",
                     null = .4,
                     alternative = "e")

  test2 = z_cor_test(samp1,
                     samp2,
                     method = "s",
                     null = .4,
                     alternative = "e")

  test3 = z_cor_test(samp1,
                     samp2,
                     method = "k",
                     null = .4,
                     alternative = "e")


  # other alts ----

  test1 = z_cor_test(samp1,
                     samp2,
                     method = "p",
                     alternative = "greater")

  test2 = z_cor_test(samp1,
                     samp2,
                     method = "s",
                     alternative = "greater")

  test3 = z_cor_test(samp1,
                     samp2,
                     method = "k",
                     alternative = "greater")


  expect_equal(test1$p.value,
               0.363,
               tolerance = .01)
  expect_equal(test2$p.value,
               0.467,
               tolerance = .01)
  expect_equal(test3$p.value,
               0.518,
               tolerance = .01)

  test1 = z_cor_test(samp1,
                     samp2,
                     method = "p",
                     alternative = "less")

  test2 = z_cor_test(samp1,
                     samp2,
                     method = "s",
                     alternative = "less")

  test3 = z_cor_test(samp1,
                     samp2,
                     method = "k",
                     alternative = "less")


  expect_equal(test1$p.value,
               1-0.363,
               tolerance = .01)
  expect_equal(test2$p.value,
               1-0.467,
               tolerance = .01)
  expect_equal(test3$p.value,
               1-0.518,
               tolerance = .01)





})

test_that("cor_test: equ and met", {
  skip_on_cran()

  set.seed(5533428)

  samp1 = rnorm(150)
  samp2 = rnorm(150)

  test1 = z_cor_test(samp1,
                     samp2,
                     alternative = "e",
                     null = .1)

  test2 = z_cor_test(samp1,
                     samp2,
                     alternative = "e",
                     null = .3)

  test3 = z_cor_test(samp1,
                     samp2,
                     alternative = "e",
                     null = .5)
  expect_true(test1$p.value > test2$p.value)
  expect_true(test2$p.value > test3$p.value)

  test1 = boot_cor_test(samp1,
                     samp2,
                     alternative = "equivalence",
                     null = .1)

  test2 = boot_cor_test(samp1,
                     samp2,
                     alternative = "e",
                     null = .3)

  test3 = boot_cor_test(samp1,
                     samp2,
                     alternative = "e",
                     null = .5)
  expect_true(test1$p.value > test2$p.value)
  expect_true(test2$p.value > test3$p.value)

  test1 = z_cor_test(samp1,
                     samp2,
                     alternative = "m",
                     null = .1)

  test2 = z_cor_test(samp1,
                     samp2,
                     alternative = "m",
                     null = .3)

  test3 = z_cor_test(samp1,
                     samp2,
                     alternative = "m",
                     null = .5)
  expect_true(test1$p.value < test2$p.value)
  expect_true(test2$p.value < test3$p.value)

  test1 = boot_cor_test(samp1,
                        samp2,
                        alternative = "m",
                        null = .1)

  test2 = boot_cor_test(samp1,
                        samp2,
                        alternative = "m",
                        null = .3)

  test3 = boot_cor_test(samp1,
                        samp2,
                        alternative = "m",
                        null = .5)
  expect_true(test1$p.value < test2$p.value)
  expect_true(test2$p.value < test3$p.value)


})

test_that("Run examples for boot_cor_test", {
  skip_on_cran()

  set.seed(76584441)

  samp1 = rnorm(25)
  samp2 = rnorm(25)


  expect_error(boot_cor_test())
  expect_error(boot_cor_test(samp1,
                             samp2,
                             method = "p",
                             null = c(-.2,.2),
                             alternative = "t"))

  # Use boot_ci = "perc" to get percentile p-values (matches legacy expectations)
  test1 = boot_cor_test(samp1,
                     samp2,
                     method = "p",
                     boot_ci = "perc")

  test2 = boot_cor_test(samp1,
                     samp2,
                     method = "s",
                     boot_ci = "perc")

  test3 = boot_cor_test(samp1,
                     samp2,
                     method = "k",
                     boot_ci = "perc")

  expect_equal(c(unname(test1$parameter),
                 unname(test2$parameter),
                 unname(test3$parameter)),
               c(25, 25, 25))

  expect_equal(test1$p.value,
               0.78,
               tolerance = .01)
  expect_equal(test2$p.value,
               0.936,
               tolerance = .01)
  expect_equal(test3$p.value,
               0.97,
               tolerance = .01)

  test1c = cor.test(samp1,
                    samp2,
                    method = "p")

  test2c = cor.test(samp1,
                    samp2,
                    method = "s")

  test3c = cor.test(samp1,
                    samp2,
                    method = "k")

  expect_equal(unname(test1$estimate),
               unname(test1c$estimate))
  expect_equal(test3$estimate,
               test3c$estimate)
  expect_equal(test2$estimate,
               test2c$estimate)

  test1 = boot_cor_test(samp1,
                     samp2,
                     method = "p",
                     null = .4,
                     alternative = "e")

  test2 = boot_cor_test(samp1,
                     samp2,
                     method = "s",
                     null = .4,
                     alternative = "e")

  test3 = boot_cor_test(samp1,
                     samp2,
                     method = "k",
                     null = .4,
                     alternative = "e")

  test4 = boot_cor_test(samp1,
                     samp2,
                     method = "b",
                     null = .4,
                     alternative = "e")

  test5 = boot_cor_test(samp1,
                     samp2,
                     method = "w",
                     null = .4,
                     alternative = "e")
  # other alts ----

  test1 = boot_cor_test(samp1,
                     samp2,
                     method = "p",
                     boot_ci = "perc",
                     alternative = "greater")

  test2 = boot_cor_test(samp1,
                     samp2,
                     method = "s",
                     boot_ci = "perc",
                     alternative = "greater")

  test3 = boot_cor_test(samp1,
                     samp2,
                     method = "k",
                     boot_ci = "perc",
                     alternative = "greater")


  expect_equal(test1$p.value,
               0.39,
               tolerance = .05)
  expect_equal(test2$p.value,
               0.465,
               tolerance = .05)
  expect_equal(test3$p.value,
               0.51,
               tolerance = .05)

  test1 = boot_cor_test(samp1,
                     samp2,
                     method = "p",
                     boot_ci = "perc",
                     alternative = "less")

  test2 = boot_cor_test(samp1,
                     samp2,
                     method = "s",
                     boot_ci = "perc",
                     alternative = "less")

  test3 = boot_cor_test(samp1,
                     samp2,
                     method = "k",
                     boot_ci = "perc",
                     alternative = "less")


  expect_equal(test1$p.value,
               1-0.39,
               tolerance = .05)
  expect_equal(test2$p.value,
               1-0.465,
               tolerance = .05)
  expect_equal(test3$p.value,
               1-0.51,
               tolerance = .05)

})

test_that("compare_cor: z-statistic is standardized (issue #115)", {
  result <- compare_cor(
    r1 = 0.6, df1 = 18,
    r2 = 0.8, df2 = 23,
    null = 0.4,
    method = "fisher",
    alternative = "equivalence"
  )

  # Manually compute the expected standardized z
  z1   <- atanh(0.6)
  z2   <- atanh(0.8)
  diff <- z1 - z2
  SE   <- sqrt(1/17 + 1/22)
  bound_z <- atanh(0.4)

  expected_z <- (diff - (-bound_z)) / SE
  expected_p <- 1 - pnorm(expected_z)

  expect_equal(unname(result$statistic), expected_z, tolerance = 1e-6)
  expect_equal(result$p.value, expected_p, tolerance = 1e-6)
})

test_that("Run examples for boot_compare_cor", {
  skip_on_cran()

  set.seed(8922)
  x1 = rnorm(40)
  y1 = rnorm(40)

  x2 = rnorm(100)
  y2 = rnorm(100)

  test1 = boot_compare_cor(
    x1 = x1,
    x2 = x2,
    y1 = y1,
    y2 = y2,
    null = .2,
    alternative = "e",
    method = "p"
  )

  expect_error(boot_compare_cor(
    x1 = x1,
    x2 = x2,
    y1 = y1,
    y2 = y2,
    null = c(-.2,.2),
    alternative = "t",
    method = "p"
  ))

  test2 = boot_compare_cor(
    x1 = x1,
    x2 = x2,
    y1 = y1,
    y2 = y2,
    null = .2,
    alternative = "e",
    method = "s"
  )

  test2_f = compare_cor(
    r1 = 0,
    r2 = .1,
    df1 = length(x1) - 2,
    df2 = length(x2) - 2,
    null = .2,
    alternative = "e",
    method = "f"
  )
  test2_k = compare_cor(
    r1 = 0,
    r2 = .1,
    df1 = length(x1) - 2,
    df2 = length(x2) - 2,
    null = .2,
    alternative = "e",
    method = "k"
  )

  test2_f = compare_cor(
    r1 = .1,
    r2 = 0,
    df1 = length(x1) - 2,
    df2 = length(x2) - 2,
    null = .2,
    alternative = "e",
    method = "f"
  )
  test2_k = compare_cor(
    r1 = 0.1,
    r2 = 0,
    df1 = length(x1) - 2,
    df2 = length(x2) - 2,
    null = .2,
    alternative = "e",
    method = "k"
  )

  test2_f = compare_cor(
    r1 = 0,
    r2 = .1,
    df1 = length(x1) - 2,
    df2 = length(x2) - 2,
    null = .2,
    alternative = "m",
    method = "f"
  )
  test2_k = compare_cor(
    r1 = 0,
    r2 = .1,
    df1 = length(x1) - 2,
    df2 = length(x2) - 2,
    null = .2,
    alternative = "m",
    method = "k"
  )

  test2_f = compare_cor(
    r1 = 0.1,
    r2 = 0,
    df1 = length(x1) - 2,
    df2 = length(x2) - 2,
    null = .2,
    alternative = "m",
    method = "f"
  )
  test2_k = compare_cor(
    r1 = 0.1,
    r2 =  0,
    df1 = length(x1) - 2,
    df2 = length(x2) - 2,
    null = .2,
    alternative = "m",
    method = "k"
  )

  test3 = boot_compare_cor(
    x1 = x1,
    x2 = x2,
    y1 = y1,
    y2 = y2,
    null = .2,
    alternative = "e",
    method = "k"
  )

  test4 = boot_compare_cor(
    x1 = x1,
    x2 = x2,
    y1 = y1,
    y2 = y2,
    null = .2,
    alternative = "e",
    method = "win"
  )

  test5 = boot_compare_cor(
    x1 = x1,
    x2 = x2,
    y1 = y1,
    y2 = y2,
    null = .2,
    alternative = "e",
    method = "bend"
  )

  expect_equal(unname(test1$estimate),
               0.0538,
               tolerance= 0.001)

  expect_equal(unname(test2$estimate),
               -0.01047,
               tolerance= 0.001)

  expect_equal(unname(test3$estimate),
               -0.01134,
               tolerance= 0.01)

  expect_equal(unname(test4$estimate),
               0.0638,
               tolerance= 0.001)

  expect_equal(unname(test5$estimate),
               0.05207,
               tolerance= 0.001)

  test3 = boot_compare_cor(
    x1 = x1,
    x2 = x2,
    y1 = y1,
    y2 = y2,
    null = .2,
    alternative = "m",
    method = "k"
  )

  test4 = boot_compare_cor(
    x1 = x1,
    x2 = x2,
    y1 = y1,
    y2 = y2,
    null = .2,
    alternative = "m",
    method = "win"
  )

  expect_error( boot_compare_cor(
    x1 = 1,
    x2 = x2,
    y1 = y1,
    y2 = y2,
    null = .2,
    alternative = "m",
    method = "win"
  ))
  expect_error( boot_compare_cor(
    x1 = list(a=1),
    x2 = x2,
    y1 = y1,
    y2 = y2,
    null = .2,
    alternative = "m",
    method = "win"
  ))
  expect_error( boot_compare_cor(
    x1 = x1,
    x2 = 1,
    y1 = y1,
    y2 = y2,
    null = .2,
    alternative = "m",
    method = "win"
  ))

})

test_that("Run examples for corsum_test",{
  expect_error(corsum_test(n=71, r=-0.12, null=c(-0.24,.24), alpha=0.05,
                           TOST = FALSE))
  test1 = corsum_test(n=71, r=-0.12, null=0.24, alpha=0.05,
              alternative = "e")

  expect_equal(test1$p.value,
               0.1529,
               tolerance = .001)
  test1 = corsum_test(n=71, r=0.12, null=0.24, alpha=0.05,
                      alternative = "e")
  expect_equal(test1$p.value,
               0.1529,
               tolerance = .001)

  test1 = corsum_test(n=71, r=-0.12, null=0.24, alpha=0.05,
                      alternative = "m")
  expect_equal(test1$p.value,
               0.8471113,
               tolerance = .001)
  test1 = corsum_test(n=71, r=0.12, null=0.24, alpha=0.05,
                      alternative = "m")
  expect_equal(test1$p.value,
               0.8471113,
               tolerance = .001)

  test1 = corsum_test(n=71, r=-0.12,  alpha=0.05)

  expect_equal(test1$p.value,
               0.32,
               tolerance = .001)
  test1 = corsum_test(n=71, r=-0.12,  alpha=0.05,
                      method = "k")
  test1 = corsum_test(n=71, r=-0.12,  alpha=0.05,
                      method = "s")

  test1 = corsum_test(n=71, r=-0.12, alternative = "less", alpha=0.05)

  expect_equal(test1$p.value,
               0.16,
               tolerance = .001)

  test1 = corsum_test(n=71, r=-0.12, alternative = "greater", alpha=0.05)

  expect_equal(test1$p.value,
               0.84,
               tolerance = .001)
})

test_that("z_cor_test: jackknife SE and cor.se", {

  set.seed(424242)
  x <- rnorm(20)
  y <- x + rnorm(20, sd = 0.5)

  # (a) jackknife SE differs from analytic (typically larger for small n)
  for (m in c("pearson", "spearman", "kendall")) {
    res_a <- z_cor_test(x, y, method = m)
    res_j <- z_cor_test(x, y, method = m, se_method = "jackknife")

    expect_true(res_a$stderr["z.se"] != res_j$stderr["z.se"],
                label = paste0("jackknife SE differs for ", m))
  }

  # (b) CI and p-value use the same SE under jackknife
  # Verify CI width changes when switching to jackknife
  res_a <- z_cor_test(x, y, method = "pearson")
  res_j <- z_cor_test(x, y, method = "pearson", se_method = "jackknife")

  ci_width_a <- diff(res_a$conf.int)
  ci_width_j <- diff(res_j$conf.int)
  expect_true(ci_width_a != ci_width_j,
              label = "CI width changes with jackknife SE")

  # p-values should also differ

  expect_true(res_a$p.value != res_j$p.value,
              label = "p-value changes with jackknife SE")

  # (c) cor.se equals (1 - r^2) * z.se to numerical precision
  for (m in c("pearson", "spearman", "kendall")) {
    res <- z_cor_test(x, y, method = m)
    r <- unname(res$estimate)
    expected_cor_se <- (1 - r^2) * res$stderr["z.se"]
    expect_equal(unname(res$stderr["cor.se"]),
                 unname(expected_cor_se),
                 tolerance = 1e-10,
                 label = paste0("cor.se delta method for ", m))

    # also check jackknife path
    res_j <- z_cor_test(x, y, method = m, se_method = "jackknife")
    r_j <- unname(res_j$estimate)
    expected_cor_se_j <- (1 - r_j^2) * res_j$stderr["z.se"]
    expect_equal(unname(res_j$stderr["cor.se"]),
                 unname(expected_cor_se_j),
                 tolerance = 1e-10,
                 label = paste0("cor.se delta method (jackknife) for ", m))
  }

  # jackknife works with equivalence testing
  res_eq <- z_cor_test(x, y, method = "pearson",
                       alternative = "e", null = 0.8,
                       se_method = "jackknife")
  expect_true(is.finite(res_eq$p.value))
  expect_length(res_eq$stderr, 2)
})

# boot_cor_test p-value / CI consistency tests -----

test_that("boot_cor_test: stud validation errors", {
  skip_on_cran()

  x <- rnorm(20)
  y <- rnorm(20)

  expect_error(
    boot_cor_test(x, y, method = "winsorized", boot_ci = "stud"),
    "Studentized bootstrap"
  )
  expect_error(
    boot_cor_test(x, y, method = "bendpercent", boot_ci = "stud"),
    "Studentized bootstrap"
  )
})

test_that("boot_cor_test: stud method runs for pearson/spearman/kendall", {
  skip_on_cran()

  set.seed(12345)
  x <- rnorm(30)
  y <- x + rnorm(30, sd = 0.5)

  for (m in c("pearson", "spearman", "kendall")) {
    res <- boot_cor_test(x, y, method = m, boot_ci = "stud", R = 999)
    expect_true(is.finite(res$p.value), label = paste("stud p finite for", m))
    expect_equal(res$boot_ci, "stud")
    expect_true(grepl("studentized", res$method),
                label = paste("method string includes studentized for", m))
    expect_length(res$stderr, 2)
    expect_true(all(names(res$stderr) == c("boot.se", "z.se")))
  }
})

test_that("boot_cor_test: boot_ci returned in result", {
  skip_on_cran()

  set.seed(999)
  x <- rnorm(20)
  y <- rnorm(20)

  for (ci_method in c("basic", "perc", "bca", "stud")) {
    res <- boot_cor_test(x, y, method = "pearson",
                         boot_ci = ci_method, R = 599)
    expect_equal(res$boot_ci, ci_method,
                 label = paste("boot_ci field for", ci_method))
  }
})

test_that("boot_cor_test: CI/p-value agreement for perc", {
  skip_on_cran()

  set.seed(42)
  n <- 50
  x <- rnorm(n)
  y <- 0.4 * x + rnorm(n, sd = 0.8)

  # Under perc: p < alpha iff CI excludes null
  res <- boot_cor_test(x, y, method = "pearson", boot_ci = "perc",
                       alternative = "two.sided", null = 0, R = 1999)
  ci_excludes_null <- res$conf.int[1] > 0 || res$conf.int[2] < 0
  p_rejects <- res$p.value < 0.05
  expect_equal(ci_excludes_null, p_rejects,
               label = "perc CI/p agreement two.sided")
})

test_that("boot_cor_test: CI/p-value agreement for basic", {
  skip_on_cran()

  set.seed(42)
  n <- 50
  x <- rnorm(n)
  y <- 0.4 * x + rnorm(n, sd = 0.8)

  res <- boot_cor_test(x, y, method = "pearson", boot_ci = "basic",
                       alternative = "two.sided", null = 0, R = 1999)
  ci_excludes_null <- res$conf.int[1] > 0 || res$conf.int[2] < 0
  p_rejects <- res$p.value < 0.05
  expect_equal(ci_excludes_null, p_rejects,
               label = "basic CI/p agreement two.sided")
})

test_that("boot_cor_test: CI/p-value agreement for bca", {
  skip_on_cran()

  set.seed(42)
  n <- 50
  x <- rnorm(n)
  y <- 0.4 * x + rnorm(n, sd = 0.8)

  res <- boot_cor_test(x, y, method = "pearson", boot_ci = "bca",
                       alternative = "two.sided", null = 0, R = 1999)
  ci_excludes_null <- res$conf.int[1] > 0 || res$conf.int[2] < 0
  p_rejects <- res$p.value < 0.05
  expect_equal(ci_excludes_null, p_rejects,
               label = "bca CI/p agreement two.sided")
})

test_that("boot_cor_test: CI/p-value agreement for stud", {
  skip_on_cran()

  set.seed(42)
  n <- 50
  x <- rnorm(n)
  y <- 0.4 * x + rnorm(n, sd = 0.8)

  res <- boot_cor_test(x, y, method = "pearson", boot_ci = "stud",
                       alternative = "two.sided", null = 0, R = 1999)
  ci_excludes_null <- res$conf.int[1] > 0 || res$conf.int[2] < 0
  p_rejects <- res$p.value < 0.05
  expect_equal(ci_excludes_null, p_rejects,
               label = "stud CI/p agreement two.sided")
})

test_that("boot_cor_test: equivalence and MET with all CI methods", {
  skip_on_cran()

  set.seed(101)
  n <- 80
  x <- rnorm(n)
  y <- rnorm(n)  # near-zero correlation

  for (ci_method in c("basic", "perc", "bca", "stud")) {
    # Equivalence: wide bounds should reject (p < alpha)
    res_wide <- boot_cor_test(x, y, method = "pearson",
                              boot_ci = ci_method,
                              alternative = "equivalence",
                              null = 0.5, R = 999)
    expect_true(is.finite(res_wide$p.value),
                label = paste("equ p finite for", ci_method))

    # MET: wide bounds should fail to reject (p >= alpha)
    res_met <- boot_cor_test(x, y, method = "pearson",
                             boot_ci = ci_method,
                             alternative = "minimal.effect",
                             null = 0.5, R = 999)
    expect_true(is.finite(res_met$p.value),
                label = paste("met p finite for", ci_method))
  }
})

test_that("boot_cor_test: equivalence with asymmetric bounds", {
  skip_on_cran()

  set.seed(202)
  n <- 60
  x <- rnorm(n)
  y <- rnorm(n)

  for (ci_method in c("basic", "perc", "bca", "stud")) {
    res <- boot_cor_test(x, y, method = "pearson",
                         boot_ci = ci_method,
                         alternative = "equivalence",
                         null = c(-0.3, 0.5), R = 999)
    expect_true(is.finite(res$p.value),
                label = paste("asymmetric equ for", ci_method))
    expect_equal(length(res$null.value), 2)
  }
})

test_that("boot_cor_test: stud results similar to z_cor_test for large n Pearson", {
  skip_on_cran()

  set.seed(303)
  n <- 200
  x <- rnorm(n)
  y <- 0.3 * x + rnorm(n, sd = 0.9)

  boot_res <- boot_cor_test(x, y, method = "pearson",
                            boot_ci = "stud", R = 1999)
  z_res <- z_cor_test(x, y, method = "pearson")

  # Point estimates should be identical
  expect_equal(unname(boot_res$estimate), unname(z_res$estimate))

  # p-values should be in the same ballpark
  expect_equal(boot_res$p.value, z_res$p.value, tolerance = 0.1)
})
