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
    expect_equal(res$boot_scale, "z")

    res_r <- boot_cor_test(x, y, method = m, boot_ci = "stud", R = 499,
                           boot_scale = "r")
    expect_true(all(names(res_r$stderr) == c("boot.se", "r.se")))
    expect_equal(res_r$boot_scale, "r")
    expect_true(all(res_r$conf.int >= -1 & res_r$conf.int <= 1))
  }
})

# Influence-function SEs for studentized bootstrap -----

test_that(".cor_se Pearson is the HC4-corrected fourth-moment (ADF) SE", {
  set.seed(11)
  n <- 60
  x <- rexp(n)
  y <- 0.5 * x + rt(n, 5)
  zx <- (x - mean(x)) / sqrt(mean((x - mean(x))^2))
  zy <- (y - mean(y)) / sqrt(mean((y - mean(y))^2))
  r <- cor(x, y)

  # uncorrected influence values reproduce the ADF variance formula
  psi0 <- TOSTER:::.cor_if_pearson(x, y)
  m <- function(a, b) mean(zx^a * zy^b)
  v <- (r^2 / 4) * (m(4, 0) + m(0, 4) + 2 * m(2, 2)) -
    r * (m(3, 1) + m(1, 3)) + m(2, 2)
  expect_equal(sqrt(sum(psi0^2)) / n, sqrt(v / n))

  # HC4 leverage correction
  X <- cbind(x, y)
  h <- 1 / n + mahalanobis(X, colMeans(X), cov(X)) / (n - 1)
  expect_equal(sum(h), 3)
  d <- pmin(4, n * h / 3)
  expect_equal(TOSTER:::.cor_se(x, y, "pearson"),
               sqrt(sum((psi0 * (1 - h)^(-d / 2))^2)) / n)
  expect_gt(TOSTER:::.cor_se(x, y, "pearson"), sqrt(v / n))
})

test_that(".cor_se Spearman without ties matches Croux-Dehon influence function", {
  set.seed(12)
  n <- 50
  x <- rnorm(n)
  y <- x + rnorm(n)
  u <- (rank(x) - 0.5) / n
  v <- (rank(y) - 0.5) / n
  sx <- sapply(x, function(xi) sum(v[x > xi]) + v[x == xi] / 2)
  sy <- sapply(y, function(yi) sum(u[y > yi]) + u[y == yi] / 2)
  h <- u * v + sx / n + sy / n
  psi <- 12 * (h - mean(h))
  expect_equal(TOSTER:::.cor_se(x, y, "spearman"), sqrt(sum(psi^2)) / n,
               tolerance = 1e-3)
})

test_that(".cor_se Kendall without ties matches U-statistic variance", {
  set.seed(13)
  n <- 40
  x <- rnorm(n)
  y <- x + rnorm(n)
  s <- sign(outer(x, x, "-")) * sign(outer(y, y, "-"))
  a <- rowSums(s) / (n - 1)
  expect_equal(TOSTER:::.cor_se(x, y, "kendall"),
               sqrt(4 / n * mean((a - mean(a))^2)))
})

test_that(".cor_se agrees with the jackknife SE, with and without ties", {
  jack_se <- function(x, y, m) {
    n <- length(x)
    j <- vapply(seq_len(n), function(i) cor(x[-i], y[-i], method = m), 0)
    sqrt((n - 1) / n * sum((j - mean(j))^2))
  }
  set.seed(14)
  n <- 300
  x1 <- rnorm(n)
  y1 <- 0.5 * x1 + rnorm(n) * sqrt(1 + x1^2)
  x2 <- sample(1:5, n, TRUE)
  y2 <- pmin(5, x2 + sample(0:2, n, TRUE))
  # Pearson without the HC4 leverage correction, which deliberately inflates it
  se_if <- function(x, y, m) {
    if (m == "pearson") {
      sqrt(sum(TOSTER:::.cor_if_pearson(x, y)^2)) / length(x)
    } else {
      TOSTER:::.cor_se(x, y, m)
    }
  }
  for (m in c("pearson", "spearman", "kendall")) {
    expect_equal(se_if(x1, y1, m), jack_se(x1, y1, m),
                 tolerance = 0.05, label = paste("continuous", m))
    expect_equal(se_if(x2, y2, m), jack_se(x2, y2, m),
                 tolerance = 0.05, label = paste("tied", m))
  }
  expect_gt(TOSTER:::.cor_se(x1, y1, "pearson"), se_if(x1, y1, "pearson"))
})

test_that(".cor_se reduces to normal-theory values under independence", {
  set.seed(15)
  n <- 4000
  x <- rnorm(n)
  y <- rnorm(n)
  expect_equal(TOSTER:::.cor_se(x, y, "pearson"), 1 / sqrt(n), tolerance = 0.05)
  expect_equal(TOSTER:::.cor_se(x, y, "spearman"), 1 / sqrt(n), tolerance = 0.05)
  expect_equal(TOSTER:::.cor_se(x, y, "kendall"), sqrt(4 / (9 * n)),
               tolerance = 0.05)
})

test_that("boot_cor_test: stud is studentized, not basic on the z scale", {
  skip_on_cran()

  set.seed(16)
  n <- 40
  x <- rnorm(n)
  y <- 0.4 * x + rnorm(n) * sqrt(1 + x^2)

  for (m in c("pearson", "spearman", "kendall")) {
    set.seed(1)
    res_s <- boot_cor_test(x, y, method = m, boot_ci = "stud", R = 999)
    set.seed(1)
    res_b <- boot_cor_test(x, y, method = m, boot_ci = "basic", R = 999)
    expect_equal(res_s$boot_res, res_b$boot_res)
    expect_false(isTRUE(all.equal(res_s$conf.int, res_b$conf.int,
                                  check.attributes = FALSE)),
                 label = paste("stud differs from basic for", m))
    r <- unname(res_s$estimate)
    expect_equal(unname(res_s$stderr["z.se"]),
                 TOSTER:::.cor_se(x, y, m) / (1 - r^2))
  }
})

test_that("boot_cor_test: boot_scale only affects basic and stud", {
  skip_on_cran()

  set.seed(17)
  x <- rnorm(30)
  y <- 0.5 * x + rnorm(30)

  for (ci_method in c("perc", "bca")) {
    set.seed(2)
    res_z <- boot_cor_test(x, y, boot_ci = ci_method, R = 599, boot_scale = "z")
    set.seed(2)
    res_r <- boot_cor_test(x, y, boot_ci = ci_method, R = 599, boot_scale = "r")
    expect_equal(res_z$conf.int, res_r$conf.int)
    expect_equal(res_z$p.value, res_r$p.value)
    expect_equal(res_z$boot_scale, "r")
  }

  # basic on the r scale reproduces the untransformed basic interval
  set.seed(3)
  res_r <- boot_cor_test(x, y, boot_ci = "basic", R = 599, boot_scale = "r")
  expect_equal(as.numeric(res_r$conf.int),
               TOSTER:::basic(res_r$boot_res, cor(x, y), 0.05))

  # basic on the z scale is the back-transformed z-scale interval
  set.seed(3)
  res_z <- boot_cor_test(x, y, boot_ci = "basic", R = 599, boot_scale = "z")
  expect_equal(as.numeric(res_z$conf.int),
               tanh(TOSTER:::basic(atanh(res_z$boot_res), atanh(cor(x, y)), 0.05)))
  expect_true(all(res_z$conf.int >= -1 & res_z$conf.int <= 1))
})

test_that("boot_cor_test: CI/p-value agreement on both scales", {
  skip_on_cran()

  set.seed(18)
  n <- 50
  x <- rnorm(n)
  y <- 0.3 * x + rnorm(n)

  run <- function(...) {
    set.seed(4)
    boot_cor_test(x, y, R = 999, ...)
  }
  for (sc in c("z", "r")) for (ci_method in c("basic", "stud")) {
    ci <- run(boot_ci = ci_method, boot_scale = sc)$conf.int
    # nulls just inside and just outside each limit
    for (nv in c(ci[1] + c(-1, 1) * 0.01, ci[2] + c(-1, 1) * 0.01)) {
      p <- run(boot_ci = ci_method, boot_scale = sc, null = nv)$p.value
      expect_equal(nv < ci[1] || nv > ci[2], p < 0.05,
                   label = paste(ci_method, sc, round(nv, 3)))
    }
  }
})

test_that("boot_cor_test: boot_ci = 'auto' picks stud for Pearson, bca otherwise", {
  skip_on_cran()

  set.seed(19)
  x <- rnorm(25)
  y <- 0.5 * x + rnorm(25)

  expected <- c(pearson = "stud", spearman = "bca", kendall = "bca",
                winsorized = "bca", bendpercent = "bca")
  for (m in names(expected)) {
    res <- boot_cor_test(x, y, method = m, R = 199)
    expect_equal(res$boot_ci, unname(expected[m]), label = paste("auto for", m))
  }

  # auto gives the same result as requesting the selected method directly
  set.seed(5)
  res_auto <- boot_cor_test(x, y, R = 199)
  set.seed(5)
  res_stud <- boot_cor_test(x, y, boot_ci = "stud", R = 199)
  expect_equal(res_auto$conf.int, res_stud$conf.int)
  expect_equal(res_auto$p.value, res_stud$p.value)

  # an explicit choice overrides auto
  expect_equal(boot_cor_test(x, y, boot_ci = "bca", R = 199)$boot_ci, "bca")
})

test_that("boot_cor_test: stud errors for perfect correlation", {
  skip_on_cran()
  x <- 1:20
  expect_error(boot_cor_test(x, 2 * x, boot_ci = "stud", R = 99),
               "undefined")
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

# z_cor_test missing data handling -----

test_that("z_cor_test handles missing data correctly", {

  set.seed(54321)
  x <- c(1.1, 2.3, 3.0, NA, 5.2, 6.1, 7.8, 8.0, 9.4, 10.7)
  y <- c(2.2, 3.8, NA, 7.9, 10.5, 11.6, 14.3, 15.1, 17.9, 21.0)

  # Should not error with NAs present
  for (m in c("pearson", "spearman", "kendall")) {
    res <- z_cor_test(x, y, method = m)
    expect_true(is.finite(res$p.value),
                label = paste("finite p with NAs for", m))
    expect_true(is.finite(unname(res$estimate)),
                label = paste("finite estimate with NAs for", m))
    # N should reflect complete cases only (8 of 10)
    expect_equal(unname(res$parameter), 8,
                 label = paste("N reflects complete cases for", m))
  }

  # Results should match manually removing NAs
  complete <- complete.cases(x, y)
  x_clean <- x[complete]
  y_clean <- y[complete]

  res_na <- z_cor_test(x, y, method = "pearson")
  res_clean <- z_cor_test(x_clean, y_clean, method = "pearson")

  expect_equal(res_na$p.value, res_clean$p.value)
  expect_equal(unname(res_na$estimate), unname(res_clean$estimate))
  expect_equal(unname(res_na$statistic), unname(res_clean$statistic))

  # Jackknife SE should also work with NAs
  res_jk <- z_cor_test(x, y, method = "pearson", se_method = "jackknife")
  expect_true(is.finite(res_jk$p.value),
              label = "jackknife with NAs gives finite p")
  expect_equal(unname(res_jk$parameter), 8)

  # Equivalence test with NAs
  res_eq <- z_cor_test(x, y, method = "pearson",
                       alternative = "equivalence", null = 0.99)
  expect_true(is.finite(res_eq$p.value),
              label = "equivalence with NAs gives finite p")
})

test_that("z_cor_test errors when |r| = 1 with jackknife SE", {
  # Perfect positive correlation: atanh(1) = Inf -> jackknife SE becomes NaN
  x <- 1:10
  y <- 2 * (1:10)

  # Analytic path still works (returns Inf z / p ~ 0)
  expect_silent(z_cor_test(x, y, method = "pearson"))

  # Jackknife should error with a clear message
  expect_error(
    z_cor_test(x, y, method = "pearson", se_method = "jackknife"),
    regexp = "Jackknife SE cannot be computed"
  )

  # Perfect negative correlation also triggers the guard
  expect_error(
    z_cor_test(x, -y, method = "pearson", se_method = "jackknife"),
    regexp = "Jackknife SE cannot be computed"
  )

  # Spearman/Kendall on strictly monotonic data also yield |r| = 1
  expect_error(
    z_cor_test(x, y, method = "spearman", se_method = "jackknife"),
    regexp = "Jackknife SE cannot be computed"
  )
  expect_error(
    z_cor_test(x, y, method = "kendall", se_method = "jackknife"),
    regexp = "Jackknife SE cannot be computed"
  )
})

test_that("z_cor_test errors when too few complete cases for SE formula", {
  # Pearson/Spearman require n >= 4; Kendall requires n >= 5.

  # --- Raw sample too small ---
  x3 <- c(1.1, 2.3, 4.7)
  y3 <- c(0.5, 1.9, 3.2)
  expect_error(
    z_cor_test(x3, y3, method = "pearson"),
    regexp = "Not enough complete cases"
  )
  expect_error(
    z_cor_test(x3, y3, method = "spearman"),
    regexp = "Not enough complete cases"
  )

  # n = 4 is below Kendall's threshold but OK for Pearson/Spearman
  x4 <- c(1.1, 2.3, 4.7, 6.8)
  y4 <- c(0.5, 1.9, 3.2, 5.1)
  expect_error(
    z_cor_test(x4, y4, method = "kendall"),
    regexp = "Not enough complete cases"
  )
  expect_silent(z_cor_test(x4, y4, method = "pearson"))
  expect_silent(z_cor_test(x4, y4, method = "spearman"))

  # --- NAs reducing complete cases below threshold ---
  # 10 obs but only 3 complete -> fails Pearson/Spearman
  xn <- c(1.1, 2.3, 4.7, NA, NA, NA, NA, NA, NA, NA)
  yn <- c(0.5, 1.9, 3.2, NA, NA, NA, NA, NA, NA, NA)
  expect_error(
    z_cor_test(xn, yn, method = "pearson"),
    regexp = "Not enough complete cases"
  )

  # Error message mentions method and observed N
  err <- tryCatch(z_cor_test(x3, y3, method = "kendall"),
                  error = function(e) conditionMessage(e))
  expect_match(err, "kendall")
  expect_match(err, "3 complete observation")
  expect_match(err, "at least 5")
})
