#context("Run Examples for boot_t_TOST")

# need hush function to run print through examples

hush = function(code) {
  sink("NUL") # use /dev/null in UNIX
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

  expect_equal(test1$estimate,
               test1c$estimate)
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

  set.seed(76584441)

  samp1 = rnorm(25)
  samp2 = rnorm(25)


  expect_error(boot_cor_test())
  expect_error(boot_cor_test(samp1,
                             samp2,
                             method = "p",
                             null = c(-.2,.2),
                             alternative = "t"))

  test1 = boot_cor_test(samp1,
                     samp2,
                     method = "p")

  test2 = boot_cor_test(samp1,
                     samp2,
                     method = "s")

  test3 = boot_cor_test(samp1,
                     samp2,
                     method = "k")

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

  expect_equal(test1$estimate,
               test1c$estimate)
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
                     alternative = "greater")

  test2 = boot_cor_test(samp1,
                     samp2,
                     method = "s",
                     alternative = "greater")

  test3 = boot_cor_test(samp1,
                     samp2,
                     method = "k",
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
                     alternative = "less")

  test2 = boot_cor_test(samp1,
                     samp2,
                     method = "s",
                     alternative = "less")

  test3 = boot_cor_test(samp1,
                     samp2,
                     method = "k",
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

test_that("Run examples for boot_compare_cor", {

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
  test1 = corsum_test(n=71, r=-0.12, null=0.24, alpha=0.05,
                      alternative = "m")

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
