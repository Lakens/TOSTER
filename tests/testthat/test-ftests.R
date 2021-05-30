
test_that("Equivalence F-tests",{

  #############################################################
  ## This reads in the data, and formats it appropriately:
  data('hawthorne')
  side_data = hawthorne

  t1 = suppressMessages(equ_ftest(Fstat = 0.49348,
                 df1 = 2,
                 df2 = 4577,
                 eqbound = .01))
  lamb = 4580*.01 / (1-.01)
  pval = pf(
    0.49348,
    df1 = 2,
    df2 = 4577,
    ncp = lamb,
    lower.tail = TRUE
  )

  expect_equal(t1$p.value, pval)
  expect_equal(t1$p.value, 1.127e-09)

  lmmodel <- lm(totaldrinking.diff  ~ group , data= side_data)

  expect_error(suppressMessages(equ_anova(lmmodel)))

  anovamodel = anova(lmmodel)

  t2 = suppressMessages(equ_anova(anovamodel , eqbound = .01)$p.equ)
  expect_equal(t2, t1$p.value)

  aovmodel = aov(totaldrinking.diff  ~ group , data= side_data)

  t3 = suppressMessages(equ_anova(aovmodel, eqbound = .01)$p.equ)
  expect_equal(t3, t1$p.value)

  p1 = plot_pes(Fstat = 1.5,
           df1 = 2,
           df2 = 45)

  p2 = plot_pes(0.49348,
                df1 = 2,
                df2 = 4577,
                type = "cd")

  p3 = plot_pes(Fstat = 1.5,
                df1 = 2,
                df2 = 45,
                type = "c")

  p4 = plot_pes(0.49348,
                df1 = 2,
                df2 = 4577)

  # Test power function

  t1 = suppressMessages(equ_ftest(Fstat = 0.49348,
                 df1 = 2,
                 df2 = 457,
                 eqbound = .01))
  t2 = suppressMessages(power_eq_f(.0623,
                  df1 = 2,
                  df2 = 457,
                  eqbound = .01))

  expect_equal(round(t1$p.value,3),
               round(t2$sig.level,3))

  expect_equal(round(t2$power,2),
               .39)

  ##  get aov tests
  op <- options(contrasts = c("contr.helmert", "contr.poly"))
  npk.aov <- suppressMessages(aov(yield ~ block + N*P*K, npk))

  t1 =  suppressMessages(equ_anova(npk.aov, eqbound = .9))
  npk.aovE <- aov(yield ~  N*P*K + Error(block), npk)
  t2 =  suppressMessages(equ_anova(npk.aovE, eqbound = .01))

})



