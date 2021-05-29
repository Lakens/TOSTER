
test_that("Equivalence F-tests",{

  #############################################################
  ## This reads in the data, and formats it appropriately:
  data('hawthorne')
  side_data = hawthorne

  t1 = equ_ftest(Fstat = 0.49348,
                 df1 = 2,
                 df2 = 4577,
                 eqbound = .01)
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

  expect_error(equ_anova(lmmodel))

  anovamodel = anova(lmmodel)
  t2 = equ_anova(anovamodel , eqbound = .01)$p.equ
  expect_equal(t2, t1$p.value)

  aovmodel = aov(totaldrinking.diff  ~ group , data= side_data)

  t3 = equ_anova(aovmodel, eqbound = .01)$p.equ
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

})



