library(tidyverse)
library(TOSTER)
library(nparcomp)
library(brunnermunzel)
library(lawstat)
library(testthat)

set.seed(17766)
tester_equal = function(x,y,tol=0.01){
  test_val = abs(unname(x) - unname(y))
  expect_equal(test_val, 0,
               tolerance = tol)
}
extract_bm_2 = function(x){
  return(x$Analysis)
}

nsim = 100

# Gaussian 2-sample ------

for(i in 1:nsim){

  n_samp1 = sample(30:100, 1)
  n_samp2 = sample(30:100, 1)

  mu_samp1 = sample(-2:2, 1)
  mu_samp2 = sample(-2:2, 1)

  samp1 = rnorm(n_samp1, mu_samp1)
  samp2 = rnorm(n_samp2, mu_samp2)

  dat = data.frame(group = c(rep("x", n_samp1),
                             rep("y", n_samp2)),
                   y = c(samp1, samp2))

  toster_res = brunner_munzel(x = samp2,
                              y = samp1)

  npar_res = npar.t.test(y ~ group,
                         data = dat,
                         info = FALSE,
                         method = "t.app")

  npar2_res = extract_bm_2(npar_res)

  law_res = lawstat::brunner.munzel.test(x = samp1,
                                         y = samp2)


  tester_equal(law_res$p.value, toster_res$p.value)
  tester_equal(law_res$estimate, toster_res$estimate)
  tester_equal(ifelse(law_res$conf.int[1] < 0,
                      0,
                      law_res$conf.int[1] ), toster_res$conf.int[1])
  tester_equal(ifelse(law_res$conf.int[2] > 1,
                      1,
                      law_res$conf.int[2]), toster_res$conf.int[2])

  bm_res = brunnermunzel.test(x = samp1,y = samp2)

  tester_equal(bm_res$p.value, toster_res$p.value)
  tester_equal(bm_res$estimate, toster_res$estimate)
  tester_equal(ifelse(bm_res$conf.int[1] < 0,
                      0,
                      bm_res$conf.int[1] ), toster_res$conf.int[1])
  tester_equal(ifelse(bm_res$conf.int[2] > 1,
                      1,
                      bm_res$conf.int[2]), toster_res$conf.int[2])

}


