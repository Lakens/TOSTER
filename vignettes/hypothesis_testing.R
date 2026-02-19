## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(TOSTER)

## -----------------------------------------------------------------------------
simple_htest(extra ~ group,
             data = sleep,
             mu = 0,
             alternative = "two.sided")

## -----------------------------------------------------------------------------
t.test(extra ~ group, data = sleep)

## -----------------------------------------------------------------------------
simple_htest(extra ~ group,
             data = sleep,
             test = "wilcox.test",
             mu = 0,
             alternative = "two.sided")

## -----------------------------------------------------------------------------
simple_htest(extra ~ group,
             data = sleep,
             mu = 2,
             alternative = "equivalence")

## -----------------------------------------------------------------------------
simple_htest(extra ~ group,
             data = sleep,
             mu = c(-1, 3),
             alternative = "equivalence")

## -----------------------------------------------------------------------------
simple_htest(extra ~ group,
             data = sleep,
             mu = -1,
             alternative = "greater")

## -----------------------------------------------------------------------------
simple_htest(extra ~ group,
             data = sleep,
             mu = 1,
             alternative = "greater")

## -----------------------------------------------------------------------------
brunner_munzel(extra ~ group, data = sleep)

## -----------------------------------------------------------------------------
boot_cor_test(mtcars$mpg, mtcars$hp,
           method = "pearson",
           alternative = "two.sided",
           null = 0)

## -----------------------------------------------------------------------------
set.seed(2101)
boot_t_test(extra ~ group,
            data = sleep,
            mu = 0,
            alternative = "two.sided",
            R = 999)

## -----------------------------------------------------------------------------
set.seed(8251)
perm_t_test(extra ~ group,
            data = sleep,
            mu = 0,
            alternative = "two.sided",
            R = 999)

## -----------------------------------------------------------------------------
tost_res <- t_TOST(extra ~ group,
                   data = sleep,
                   eqb = 2)
as_htest(tost_res)

## -----------------------------------------------------------------------------
res_t <- simple_htest(extra ~ group, data = sleep, mu = 0)
df_htest(res_t)

## -----------------------------------------------------------------------------
describe_htest(res_t)

## -----------------------------------------------------------------------------
res_equiv <- simple_htest(extra ~ group, data = sleep, 
                           mu = 2, alternative = "equivalence")
describe_htest(res_equiv)

## ----fig.width=6, fig.height=3------------------------------------------------
plot_htest_est(res_t)

## ----fig.width=6, fig.height=3------------------------------------------------
plot_htest_est(res_equiv)

## ----fig.width=6, fig.height=3------------------------------------------------
plot_htest_est(res_equiv, describe = FALSE)

## ----fig.width=6, fig.height=3------------------------------------------------
# Run equivalence test
result <- simple_htest(extra ~ group,
                       data = sleep,
                       mu = 2,
                       alternative = "equivalence")

# Tabulate
df_htest(result)

# Describe
describe_htest(result)

# Visualize
plot_htest_est(result)

