## -----------------------------------------------------------------------------
data('sleep')
library(TOSTER)

test1 = wilcox_TOST(formula = extra ~ group,
                      data = sleep,
                      paired = FALSE,
                      eqb = .5)


print(test1)

## -----------------------------------------------------------------------------
# Rank biserial
wilcox_TOST(formula = extra ~ group,
                      data = sleep,
                      paired = FALSE,
                      ses = "r",
                      eqb = .5)

# Odds

wilcox_TOST(formula = extra ~ group,
                      data = sleep,
                      paired = FALSE,
                      ses = "o",
                      eqb = .5)

# Concordance

wilcox_TOST(formula = extra ~ group,
                      data = sleep,
                      paired = FALSE,
                      ses = "c",
                      eqb = .5)


## -----------------------------------------------------------------------------
# studentized test
brunner_munzel(formula = extra ~ group,
                      data = sleep,
                      paired = FALSE)
# permutation
brunner_munzel(formula = extra ~ group,
                      data = sleep,
                      paired = FALSE,
               perm = TRUE)

## -----------------------------------------------------------------------------
# permutation based Brunner-Munzel test of equivalence
simple_htest(formula = extra ~ group,
             test = "brunner",
             data = sleep,
             paired = FALSE,
             alternative = "equ",
             mu = .7,
             perm = TRUE)


## -----------------------------------------------------------------------------
data('sleep')

test1 = boot_t_TOST(formula = extra ~ group,
                      data = sleep,
                      paired = TRUE,
                      eqb = .5,
                    R = 499)


print(test1)

plot(test1)

## -----------------------------------------------------------------------------
log_TOST(
  mpg ~ am,
  data = mtcars
)

## -----------------------------------------------------------------------------
boot_log_TOST(
  mpg ~ am,
  data = mtcars,
  R = 499
)

## -----------------------------------------------------------------------------
ses_calc(formula = extra ~ group,
         data = sleep,
         paired = TRUE,
         ses = "r")

# Setting bootstrap replications low to
## reduce compiling time of vignette
boot_ses_calc(formula = extra ~ group,
         data = sleep,
         paired = TRUE,
         R = 199,
         boot_ci = "perc", # recommend percentile bootstrap for paired SES
         ses = "r") 

