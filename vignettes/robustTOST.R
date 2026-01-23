## -----------------------------------------------------------------------------
data('sleep')
library(TOSTER)

test1 = wilcox_TOST(formula = extra ~ group,
                      data = sleep,
                      eqb = .5)


print(test1)

## -----------------------------------------------------------------------------
# Rank biserial
wilcox_TOST(formula = extra ~ group,
                      data = sleep,
                      ses = "r",
                      eqb = .5)

# Odds

wilcox_TOST(formula = extra ~ group,
                      data = sleep,
                      ses = "o",
                      eqb = .5)

# Concordance

wilcox_TOST(formula = extra ~ group,
                      data = sleep,
                      ses = "c",
                      eqb = .5)


## -----------------------------------------------------------------------------
# Default studentized test (t-approximation)
brunner_munzel(formula = extra ~ group,
               data = sleep)

## -----------------------------------------------------------------------------
# Permutation test (recommended for small samples)
brunner_munzel(formula = extra ~ group,
               data = sleep,
               test_method = "perm")

## -----------------------------------------------------------------------------
# Logit transformation for range-preserving CIs
brunner_munzel(formula = extra ~ group,
               data = sleep,
               test_method = "logit")

## -----------------------------------------------------------------------------
# Equivalence test: is the relative effect between 0.3 and 0.7?
brunner_munzel(formula = extra ~ group,
               data = sleep,
               alternative = "equivalence",
               mu = c(0.3, 0.7))

## -----------------------------------------------------------------------------
# Permutation-based equivalence test
brunner_munzel(formula = extra ~ group,
               data = sleep,
               alternative = "equivalence",
               mu = c(0.3, 0.7),
               test_method = "perm")

## -----------------------------------------------------------------------------
# Minimal effect test: is the relative effect outside 0.4 to 0.6?
brunner_munzel(formula = extra ~ group,
               data = sleep,
               alternative = "minimal.effect",
               mu = c(0.4, 0.6))

## -----------------------------------------------------------------------------
# Paired samples test
brunner_munzel(x = sleep$extra[sleep$group == 1],
               y = sleep$extra[sleep$group == 2],
               paired = TRUE)

## -----------------------------------------------------------------------------
# Paired samples with permutation test
brunner_munzel(x = sleep$extra[sleep$group == 1],
               y = sleep$extra[sleep$group == 2],
               paired = TRUE,
               test_method = "perm")

## -----------------------------------------------------------------------------
# Test if the relative effect differs from 0.3
brunner_munzel(formula = extra ~ group,
               data = sleep,
               mu = 0.3)

## -----------------------------------------------------------------------------
# Example: test non-inferiority with lower bound of 0.35
# (i.e., the new treatment should not be substantially worse)
brunner_munzel(formula = extra ~ group,
               data = sleep,
               alternative = "greater",
               mu = 0.35)

## -----------------------------------------------------------------------------
data('sleep')

# Two-sample permutation t-test
perm_result <- perm_t_test(extra ~ group, 
                           data = sleep,
                           R = 1999)
perm_result

## -----------------------------------------------------------------------------
# Simulate data with outliers
set.seed(42)
x <- c(rnorm(18, mean = 0), 8, 12)  # Two outliers
y <- rnorm(20, mean = 0)

# Standard permutation test (sensitive to outliers)
perm_t_test(x, y, R = 999)

# Trimmed permutation test (robust to outliers)
perm_t_test(x, y, tr = 0.1, R = 999)

## -----------------------------------------------------------------------------
# Equivalence test: is the effect within ±3 units?
perm_t_test(extra ~ group, 
            data = sleep,
            alternative = "equivalence",
            mu = c(-3, 3),
            R = 999)

## -----------------------------------------------------------------------------
# Small paired sample - exact permutations will be computed
before <- c(5.1, 4.8, 6.2, 5.7, 6.0, 5.5, 4.9, 5.8)
after <- c(5.6, 5.2, 6.7, 6.1, 6.5, 5.8, 5.3, 6.2)

perm_t_test(x = before, y = after,
            paired = TRUE,
            alternative = "less",
            R = 999)

## -----------------------------------------------------------------------------
data('sleep')

# For paired tests with bootstrap methods, use separate vectors
test1 = boot_t_TOST(x = sleep$extra[sleep$group == 1],
                    y = sleep$extra[sleep$group == 2],
                      paired = TRUE,
                      eqb = .5,
                    R = 499)


print(test1)

plot(test1)

## -----------------------------------------------------------------------------
# Same equivalence test using both methods
data('sleep')

# Bootstrap approach
boot_result <- boot_t_test(extra ~ group, 
                           data = sleep,
                           alternative = "equivalence",
                           mu = c(-2, 2),
                           R = 999)

# Permutation approach  
perm_result <- perm_t_test(extra ~ group,
                           data = sleep,
                           alternative = "equivalence", 
                           mu = c(-2, 2),
                           R = 999)

# Compare results
boot_result
perm_result

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
# For paired tests, use separate vectors
ses_calc(x = sleep$extra[sleep$group == 1],
         y = sleep$extra[sleep$group == 2],
         paired = TRUE,
         ses = "r")

# Setting bootstrap replications low to
## reduce compiling time of vignette
boot_ses_calc(x = sleep$extra[sleep$group == 1],
              y = sleep$extra[sleep$group == 2],
         paired = TRUE,
         R = 199,
         boot_ci = "perc", # recommend percentile bootstrap for paired SES
         ses = "r") 

