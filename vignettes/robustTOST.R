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
data('sleep')

test1 = boot_t_TOST(formula = extra ~ group,
                      data = sleep,
                      paired = TRUE,
                      eqb = .5,
                    R = 999)


print(test1)

plot(test1)

