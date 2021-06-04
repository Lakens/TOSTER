## -----------------------------------------------------------------------------
data('sleep')

test1 = wilcox_TOST(formula = extra ~ group,
                      data = sleep,
                      paired = FALSE,
                      low_eqbound = -.5,
                      high_eqbound = .5)


print(test1)

## -----------------------------------------------------------------------------
data('sleep')

test1 = boot_t_TOST(formula = extra ~ group,
                      data = sleep,
                      paired = TRUE,
                      low_eqbound = -.5,
                      high_eqbound = .5,
                    R = 999)


print(test1)

plot(test1)

