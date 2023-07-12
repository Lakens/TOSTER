library(TOSTER)
library(nparcomp)
library(tidyverse)
sleep2 = sleep|>
  mutate(group = factor(group,
                        levels = c(2,1)))
test1 = brunner_munzel(data = sleep, extra ~ group)
test1_comp = npar.t.test(data = sleep2, extra ~ group, method = "t.app")
summary(test1_comp)
test1

test2_t = brunner_munzel(data = sleep, extra ~ group,
                       paired = TRUE,
                       perm = FALSE)

test2_p = brunner_munzel(data = sleep, extra ~ group,
                       paired = TRUE,
                       perm = TRUE)

test2_comp = npar.t.test.paired(data = sleep2, extra ~ group)
summary(test2_comp)

npar.t.test.paired(data = sleep2, extra ~ group)


test1_comp = npar.t.test(data = sleep2, extra ~ group, method = "permu",
                         alternative = "two.sided",
                         nperm = 20000)
test1_comp$Analysis
test1 = brunner_munzel(data = sleep,
                       x=x,
                       y=y,
                       perm = TRUE,
                       max_n_perm = 20000)
test1
