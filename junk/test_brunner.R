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

test2 = brunner_munzel(data = sleep, extra ~ group,
                       paired = TRUE)
wilcox.test(data = sleep, extra ~ group,
               paired = TRUE)
test2_comp = npar.t.test.paired(data = sleep2, extra ~ group)
summary(test2_comp)
