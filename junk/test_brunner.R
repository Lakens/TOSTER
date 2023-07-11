
test1 = brunner_munzel(data = sleep, extra ~ group)
test1_comp = npar.t.test(data = sleep, extra ~ group, method = "t.app")
summary(test1_comp)
test1

test2 = brunner_munzel(data = sleep, extra ~ group,
                       paired = TRUE)
test2_comp = npar.t.test.paired(data = sleep, extra ~ group)
summary(test2_comp)
