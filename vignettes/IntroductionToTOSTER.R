## ---- fig.width=6--------------------------------------------------------
library("TOSTER")
TOSTmeta(ES = 0.06, se = 0.003, low_eqbound_d=-0.1, high_eqbound_d=0.1, alpha=0.05)

## ---- fig.width=6--------------------------------------------------------
TOSTtwo.raw(m1=5.25,m2=5.22,sd1=0.95,sd2=0.83,n1=95,n2=89,low_eqbound=-0.384, high_eqbound=0.384, alpha = 0.05, var.equal=TRUE)

## ---- fig.width=6--------------------------------------------------------
TOSTtwo(m1=100.64,m2=100.48,sd1=14.1,sd2=14.9,n1=39343,n2=40033,low_eqbound_d=-0.05, high_eqbound_d=0.05, alpha = 0.05, var.equal=FALSE)

## ------------------------------------------------------------------------
powerTOSTone(alpha=0.05, statistical_power=0.8, low_eqbound_d=-0.68, high_eqbound_d=0.68)

## ---- fig.width=6--------------------------------------------------------
TOSTone(m=5.71,mu=6,sd=1.79,n=20,low_eqbound_d=-0.68, high_eqbound_d=0.68, alpha=0.05)

## ---- fig.width=6--------------------------------------------------------
powerTOSTr(alpha=0.05, statistical_power=0.8, low_eqbound_r=-0.24, high_eqbound_r=0.24)

## ---- fig.width=6--------------------------------------------------------
TOSTr(n=71, r=-0.06, low_eqbound_r=-0.24, high_eqbound_r=0.24, alpha=0.05)

