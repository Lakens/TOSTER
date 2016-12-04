## ---- fig.width=6--------------------------------------------------------
library("TOSTER")
TOSTtwo.raw(m1=5.25,m2=5.22,sd1=0.95,sd2=0.83,n1=95,n2=89,low_eqbound=-0.384, high_eqbound=0.384, alpha = 0.05, var.equal=TRUE)

## ---- fig.width=6--------------------------------------------------------
library("TOSTER")
TOSTtwo(m1=100.64,m2=100.48,sd1=14.1,sd2=14.9,n1=39343,n2=40033,low_eqbound_d=-0.05, high_eqbound_d=0.05, alpha = 0.05, var.equal=FALSE)

## ------------------------------------------------------------------------
library("TOSTER")
powerTOSTone(alpha=0.05, statistical_power=0.8, low_eqbound_d=-0.68, high_eqbound_d=0.68)

## ---- fig.width=6--------------------------------------------------------
library("TOSTER")
TOSTone(m=6.05,mu=6,sd=1.5,n=21,low_eqbound_d=-0.68, high_eqbound_d=0.68, alpha=0.05)

TOSTone(m=6.68,mu=6,sd=1.7,n=19,low_eqbound_d=-0.68, high_eqbound_d=0.68, alpha=0.05)

TOSTone(m=5.95,mu=6,sd=1.4,n=20,low_eqbound_d=-0.68, high_eqbound_d=0.68, alpha=0.05)

TOSTone(m=6.45,mu=6,sd=1.91,n=20,low_eqbound_d=-0.68, high_eqbound_d=0.68, alpha=0.05)

TOSTone(m=5.71,mu=6,sd=1.79,n=20,low_eqbound_d=-0.68, high_eqbound_d=0.68, alpha=0.05)

## ---- fig.width=6--------------------------------------------------------
library("TOSTER")
TOSTone.raw(m=6.05,mu=6,sd=1.5,n=21,low_eqbound=-1, high_eqbound=1, alpha=0.05)

TOSTone.raw(m=6.68,mu=6,sd=1.7,n=19,low_eqbound=-1, high_eqbound=1, alpha=0.05)

TOSTone.raw(m=5.95,mu=6,sd=1.4,n=20,low_eqbound=-1, high_eqbound=1, alpha=0.05)

TOSTone.raw(m=6.45,mu=6,sd=1.91,n=20,low_eqbound=-1, high_eqbound=1, alpha=0.05)

TOSTone.raw(m=5.71,mu=6,sd=1.79,n=20,low_eqbound=-1, high_eqbound=1, alpha=0.05)

