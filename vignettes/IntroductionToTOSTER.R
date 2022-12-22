## ---- fig.width=6-------------------------------------------------------------
library("TOSTER")
TOSTmeta(ES = 0.06, se = 0.003, low_eqbound_d=-0.1, high_eqbound_d=0.1, alpha=0.05)

## ---- fig.width=6-------------------------------------------------------------
# OLD CODE
#TOSTtwo(m1=5.25,m2=5.22,sd1=0.95,sd2=0.83,n1=95,n2=89,low_eqbound_d=-0.48, high_eqbound=0.48, alpha = 0.05, var.equal=TRUE)
TOSTtwo.raw(m1=5.25,m2=5.22,sd1=0.95,sd2=0.83,n1=95,n2=89,low_eqbound=-0.429, high_eqbound=0.429, alpha = 0.05, var.equal=TRUE)

# NEW CODE
tsum_TOST(m1=5.25,m2=5.22,sd1=0.95,sd2=0.83,n1=95,n2=89,eqb=0.429, alpha = 0.05, var.equal=TRUE)

## ---- fig.width=6-------------------------------------------------------------
# OLD CODE
TOSTtwo(m1=100.64,m2=100.48,sd1=14.1,sd2=14.9,n1=39343,n2=40033,low_eqbound_d=-0.05, high_eqbound_d=0.05, alpha = 0.05, var.equal=FALSE)

# NEW CODE
tsum_TOST(m1=100.64,m2=100.48,sd1=14.1,sd2=14.9,n1=39343,n2=40033,eqb=0.05, alpha = 0.05, var.equal=FALSE,
          eqbound_type = "SMD")

## -----------------------------------------------------------------------------
# OLD CODE
powerTOSTone(alpha=0.05, statistical_power=0.8, low_eqbound_d=-0.68, high_eqbound_d=0.68)

# NEW CODE -- note there is a minor discrepancy
# because the new function uses a different solution for power
power_t_TOST(type = "one.sample",eqb = 0.68,
             power = 0.8,alpha=.05)

## ---- fig.width=6-------------------------------------------------------------
# OLD CODE
TOSTone(m=5.71,mu=6,sd=1.79,n=20,low_eqbound_d=-0.68, high_eqbound_d=0.68, alpha=0.05)

# NEW CODE
tsum_TOST(m1=5.71-6,sd1=1.79,n1=20,eqb=0.68, eqbound_type = "SMD")

## ---- fig.width=6-------------------------------------------------------------
# No new code for correlations (yet)
powerTOSTr(alpha=0.05, statistical_power=0.8, low_eqbound_r=-0.24, high_eqbound_r=0.24)


## ---- fig.width=6-------------------------------------------------------------
# OLD CODE
TOSTr(n=71, r=-0.12, low_eqbound_r=-0.24, high_eqbound_r=0.24, alpha=0.05)

# NEW CODE
corsum_test(n=71, r=-0.12, null=0.24, alpha=0.05,
            alternative = "equivalence")

