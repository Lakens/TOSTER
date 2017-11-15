## ------------------------------------------------------------------------
library(knitr)
opts_chunk$set(tidy.opts=list(width.cutoff=60),tidy=TRUE)

## ------------------------------------------------------------------------
require(TOSTER)
TOSTtwo.raw(m1=1.5679,m2=1.6764,sd1=0.4285,sd2=0.4748,n1=64,n2=70,low_eqbound=-0.3136, high_eqbound=0.3136,alpha=0.05, var.equal=TRUE)


## ------------------------------------------------------------------------
TOSTtwo.raw(m1=89.3,m2=87.7,sd1=1.9,sd2=1.3,n1=12,n2=12,low_eqbound=-3.7, high_eqbound=3.7,alpha=0.05, var.equal=TRUE)


## ------------------------------------------------------------------------
TOSTtwo.raw(m1=7.1, m2=6.9, sd1=sqrt(2.0), sd2=sqrt(2.2), n1=100, n2=100,low_eqbound=-0.5, high_eqbound=0.5,alpha=0.05, var.equal=TRUE)


## ------------------------------------------------------------------------
#Example 1
TOSTtwo.raw(m1=459.09, m2=402.61, sd1=47.53, sd2=38.42, n1=5, n2=5, low_eqbound=-15, high_eqbound=15,alpha=0.05, var.equal=TRUE)

#Example 2
TOSTtwo.raw(m1=407.24, m2=402.98, sd1=41.77, sd2=41.79, n1=150, n2=150, low_eqbound=-15, high_eqbound=15,alpha=0.05, var.equal=TRUE)


## ------------------------------------------------------------------------
TOSTtwo.raw(m1=2.5,m2=2.3,sd1=1.3, sd2=1.1,n1=86,n2=102,low_eqbound=-0.4, high_eqbound=0.4,alpha=0.05)

## ------------------------------------------------------------------------
TOSTtwo.raw(m1=2.5,m2=2.3,sd1=1.3, sd2=1.1,n1=86,n2=102,low_eqbound=-0.4, high_eqbound=0.4,alpha=0.05, var.equal=TRUE)

## ------------------------------------------------------------------------
TOSTtwo.raw(m1=2.5,m2=2.3,sd1=1.3, sd2=1.1,n1=86,n2=102,low_eqbound=-0.4, high_eqbound=0.4,alpha=0.05)

## ------------------------------------------------------------------------
require(equivalence)
data(ufc)
#Remove missing data
ufc<-ufc[complete.cases(ufc[,9:10]),]
      
tost(ufc$Height.m.p, ufc$Height.m, epsilon = 1, paired = TRUE)

TOSTpaired.raw(n=length(ufc$Height.m.p),m1=mean(ufc$Height.m.p,na.rm = TRUE), m2=mean(ufc$Height.m,na.rm = TRUE), sd1=sd(ufc$Height.m.p,na.rm = TRUE), sd2=sd(ufc$Height.m,na.rm = TRUE), r12=cor(ufc$Height.m,ufc$Height.m.p, use="pairwise.complete.obs"), low_eqbound=-1,high_eqbound=1, alpha=0.05)

## ------------------------------------------------------------------------
morning<-c(3,4,4,5,4)
evening<-c(1,4,1,4,5)

TOSTpaired.raw(n=5,m1=mean(morning), m2=mean(evening), sd1=sd(morning), sd2=sd(evening), r12=cor(morning,evening), low_eqbound=-0.5,high_eqbound=0.5, alpha=0.05)

## ------------------------------------------------------------------------
DV1<-c(4,2,4,3,5,4,3,4,5,4,2,3,4,5,5)
DV2<-c(3,4,2,3,5,4,3,4,5,3,4,3,2,3,4)

TOSTpaired.raw(n=length(DV1),m1=mean(DV1), m2=mean(DV2), sd1=sd(DV1), sd2=sd(DV2), r12=cor(DV1,DV2), low_eqbound=-0.8,high_eqbound=0.8, alpha=0.05)


## ------------------------------------------------------------------------
TOSTone.raw(m=-0.3,mu=0.0,sd=1,n=100,low_eqbound=-0.5, high_eqbound=0.5, alpha=0.05)

## ------------------------------------------------------------------------
TOSTone(m=-0.3,mu=0.0,sd=1,n=100,low_eqbound_d=-0.5, high_eqbound_d=0.5, alpha=0.05)

## ------------------------------------------------------------------------
TOSTone.raw(m=-1,mu=0.0,sd=1.2,n=100,low_eqbound=-0.8, high_eqbound=0.8, alpha=0.05)

## ------------------------------------------------------------------------
TOSTr(n=50, r = 0.1, low_eqbound_r=-0.2, high_eqbound_r=0.2, alpha=0.05)


## ------------------------------------------------------------------------
#Running an original two t test procedure for equivalence
#######
equivint<-0.3
corxy<-0.02
n<-100
alpha<-0.05

zei <-log((1 + equivint)/(1 - equivint))/2
zcorxy <-log((1 + corxy)/(1-corxy))/2
equivt1_fz <-(zcorxy - zei)/(1/sqrt(n - 3))
pvalue1_fz <-pnorm(equivt1_fz)
equivt2_fz <-(zcorxy + zei)/(1/sqrt(n - 3))
pvalue2_fz <-1 - pnorm(equivt2_fz)
ifelse (pvalue1_fz <= alpha & pvalue2_fz <= alpha, decis_fz <-"The null hypothesis that the correlation between var1 and var2 falls outside of the equivalence interval can be rejected.", decis_fz <-"The null hypothesis that the correlation between var1 and var2 falls outside of the equivalence interval cannot be rejected.")
pvalue1_fz
pvalue2_fz

## ------------------------------------------------------------------------
TOSTr(n=100, r = 0.02, low_eqbound_r=-0.3, high_eqbound_r=0.3, alpha=0.05)


## ------------------------------------------------------------------------
TOSTmeta(ES=0.12, var=0.0081, se=0.09, low_eqbound_d=-0.2, high_eqbound_d=0.2, alpha=0.05)

## ------------------------------------------------------------------------
TOSTtwo.prop(prop1 =  133/262,
             prop2 =  136/265,
             n1 = 262,
             n2 = 265, 
             low_eqbound = -0.12,
             high_eqbound = 0.12, 
             alpha = 0.025, 
             plot = TRUE)

## ------------------------------------------------------------------------
TOSTtwo.prop(prop1 =  29/148,
             prop2 =  30/138,
             n1 = 148,
             n2 = 138, 
             low_eqbound = -0.15,
             high_eqbound = 0.15, 
             alpha = 0.05, 
             plot = TRUE)

## ------------------------------------------------------------------------
TOSTtwo.prop(prop1 =  18/246,
             prop2 =  15/224,
             n1 = 246,
             n2 = 224, 
             low_eqbound = -0.10,
             high_eqbound = 0.10, 
             alpha = 0.05, 
             plot = TRUE)


## ------------------------------------------------------------------------
TOSTtwo.prop(prop1 =  187/583,
             prop2 =  95/328,
             n1 = 583,
             n2 = 328, 
             low_eqbound = -0.10,
             high_eqbound = 0.10, 
             alpha = 0.05, 
             plot = TRUE)


## ------------------------------------------------------------------------
TOSTtwo.prop(prop1 =  331/583,
             prop2 =  201/328,
             n1 = 583,
             n2 = 328, 
             low_eqbound = -0.10,
             high_eqbound = 0.10, 
             alpha = 0.05, 
             plot = TRUE)


