##  
# Here are the settings for the plot:
BFthres<-3
alpha1<-0.05
alpha2<-0.05
p <- K <- 12


## code to calculate/estimate intersection of two curves:
curve_intersect <- function(curve1, curve2, empirical=TRUE, domain=NULL) {
  if (!empirical & missing(domain)) {
    stop("'domain' must be provided with non-empirical curves")
  }

  if (!empirical & (length(domain) != 2 | !is.numeric(domain))) {
    stop("'domain' must be a two-value numeric vector, like c(0, 10)")
  }

  if (empirical) {
    # Approximate the functional form of both curves
    curve1_f <- approxfun(curve1$x, curve1$y, rule = 2)
    curve2_f <- approxfun(curve2$x, curve2$y, rule = 2)

    # Calculate the intersection of curve 1 and curve 2 along the x-axis
    point_x <- uniroot(function(x) curve1_f(x) - curve2_f(x),
                       c(min(curve1$x), max(curve1$x)))$root

    # Find where point_x is in curve 2
    point_y <- curve2_f(point_x)
  } else {
    # Calculate the intersection of curve 1 and curve 2 along the x-axis
    # within the given domain
    point_x <- uniroot(function(x) curve1(x) - curve2(x), domain)$root

    # Find where point_x is in curve 2
    point_y <- curve2(point_x)
  }

  return(list(x = point_x, y = point_y))
}

####

library(BayesFactor)

# Set colours #
green<-rgb(0.1, 0.8, 0.1, alpha=0.25)
blue<-rgb(0.1, 0.8, 0.2, alpha=0.25)
yellow<-rgb(0.3, 0.8, 0.2, alpha=0.25)
clear<-rgb(1, 0, 0, alpha=0.0001)
red<-rgb(1, 0, 0.1, alpha=0.25)


#################################
lseq<-function(from, to, length.out){10^(seq(log10(from), log10(to), length.out = length.out))}






resultmat<-data.frame()
Nvec<-c(lseq(30,200,50),200,240:260, 300, 500,750,1000)
BFvec<-c(round(1/BFthres,3),BFthres)

k<-1
for(i in 1:length(BFvec)){
	for(j in 1:length(Nvec)){
N <- Nvec[j]
Delta <- 0.1
BF <- BFvec[i]

matchBF<-function(x){(BF-linearReg.R2stat(N, p, x, rscale = "medium", simple = TRUE))}

R2 <- tryCatch(uniroot(matchBF, interval=c(0.0001,0.999))$root, error=function(x) 0)

resultmat[k,"BF"]<-BF
resultmat[k,"N"]<-N
resultmat[k,"R2"]<-R2
resultmat[k,"p"]<-p
resultmat[k,"Delta"]<-Delta

R2_CET<-function(x){
Fstat <- (x/K)/( (1-x)/(N-K-1) )
pval1 <- pf(Fstat,df1 = K, df2 = N-K-1, lower.tail=FALSE)
pval2 <- pf(Fstat, df1 = K, df2 = N-K-1, ncp = (N*Delta)/(1-Delta), lower.tail=TRUE ) 
if(pval1<alpha1){return(list(pval=pval1,result=1))}
if(pval1>=alpha1){return(list(pval=pval2,result=-1)) }
	}

R2_NHST<-function(x){
	Fstat <- (x/K)/( (1-x)/(N-K-1) )
	 pf(Fstat,df1 = K, df2 = N-K-1,lower.tail=FALSE)
}

R2_EQUIV<-function(x){
	Fstat <- (x/K)/( (1-x)/(N-K-1) )
	pf(Fstat, df1 = K, df2 = N-K-1, ncp = (N*Delta)/(1-Delta) ,lower.tail=TRUE) 
}


if(BF>=1){
R2_freq <- uniroot(function(x){R2_NHST(x)-alpha1}, interval=c(0,1))$root}


if(BF<1){
R2_freq <- uniroot(function(x){R2_EQUIV(x)-alpha2}, interval=c(0,1))$root}



resultmat[k, "R2_freq"] <- R2_freq

k <- k+1
print(k)
}}

library(ggplot2)



######################
x <-data.frame(N=resultmat[resultmat$BF==unique(resultmat $BF)[1],]$N)
x$R2<-resultmat[resultmat $BF==unique(resultmat $BF)[1],]$R2_freq
x$x2<-resultmat[resultmat $BF==unique(resultmat $BF)[2],]$R2_freq

interxect_x<-round(curve_intersect(data.frame(x=x$N, y=x$R2), data.frame(x=x$N, y=x$x2) )$x)


equivplot<- ggplot(x, aes(x=N, y=R2))      +
    geom_ribbon(data=subset(x, 0 <= N & N <= interxect_x), 
          aes(ymin=R2,ymax=x2), fill="yellow", alpha="0.5") +
 geom_ribbon(data=subset(x, 0 <= (N+1) & (N-1) <= interxect_x+1), 
          aes(ymin=0,ymax= R2), fill="red", alpha="0.5") +
 geom_ribbon(data=subset(x, interxect_x <= (N) & (N-1) <=10000), 
          aes(ymin=0,ymax= x2), fill="red", alpha="0.5")  +
 geom_ribbon(data=subset(x, interxect_x <= N & N <=10000), 
          aes(ymin= x2,ymax= R2), fill="palegreen1", alpha="0.5")+
 geom_ribbon(data=subset(x, interxect_x <= (N+1) & (N-3) <=10000), 
          aes(ymin= R2,ymax= 1), fill="green", alpha="0.5")+
 geom_ribbon(data=subset(x,  0<= (N+2) & (N-2) <=interxect_x), 
          aes(ymin= x2,ymax= 1), fill="green", alpha="0.5")+ 
    geom_line(aes(y = R2)) + 
    geom_line(aes(y = x2))+    geom_line(aes(y = 0.10), lty=3)+ 
    scale_y_log10(breaks=c(.001,.1,0.5,1), limits=c(0.0001,1))+ scale_x_log10(breaks=c(10,30,50,100, interxect_x ,500,1000,5000,10000))
                     
       
       
 




#####



x <-data.frame(N=resultmat[resultmat$BF==unique(resultmat $BF)[1],]$N)
x$R2<-resultmat[resultmat$BF==unique(resultmat $BF)[1],]$R2
x$x2<-resultmat[resultmat$BF==unique(resultmat $BF)[2],]$R2

BFplot<- ggplot(x, aes(x=N, y=R2))      +    geom_line(aes(y = R2)) + 
    geom_line(aes(y = x2)) + 
    scale_y_log10(breaks=c(.001,.1,0.5,1), limits=c(0.0001,1))+ scale_x_log10(breaks=c(10,30,50,100,500,1000,5000,10000)) +
    geom_ribbon(data=subset(x, 0 <= (N+1) & (N-1) <= 10000), 
          aes(ymin=x2,ymax=1), fill="green", alpha="0.5") +
 geom_ribbon(data=subset(x, 0 <= N & N <= 10000), 
          aes(ymin=R2,ymax= x2), fill="yellow", alpha="0.5") +
 geom_ribbon(data=subset(x, 0 <= (N+1) & (N-1) <=10000), 
          aes(ymin=0,ymax= R2), fill="red", alpha="0.5") 
                     
       
       
       require(gridExtra)

       grid.arrange(equivplot, BFplot, ncol=2)

		 
## This behaviour may seem odd but:
#plot(apply(cbind(lseq(50,500,90)), 1, function(x){ci.R2(R2=0.12,N=x, K=12, alpha.upper=0.05, alpha.lower=0.0, Random.Regressors=FALSE)$Upper.Conf.Limit.R2})~lseq(50,500,90), ylim=c(0,0.20)); abline(0.12,0)