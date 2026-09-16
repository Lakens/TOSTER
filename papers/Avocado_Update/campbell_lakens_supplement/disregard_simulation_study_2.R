set.seed(123)


# Functions for rounding:
ceiling_dec <- function(x, level=1) round((x + 5*10^(-level-1))/0.5, level)*0.5

# number of simulations:
nSim<-10000

resultsmat <- data.frame(
JJJindex=NA,
trueP2_num= NA, 
power=NA, 
as.matrix(expand.grid(Nvar=c(2,4),sigma2=c(0.4,0.5,1,1.0001), Nsample=c(60, 180, 450, 1000))))

resultsmat$trueP2 <- NA

resultsmat$alpha_sig <- 0.05

dim(resultsmat)[1]


resultsmat<-resultsmat[order(resultsmat$Nvar, resultsmat$sigma2),]
pval_list<-list()
problist<-list()

resultsmat<-resultsmat[sample(1:dim(resultsmat)[1]),]

for(jjj in 1:dim(resultsmat)[1]){

print(jjj)

resultsmat[jjj, "JJJindex"] <- paste("jjj",jjj,sep="")
sigma2 	<- resultsmat[jjj, "sigma2"] 
N 		<- resultsmat[jjj, "Nsample"]
nVar 	<- resultsmat[jjj, "Nvar"]

basematrix <- data.frame(expand.grid(
				X1=c(0,1),
				X2=c(0,1),
				X3=c(0,1),
				X4=c(0,1)))

X1 <- rep(basematrix$X1,1000000)
X2 <- rep(basematrix$X2,1000000)
X3 <- rep(basematrix$X3,1000000)
X4 <- rep(basematrix$X4,1000000)

epsilon <- rnorm(length(X1), 0, sqrt(sigma2))


if(nVar==2){
X  <- as.matrix(cbind(1, X1, X2))
betavec <- c(0.0, 0.2, 0.3)
if(sigma2==1.0001){betavec <- c(0.0, 0.0, 0.0)}
Y <- c(X%*%betavec) + epsilon
sigmaXY <-(c(cov(X[,-1], Y)))
SIGMAX <- cov(X[,-1])
(trueP2_num <-(t(sigmaXY)%*%solve(SIGMAX)%*%sigmaXY))

resultsmat[jjj, "trueP2_num"] <- trueP2_num
resultsmat[jjj, "trueP2"] <-(resultsmat[jjj, "trueP2_num"]/var(Y))

confirmmod<-lm(Y~ X[,-1])
print(c(summary(confirmmod)$r.squared,resultsmat[jjj, "trueP2"]))

resultsmat[jjj, "trueP2_v2"] <- summary(confirmmod)$r.squared

Xmatrix_ <- X[1:N,]
}

if(nVar==4){
X  <- as.matrix(cbind(1, X1, X2, X3, X4))
betavec <- c(0.0, 0.2, 0.2, -0.1, -0.2)
if(sigma2==1.0001){betavec <- c(0.0, 0.0, 0.0, 0.0, 0.0)}
Y <- c(X%*%betavec) + epsilon
sigmaXY <- c(cov(X[,-1], Y))
SIGMAX <- cov(X[,-1])
(trueP2_num <- (t(sigmaXY)%*%solve(SIGMAX)%*%sigmaXY))

resultsmat[jjj, "trueP2_num"] <- trueP2_num
resultsmat[jjj, "trueP2"] <-(resultsmat[jjj, "trueP2_num"]/var(Y))


confirmmod<-lm(Y~ X[,-1])
print(c(summary(confirmmod)$r.squared,resultsmat[jjj, "trueP2"]))

resultsmat[jjj, "trueP2_v2"] <- summary(confirmmod)$r.squared

Xmatrix_ <- X[1:N,]
}


pval_list[[jjj]]<-list()
Deltavec<-seq(0.01,0.10,0.005)

for(iii in 1:nSim){

	epsilon <- rnorm(N, 0, sqrt(sigma2))
	Y <- Xmatrix_%*%betavec + epsilon

	lmmodel<-lm(Y~ Xmatrix_[,-1])
	Xmatrix <- model.matrix(lmmodel)

	R2 <- summary(lmmodel)$r.squared
	Fstat <- summary(lmmodel)$fstatistic[1]
	K <- dim(Xmatrix)[2] - 1
	N <- dim(Xmatrix)[1]

	pval_list[[jjj]][[iii]]<-vector() 

	for(kk in 1:length(Deltavec)){
	Delta<-Deltavec[kk]
	pval <- pf(Fstat,df1=K,df2=N-K-1,ncp=(N*Delta)/(1-Delta),lower.tail=TRUE)
	pval_list[[jjj]][[iii]][kk]<-pval
	}
	
	}

pvalmat<-as.data.frame(pval_list[[jjj]], col.names=NA, row.names=Deltavec)	
problist[[jjj]]<-rowMeans(pvalmat<resultsmat[jjj, "alpha_sig"])


#resultsmat[jjj,] <- mean(pval_vec< resultsmat[jjj, "alpha_sig"])
print(resultsmat[jjj,])
}

problist_df<-as.data.frame(problist, col.names=c(1:dim(resultsmat)[1]))
colnames(problist_df)<-paste("jjj",c(1:dim(resultsmat)[1]), sep="")
problist_df$Delta<-Deltavec
library(tidyr)
problist_long<-gather(problist_df,JJJindex,pr_less_alpha,jjj1:jjj32, factor_key=TRUE)


resultsmatall<-merge(resultsmat,problist_long, by="JJJindex",all=TRUE)

resultsmatall$trueP2_v2[resultsmatall$trueP2_v2>0.03 & resultsmatall$trueP2_v2<0.05]<-round( mean(resultsmatall$trueP2_v2[resultsmatall$trueP2_v2>0.03 & resultsmatall$trueP2_v2<0.05]),3)

resultsmatall$trueP2_v2[resultsmatall$trueP2_v2>0.06 & resultsmatall$trueP2_v2<0.07]<-round( mean(resultsmatall$trueP2_v2[resultsmatall$trueP2_v2>0.06 & resultsmatall$trueP2_v2<0.07]),3)

resultsmatall$trueP2_v2[resultsmatall$trueP2_v2>0.07]<-round(  mean(resultsmatall$trueP2_v2[resultsmatall$trueP2_v2>0.07]),3)

resultsmatall$trueP2_v2[resultsmatall$trueP2_v2<0.02]<-round(  mean(resultsmatall$trueP2_v2[resultsmatall$trueP2_v2<0.02]),3)

resultsmatall$P2 <-as.factor(round((resultsmatall$trueP2_v2),3))
resultsmatall$N <-as.factor(resultsmatall$Nsample)

resultsmatall$g <- as.factor(paste(as.character(resultsmatall$P2),as.character(resultsmatall$N), sep="_"))

resultsmatall <- transform(resultsmatall,
  Nvar = factor(Nvar, levels = sort(unique(resultsmatall$Nvar)), c( (("K = 2")),(("K = 4")))  ))

#saveRDS(resultsmatall,'~/Desktop/UBC/ThesisProposal/Rcode/lakens_project/resultsmatall.rds')

resultsmatall<- readRDS('~/Desktop/UBC/ThesisProposal/Rcode/lakens_project/resultsmatall.rds')

resultsmatall$KK <- as.numeric(resultsmatall$Nvar)*2
resultsmatall$NN <- as.numeric(as.character(resultsmatall$N))


resultsmatall$power <-NA
power_estimate<-function(Delta, K, N){
Fstat_star <- qf(0.05,df1=K,df2=N-K-1,ncp=(N*Delta)/(1-Delta),lower.tail=TRUE)
power<-pf(Fstat_star,df1=K,df2=N-K-1,lower.tail=TRUE)
return(round(power,3))}


resultsmatall[resultsmatall$P2==0,]$power <- apply(resultsmatall[resultsmatall$P2==0, c("Delta","KK","NN")],1, function(y){ power_estimate(y[1],y[2],y[3])})


####  plot with truncated axis:
library(ggplot2)
qplot(x=Delta, y= pr_less_alpha, group= g, pch=N, lty=N, col= P2, data= resultsmatall)+geom_line()+  geom_hline(yintercept = 0.05)+ facet_grid(Nvar ~ . ) + scale_x_continuous(breaks = seq(0, 0.1, by = 0.01))  + scale_y_continuous(breaks = seq(-0.01, 1, by = 0.05))+  coord_cartesian(ylim = c(0, 0.20))+ labs(x = expression(Delta), y = expression("probability of p" < alpha))+geom_line(aes(y = power, x= Delta), color = "black", linetype = "dotted")

####  plot with full axis:
library(ggplot2)
qplot(x=Delta, y= pr_less_alpha, group= g, pch=N, lty=N, col= P2, data= resultsmatall)+geom_line()+  geom_hline(yintercept = 0.05)+ facet_grid(Nvar ~. ) + scale_x_continuous(breaks = seq(0, 0.1, by = 0.01))  + scale_y_continuous(breaks = seq(0, 1, by = 0.1))+ labs(x = expression(Delta),y = expression("probability of p" < alpha))+ coord_cartesian(ylim = c(0, 1)) +geom_line(aes(y = power, x= Delta), color = "black", linetype = "dotted")




