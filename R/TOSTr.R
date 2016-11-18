#' TOST function
#' @param n number of pairs of observations
#' @param r observed correlation
#' @param low_eqbound_r lower equivalence bounds (e.g., -0.3) expressed in a correlation effect size 
#' @param high_eqbound_r upper equivalence bounds (e.g., 0.3) expressed in a correlation effect size 
#' @param alpha alpha level (default = 0.05)
#' @return Returns dataframe with p-value, and lower and upper limit of confidence interval
#' @examples
#' TOSTr(n=100, r = 0.02, low_eqbound_r=-0.3, high_eqbound_r=0.3, alpha=0.05)
#' @export
#' 

TOSTr<-function(n, r, low_eqbound_r, high_eqbound_r, alpha){
  if(missing(alpha)) {
    alpha <- 0.05
  }  
  z1<-((log((1+abs(r))/(1-abs(r)))/2)+(log((1+low_eqbound_r)/(1-low_eqbound_r))/2))/(sqrt(1/(n-3)))
  z2<-((log((1+abs(r))/(1-abs(r)))/2)+(log((1+high_eqbound_r)/(1-high_eqbound_r))/2))/(sqrt(1/(n-3)))
  ptost1<-ifelse(low_eqbound_r<r,pnorm(-abs(z1)),1-pnorm(-abs(z1)))
  ptost2<-ifelse(high_eqbound_r>r,pnorm(-abs(z2)),1-pnorm(-abs(z2)))
  ptost<-max(ptost1,ptost2)
  pttest<-2*(1-pt(abs(r)*sqrt(n-2)/sqrt(1-abs(r)^2),n-2))
  UL<-(log((1+r)/(1-r))/2)+qnorm(1-alpha)*sqrt(1/(n-3))
  LL<-(log((1+r)/(1-r))/2)-qnorm(1-alpha)*sqrt(1/(n-3))
  r_LL_90<-(exp(1)^(2*LL)-1)/(exp(1)^(2*LL)+1)
  r_UL_90<-(exp(1)^(2*UL)-1)/(exp(1)^(2*UL)+1)
  results<-data.frame(ptost,r_LL_90,r_UL_90)
  plotdata = data.frame(labels=c(paste(100*(1-2*alpha),"% CI r", sep=""),"Equivalence Range"), mean=c(r,0), lower=c(r_LL_90,low_eqbound_r), upper = c(r_UL_90,high_eqbound_r))
  plot(NA, xlim=c(.5,2.5), ylim=c(min(r_LL_90,low_eqbound_r)-0.5, max(r_UL_90,high_eqbound_r)+0.5), bty="l", xaxt="n", xlab="",ylab="Correlation")
  points(plotdata$mean[1:2], pch=19)
  r_UL_95<-(exp(1)^(2*((log((1+r)/(1-r))/2)+qnorm(0.975)*sqrt(1/(n-3))))-1)/(exp(1)^(2*((log((1+r)/(1-r))/2)+qnorm(0.975)*sqrt(1/(n-3))))+1)
  r_LL_95<-(exp(1)^(2*((log((1+r)/(1-r))/2)-qnorm(0.975)*sqrt(1/(n-3))))-1)/(exp(1)^(2*((log((1+r)/(1-r))/2)-qnorm(0.975)*sqrt(1/(n-3))))+1)
  points(1,r_UL_95,pch=10)
  points(1,r_LL_95,pch=10)
  axis(1, 1:2, plotdata$labels)
  segments(1:2,plotdata$lower[1:2],1:2,plotdata$upper[1:2])
  segments(1:1,plotdata$upper[1:1],1:2,plotdata$upper[1:1],lty=3)
  segments(2,0,0,0,lty=2)
  segments(1:1,plotdata$lower[1:1],1:2,plotdata$lower[1:1],lty=3)
  text(2, min(r_LL_90,low_eqbound_r)-0.1, paste("P-value",round(ptost, digits=3)), cex = .8)
  text(1, min(r_LL_90,high_eqbound_r)-0.1, paste("P-value",round(pttest, digits=3)), cex = .8)
  text(1.5, r, paste("r = ",round(r, digits=3)), cex = .8)
  title(main=paste("r = ",round(r,digits=3),", 95% CI [",round(r_LL_95,digits=3),";",round(r_UL_95,digits=3),"]",", p = ",round(pttest,digits=3),sep=""))
  testoutcome<-ifelse(pttest<0.05,"significant","non-significant")
  TOSToutcome<-ifelse(ptost<0.05,"significant","non-significant")
  message(cat("The NHST t-test was ",testoutcome,", p = ",pttest,sep=""))
  message(cat("The equivalence test was ",TOSToutcome,", p = ",ptost,sep=""))
  return(results)
}