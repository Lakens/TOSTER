#' TOST function for a correlations
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
  zLL90<-(log((1+r)/(1-r))/2)-qnorm(1-alpha)*sqrt(1/(n-3))
  zUL90<-(log((1+r)/(1-r))/2)+qnorm(1-alpha)*sqrt(1/(n-3))
  LL90<-(exp(1)^(2*zLL90)-1)/(exp(1)^(2*zLL90)+1)
  UL90<-(exp(1)^(2*zUL90)-1)/(exp(1)^(2*zUL90)+1)
  zLL95<-(log((1+r)/(1-r))/2)-qnorm(1-(alpha/2))*sqrt(1/(n-3))
  zUL95<-(log((1+r)/(1-r))/2)+qnorm(1-(alpha/2))*sqrt(1/(n-3))
  LL95<-(exp(1)^(2*zLL95)-1)/(exp(1)^(2*zLL95)+1)
  UL95<-(exp(1)^(2*zUL95)-1)/(exp(1)^(2*zUL95)+1)
  results<-data.frame(ptost,LL90,UL90)
  testoutcome<-ifelse(pttest<0.05,"significant","non-significant")
  TOSToutcome<-ifelse(ptost<0.05,"significant","non-significant")
  plot(NA, ylim=c(0,1), xlim=c(min(LL90,low_eqbound_r)-max(UL90-LL90, high_eqbound_r-low_eqbound_r)/10, max(UL90,high_eqbound_r)+max(UL90-LL90, high_eqbound_r-low_eqbound_r)/10), bty="l", yaxt="n", ylab="",xlab="Correlation")
  points(x=r, y=0.5, pch=15, cex=2)
  abline(v=high_eqbound_r, lty=2)
  abline(v=low_eqbound_r, lty=2)
  abline(v=0, lty=2, col="grey")
  segments(LL90,0.5,UL90,0.5, lwd=3)
  segments(LL95,0.5,UL95,0.5, lwd=1)
  title(main=paste("Equivalence bounds ",round(low_eqbound_r,digits=3)," and ",round(high_eqbound_r,digits=3),"\nr = ",round(r,digits=3)," \n TOST: 90% CI [",round(LL90,digits=3),";",round(UL90,digits=3),"] ", TOSToutcome," \n NHST: 95% CI [",round(LL95,digits=3),";",round(UL95,digits=3),"] ", testoutcome,sep=""), cex.main=1)
  message(cat("The NHST t-test was ",testoutcome,", p = ",pttest,sep=""))
  message(cat("The equivalence test was ",TOSToutcome,", p = ",ptost,sep=""))
  return(results)
}