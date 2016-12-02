#' TOST function for a correlations
#' @param n number of pairs of observations
#' @param r observed correlation
#' @param low_eqbound_r lower equivalence bounds (e.g., -0.3) expressed in a correlation effect size 
#' @param high_eqbound_r upper equivalence bounds (e.g., 0.3) expressed in a correlation effect size 
#' @param alpha alpha level (default = 0.05)
#' @return Returns TOST p-value 1, TOST p-value 2, alpha, low equivalence bound r, high equivalence bound r, Lower limit confidence interval TOST, Upper limit confidence interval TOST
#' @examples
#' TOSTr(n=100, r = 0.02, low_eqbound_r=-0.3, high_eqbound_r=0.3, alpha=0.05)
#' @section References:
#' Goertzen, J. R., & Cribbie, R. A. (2010). Detecting a lack of association: An equivalence testing approach. British Journal of Mathematical and Statistical Psychology, 63(3), 527-537. https://doi.org/10.1348/000711009X475853, formula page 531.
#' @export
#' 

TOSTr<-function(n, r, low_eqbound_r, high_eqbound_r, alpha){
  if(missing(alpha)) {
    alpha <- 0.05
  }  
  z1<-((log((1+abs(r))/(1-abs(r)))/2)-(log((1+low_eqbound_r)/(1-low_eqbound_r))/2))/(sqrt(1/(n-3)))
  z2<-((log((1+abs(r))/(1-abs(r)))/2)-(log((1+high_eqbound_r)/(1-high_eqbound_r))/2))/(sqrt(1/(n-3)))
  p1<-ifelse(low_eqbound_r<r,pnorm(-abs(z1)),1-pnorm(-abs(z1)))
  p2<-ifelse(high_eqbound_r>r,pnorm(-abs(z2)),1-pnorm(-abs(z2)))
  ptost<-max(p1,p2)
  pttest<-2*(1-pt(abs(r)*sqrt(n-2)/sqrt(1-abs(r)^2),n-2))
  zLL90<-(log((1+r)/(1-r))/2)-qnorm(1-alpha)*sqrt(1/(n-3))
  zUL90<-(log((1+r)/(1-r))/2)+qnorm(1-alpha)*sqrt(1/(n-3))
  LL90<-(exp(1)^(2*zLL90)-1)/(exp(1)^(2*zLL90)+1)
  UL90<-(exp(1)^(2*zUL90)-1)/(exp(1)^(2*zUL90)+1)
  zLL95<-(log((1+r)/(1-r))/2)-qnorm(1-(alpha/2))*sqrt(1/(n-3))
  zUL95<-(log((1+r)/(1-r))/2)+qnorm(1-(alpha/2))*sqrt(1/(n-3))
  LL95<-(exp(1)^(2*zLL95)-1)/(exp(1)^(2*zLL95)+1)
  UL95<-(exp(1)^(2*zUL95)-1)/(exp(1)^(2*zUL95)+1)
  testoutcome<-ifelse(pttest<alpha,"significant","non-significant")
  TOSToutcome<-ifelse(ptost<alpha,"significant","non-significant")
  plot(NA, ylim=c(0,1), xlim=c(min(LL90,low_eqbound_r)-max(UL90-LL90, high_eqbound_r-low_eqbound_r)/10, max(UL90,high_eqbound_r)+max(UL90-LL90, high_eqbound_r-low_eqbound_r)/10), bty="l", yaxt="n", ylab="",xlab="Correlation")
  points(x=r, y=0.5, pch=15, cex=2)
  abline(v=high_eqbound_r, lty=2)
  abline(v=low_eqbound_r, lty=2)
  abline(v=0, lty=2, col="grey")
  segments(LL90,0.5,UL90,0.5, lwd=3)
  segments(LL95,0.5,UL95,0.5, lwd=1)
  title(main=paste("Equivalence bounds ",round(low_eqbound_r,digits=3)," and ",round(high_eqbound_r,digits=3),"\nr = ",round(r,digits=3)," \n TOST: 90% CI [",round(LL90,digits=3),";",round(UL90,digits=3),"] ", TOSToutcome," \n NHST: 95% CI [",round(LL95,digits=3),";",round(UL95,digits=3),"] ", testoutcome,sep=""), cex.main=1)
  message(cat("Using alpha = ",alpha," the NHST t-test was ",testoutcome,", p = ",pttest,sep=""))
  message(cat("Using alpha = ",alpha, " the equivalence test was ",TOSToutcome,", p = ",ptost,sep=""))
  TOSTresults<-data.frame(p1,p2)
  colnames(TOSTresults) <- c("p-value 1","p-value 2")
  bound_r_results<-data.frame(low_eqbound_r,high_eqbound_r)
  colnames(bound_r_results) <- c("low bound r","high bound r")
  CIresults<-data.frame(LL90,UL90)
  colnames(CIresults) <- c(paste("Lower Limit ",100*(1-alpha*2),"% CI raw",sep=""),paste("Upper Limit ",100*(1-alpha*2),"% CI raw",sep=""))
  cat("TOST results:\n")
  print(TOSTresults)  
  cat("\n")
  cat("Equivalence bounds (r):\n")
  print(bound_r_results)  
  cat("\n")
  cat("TOST confidence interval:\n")
  print(CIresults)
  invisible(list(TOST_p1=p1,TOST_p2=p2,alpha=alpha,low_eqbound_r=low_eqbound_r,high_eqbound_r=high_eqbound_r, LL_CI_TOST=LL90,UL_CI_TOST=UL90))
}