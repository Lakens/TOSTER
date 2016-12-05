#' TOST function for a dependent t-test (Cohen's dz)
#' @param m1 mean of group 1
#' @param m2 mean of group 2
#' @param sd1 standard deviation of group 1
#' @param sd2 standard deviation of group 2
#' @param n sample size (pairs)
#' @param r12 correlation of dependent variable between group 1 and  group 2
#' @param low_eqbound_dz lower equivalence bounds (e.g., -0.5) expressed in standardized mean difference (Cohen's dz)
#' @param high_eqbound_dz upper equivalence bounds (e.g., 0.5) expressed in standardized mean difference (Cohen's dz)
#' @param alpha alpha level (default = 0.05)
#' @return Returns TOST t-value 1, TOST p-value 1, TOST t-value 2, TOST p-value 2, degrees of freedom, low equivalence bound, high equivalence bound, low equivalence bound in dz, high equivalence bound in dz, Lower limit confidence interval TOST, Upper limit confidence interval TOST
#' @examples
#' TOSTpaired(n=65,m1=5.83,m2=5.75,sd1=1.17,sd2=1.29,r12=0.75,low_eqbound_dz=-0.4,high_eqbound_dz=0.4)
#' @section References:
#' Mara, C. A., & Cribbie, R. A. (2012). Paired-Samples Tests of Equivalence. Communications in Statistics - Simulation and Computation, 41(10), 1928-1943. https://doi.org/10.1080/03610918.2011.626545, formula page 1932. Note there is a typo in the formula: n-1 should be n (personal communication, 31-8-2016)
#' @importFrom stats pnorm pt qnorm qt
#' @export
#' 

TOSTpaired<-function(n,m1,m2,sd1,sd2,r12,low_eqbound_dz, high_eqbound_dz, alpha){
  if(missing(alpha)) {
    alpha <- 0.05
  }  
  sdif<-sqrt(sd1^2+sd2^2-2*r12*sd1*sd2)
  low_eqbound<-low_eqbound_dz*sdif
  high_eqbound<-high_eqbound_dz*sdif
  se<-sdif/sqrt(n)
  t<-(m1-m2)/se
  degree_f<-n-1
  pttest<-2*pt(abs(t), degree_f, lower.tail=FALSE)
  t1<-((m1-m2)+(low_eqbound_dz*sdif))/se
  p1<-1-pt(t1, degree_f, lower.tail=FALSE)
  t2<-((m1-m2)+(high_eqbound_dz*sdif))/se
  p2<-pt(t2, degree_f, lower.tail=FALSE)
  ttost<-ifelse(abs(t1) < abs(t2), t1, t2)
  LL90<-((m1-m2)-qt(1-alpha, degree_f)*se)
  UL90<-((m1-m2)+qt(1-alpha, degree_f)*se)
  ptost<-max(p1,p2)
  results<-data.frame(t1,p1,t2,p2,degree_f,LL90,UL90)
  colnames(results) <- c("t-value 1","p-value 1","t-value 2","p-value 2","df", paste("Lower Limit ",100*(1-alpha*2),"% CI",sep=""),paste("Upper Limit ",100*(1-alpha*2),"% CI",sep=""))
  dif<-(m1-m2)
  LL95<-((m1-m2)-qt(1-(alpha/2), degree_f)*se)
  UL95<-((m1-m2)+qt(1-(alpha/2), degree_f)*se)
  testoutcome<-ifelse(pttest<alpha,"significant","non-significant")
  TOSToutcome<-ifelse(ptost<alpha,"significant","non-significant")
  plot(NA, ylim=c(0,1), xlim=c(min(LL90,low_eqbound)-max(UL90-LL90, high_eqbound-low_eqbound)/10, max(UL90,high_eqbound)+max(UL90-LL90, high_eqbound-low_eqbound)/10), bty="l", yaxt="n", ylab="",xlab="Mean Difference")
  points(x=dif, y=0.5, pch=15, cex=2)
  abline(v=high_eqbound, lty=2)
  abline(v=low_eqbound, lty=2)
  abline(v=0, lty=2, col="grey")
  segments(LL90,0.5,UL90,0.5, lwd=3)
  segments(LL95,0.5,UL95,0.5, lwd=1)
  title(main=paste("Equivalence bounds ",round(low_eqbound,digits=3)," and ",round(high_eqbound,digits=3),"\nMean difference = ",round(dif,digits=3)," \n TOST: ", 100*(1-alpha*2),"% CI [",round(LL90,digits=3),";",round(UL90,digits=3),"] ", TOSToutcome," \n NHST: ", 100*(1-alpha),"% CI [",round(LL95,digits=3),";",round(UL95,digits=3),"] ", testoutcome,sep=""), cex.main=1)
  message(cat("Using alpha = ",alpha," the NHST t-test was ",testoutcome,", t(",degree_f,") = ",t,", p = ",pttest,sep=""))
  cat("\n")
  message(cat("Using alpha = ",alpha," the equivalence test was ",TOSToutcome,", t(",degree_f,") = ",ttost,", p = ",ptost,sep=""))
  TOSTresults<-data.frame(t1,p1,t2,p2,degree_f)
  colnames(TOSTresults) <- c("t-value 1","p-value 1","t-value 2","p-value 2","df")
  bound_dz_results<-data.frame(low_eqbound_dz,high_eqbound_dz)
  colnames(bound_dz_results) <- c("low bound dz","high bound dz")
  bound_results<-data.frame(low_eqbound,high_eqbound)
  colnames(bound_results) <- c("low bound raw","high bound raw")
  CIresults<-data.frame(LL90,UL90)
  colnames(CIresults) <- c(paste("Lower Limit ",100*(1-alpha*2),"% CI raw",sep=""),paste("Upper Limit ",100*(1-alpha*2),"% CI raw",sep=""))
  cat("TOST results:\n")
  print(TOSTresults)  
  cat("\n")
  cat("Equivalence bounds (Cohen's dz):\n")
  print(bound_dz_results)  
  cat("\n")
  cat("Equivalence bounds (raw scores):\n")
  print(bound_results)  
  cat("\n")
  cat("TOST confidence interval:\n")
  print(CIresults)
  invisible(list(TOST_t1=t1,TOST_p1=p1,TOST_t2=t2,TOST_p2=p2, TOST_df=degree_f,alpha=alpha,low_eqbound=low_eqbound,high_eqbound=high_eqbound,low_eqbound_dz=low_eqbound_dz,high_eqbound_dz=high_eqbound_dz, LL_CI_TOST=LL90,UL_CI_TOST=UL90))
}