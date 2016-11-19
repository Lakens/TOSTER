#' TOST function for a dependent t-test using equivalence bounds based on Cohen's dz
#' @param m1 mean of group 1
#' @param m2 mean of group 2
#' @param sd1 standard deviation of group 1
#' @param sd2 standard deviation of group 2
#' @param n sample size (pairs)
#' @param r12 correlation of dependent variable between group 1 and  group 2
#' @param low_eqbound_dz lower equivalence bounds (e.g., -0.5) expressed in standardized mean difference (Cohen's dz)
#' @param high_eqbound_dz upper equivalence bounds (e.g., 0.5) expressed in standardized mean difference (Cohen's dz)
#' @param alpha alpha level (default = 0.05)
#' @return Returns dataframe with t-value, degrees of freedom, p-value, and lower and upper limit of confidence interval
#' @examples
#' TOSTdz(n=65,m1=5.830769,m2=5.7538462,sd1=1.1668498,sd2=1.2994081,r12=0.744992,low_eqbound_dz=-0.337341706,high_eqbound_dz=0.337341706, alpha=0.05)
#' @export
#' 

TOSTdz<-function(n,m1,m2,sd1,sd2,r12,low_eqbound_dz, high_eqbound_dz, alpha){
  if(missing(alpha)) {
    alpha <- 0.05
  }  
  sdif<-sqrt(sd1^2+sd2^2-2*r12*sd1*sd2)
  low_eqbound<-low_eqbound_dz*sdif
  high_eqbound<-high_eqbound_dz*sdif
  se<-sdif/sqrt(n)
  t<-(m1-m2)/se
  degree_f<-n-1
  pttest<-2*pt(abs(t), degree_f, lower=FALSE)
  t1<-((m1-m2)+(low_eqbound_dz*sdif))/se
  p1<-1-pt(t1, degree_f, lower=FALSE)
  t2<-((m1-m2)+(high_eqbound_dz*sdif))/se
  p2<-pt(t2, degree_f, lower=FALSE)
  ttost<-ifelse(abs(t1) < abs(t2), t1, t2)
  LL90<-((m1-m2)-qt(1-alpha, degree_f)*se)
  UL90<-((m1-m2)+qt(1-alpha, degree_f)*se)
  ptost<-max(p1,p2)
  results<-data.frame(ttost,degree_f,ptost,LL90,UL90)
  dif<-(m1-m2)
  LL95<-((m1-m2)-qt(1-(alpha/2), degree_f)*se)
  UL95<-((m1-m2)+qt(1-(alpha/2), degree_f)*se)
  testoutcome<-ifelse(pttest<0.05,"significant","non-significant")
  TOSToutcome<-ifelse(ptost<0.05,"significant","non-significant")
  plot(NA, ylim=c(0,1), xlim=c(min(LL90,low_eqbound)-max(UL90-LL90, high_eqbound-low_eqbound)/10, max(UL90,high_eqbound)+max(UL90-LL90, high_eqbound-low_eqbound)/10), bty="l", yaxt="n", ylab="",xlab="Mean Difference")
  points(x=dif, y=0.5, pch=15, cex=2)
  abline(v=high_eqbound, lty=2)
  abline(v=low_eqbound, lty=2)
  abline(v=0, lty=2, col="grey")
  segments(LL90,0.5,UL90,0.5, lwd=3)
  segments(LL95,0.5,UL95,0.5, lwd=1)
  title(main=paste("Equivalence bounds ",round(low_eqbound,digits=3)," and ",round(high_eqbound,digits=3),"\nMean difference = ",round(dif,digits=3)," \n TOST: 90% CI [",round(LL90,digits=3),";",round(UL90,digits=3),"] ", TOSToutcome," \n NHST: 95% CI [",round(LL95,digits=3),";",round(UL95,digits=3),"] ", testoutcome,sep=""), cex.main=1)
  message(cat("The NHST t-test was ",testoutcome,", t(",degree_f,") = ",t,", p = ",pttest,sep=""))
  message(cat("The equivalence test was ",TOSToutcome,", t(",degree_f,") = ",ttost,", p = ",ptost,sep=""))
  return(results)
}