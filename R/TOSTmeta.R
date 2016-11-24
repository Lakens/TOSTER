#' TOST function for meta-analysis
#' @param ES meta-analytic effect size
#' @param var meta-analytic variance
#' @param se standard error
#' @param low_eqbound_d lower equivalence bounds (e.g., -0.5) expressed in standardized mean difference (Cohen's d)
#' @param high_eqbound_d upper equivalence bounds (e.g., 0.5) expressed in standardized mean difference (Cohen's d)
#' @param alpha alpha level (default = 0.05)
#' @return Returns dataframe with Z-value, p-value, and lower and upper limit of confidence interval
#' @examples
#' TOSTmeta(ES=0.12, var=0.0081, se=0.09, low_eqbound_d=-0.2, high_eqbound_d=0.2, alpha=0.05)
#' @export
#' 


TOSTmeta<-function(ES,var,se,low_eqbound_d, high_eqbound_d, alpha){
  if(missing(alpha)) {
    alpha<-0.05
  }
  Z1<-(ES+high_eqbound_d)/se
  p1<-pnorm(Z1, lower=FALSE) 
  Z2<-(ES+low_eqbound_d)/se
  p2<-pnorm(Z2, lower=TRUE) 
  Z<-(ES/se)
  pttest<-pnorm(-abs(Z))
  LL90<-ES-qnorm(1-alpha)*(se)
  UL90<-ES+qnorm(1-alpha)*(se)
  LL95<-ES-qnorm(1-alpha/2)*(se)
  UL95<-ES+qnorm(1-alpha/2)*(se)
  ptost<-max(p1,p2) #Get highest p-value for summary TOST result
  Ztost<-ifelse(abs(Z1) < abs(Z2), Z1, Z2) #Get lowest t-value for summary TOST result
  results<-data.frame(Ztost,ptost,LL90,UL90)
  testoutcome<-ifelse(pttest<alpha,"significant","non-significant")
  TOSToutcome<-ifelse(ptost<alpha,"significant","non-significant")
  plot(NA, ylim=c(0,1), xlim=c(min(LL95,low_eqbound_d,ES)-max(UL95-LL95, high_eqbound_d-low_eqbound_d,ES)/10, max(UL95,high_eqbound_d,ES)+max(UL95-LL95, high_eqbound_d-low_eqbound_d, ES)/10), bty="l", yaxt="n", ylab="",xlab="Effect size")
  points(x=ES, y=0.5, pch=15, cex=2)
  abline(v=high_eqbound_d, lty=2)
  abline(v=low_eqbound_d, lty=2)
  abline(v=0, lty=2, col="grey")
  segments(LL90,0.5,UL90,0.5, lwd=3)
  segments(LL95,0.5,UL95,0.5, lwd=1)
  title(main=paste("Equivalence bounds ",round(low_eqbound_d,digits=3)," and ",round(high_eqbound_d,digits=3),"\nEffect size = ",round(ES,digits=3)," \n TOST: 90% CI [",round(LL90,digits=3),";",round(UL90,digits=3),"] ", TOSToutcome," \n NHST: 95% CI [",round(LL95,digits=3),";",round(UL95,digits=3),"] ", testoutcome,sep=""), cex.main=1)
  message(cat("Using alpha = ",alpha," the meta-analysis was ",testoutcome,", Z = ",Z,", p = ",pttest,sep=""))
  message(cat("Using alpha = ",alpha," the equivalence test was ",TOSToutcome,", Z = ",Ztost,", p = ",ptost,sep=""))
  return(results)
}