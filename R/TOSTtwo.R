#' TOST function for an independent t-test (Cohen's d)
#' @param m1 mean of group 1
#' @param m2 mean of group 2
#' @param sd1 standard deviation of group 1
#' @param sd2 standard deviation of group 2
#' @param n1 sample size in group 1
#' @param n2 sample size in group 2
#' @param low_eqbound_d lower equivalence bounds (e.g., -0.5) expressed in standardized mean difference (Cohen's d)
#' @param high_eqbound_d upper equivalence bounds (e.g., 0.5) expressed in standardized mean difference (Cohen's d)
#' @param alpha alpha level (default = 0.05)
#' @param var.equal logical variable indicating whether equal variances assumption is assumed to be TRUE or FALSE.  Defaults to FALSE.
#' @return Returns dataframe with t-value, degrees of freedom, p-value, and lower and upper limit of confidence interval
#' @examples
#' TOSTtwo(m1=7.83,m2=7.98,sd1=1.21,sd2=1.29,n1=400,n2=400,low_eqbound_d=-0.199897678576151, high_eqbound_d=0.199897678576151,alpha=0.05, var.equal=TRUE)
#' @export
#' 

TOSTtwo<-function(m1,m2,sd1,sd2,n1,n2,low_eqbound_d, high_eqbound_d, alpha, var.equal){
  if(missing(alpha)) {
    alpha<-0.05
  }
  if(missing(var.equal)) {
    var.equal<-FALSE
  }
  if(var.equal==TRUE) {
    sdpooled<-sqrt((((n1 - 1)*(sd1^2)) + (n2 - 1)*(sd2^2))/((n1+n2)-2)) #calculate sd pooled
    low_eqbound<-low_eqbound_d*sdpooled
    high_eqbound<-high_eqbound_d*sdpooled
    degree_f<-n1+n2-2
    t1<-(abs(m1-m2)-low_eqbound)/(sdpooled*sqrt(1/n1 + 1/n2))  #students t-test lower bound
    p1<-pt(t1, degree_f, lower=FALSE) 
    t2<-(abs(m1-m2)-high_eqbound)/(sdpooled*sqrt(1/n1 + 1/n2)) #students t-test upper bound
    p2<-pt(t2, degree_f, lower=TRUE) 
    t<-(m1-m2)/(sdpooled*sqrt(1/n1 + 1/n2))
    pttest<-2*pt(-abs(t), df=degree_f)
    LL90<-(m1-m2)-qt(1-alpha, n1+n2-2)*(sdpooled*sqrt(1/n1 + 1/n2))
    UL90<-(m1-m2)+qt(1-alpha, n1+n2-2)*(sdpooled*sqrt(1/n1 + 1/n2))
    LL95<-(m1-m2)-qt(1-(alpha/2), n1+n2-2)*(sdpooled*sqrt(1/n1 + 1/n2))
    UL95<-(m1-m2)+qt(1-(alpha/2), n1+n2-2)*(sdpooled*sqrt(1/n1 + 1/n2))
  } else {
    sdpooled<-sqrt((sd1^2 + sd2^2)/2) #calculate sd root mean squared for Welch's t-test
    low_eqbound<-low_eqbound_d*sdpooled
    high_eqbound<-high_eqbound_d*sdpooled
    degree_f<-(sd1^2/n1+sd2^2/n2)^2/(((sd1^2/n1)^2/(n1-1))+((sd2^2/n2)^2/(n2-1))) #degrees of freedom for Welch's t-test
    t1<-(abs(m1-m2)-low_eqbound)/sqrt(sd1^2/n1 + sd2^2/n2) #welch's t-test upper bound
    p1<-pt(t1, degree_f, lower=FALSE) #p-value for Welch's TOST t-test
    t2<-(abs(m1-m2)-high_eqbound)/sqrt(sd1^2/n1 + sd2^2/n2) #welch's t-test lower bound
    p2<-pt(t2, degree_f, lower=TRUE) #p-value for Welch's TOST t-test
    t<-(m1-m2)/sqrt(sd1^2/n1 + sd2^2/n2) #welch's t-test NHST
    pttest<-2*pt(-abs(t), df=degree_f) #p-value for Welch's t-test
    LL90<-(m1-m2)-qt(1-alpha, degree_f)*sqrt(sd1^2/n1 + sd2^2/n2) #Lower limit for CI Welch's t-test
    UL90<-(m1-m2)+qt(1-alpha, degree_f)*sqrt(sd1^2/n1 + sd2^2/n2) #Upper limit for CI Welch's t-test
    LL95<-(m1-m2)-qt(1-(alpha/2), degree_f)*sqrt(sd1^2/n1 + sd2^2/n2) #Lower limit for CI Welch's t-test
    UL95<-(m1-m2)+qt(1-(alpha/2), degree_f)*sqrt(sd1^2/n1 + sd2^2/n2) #Upper limit for CI Welch's t-test
  }
  ptost<-max(p1,p2) #Get highest p-value for summary TOST result
  ttost<-ifelse(abs(t1) < abs(t2), t1, t2) #Get lowest t-value for summary TOST result
  dif<-(m1-m2)
  testoutcome<-ifelse(pttest<alpha,"significant","non-significant")
  TOSToutcome<-ifelse(ptost<alpha,"significant","non-significant")
  plot(NA, ylim=c(0,1), xlim=c(min(LL90,low_eqbound)-max(UL90-LL90, high_eqbound-low_eqbound)/10, max(UL90,high_eqbound)+max(UL90-LL90, high_eqbound-low_eqbound)/10), bty="l", yaxt="n", ylab="",xlab="Mean Difference")
  points(x=dif, y=0.5, pch=15, cex=2)
  abline(v=high_eqbound, lty=2)
  abline(v=low_eqbound, lty=2)
  abline(v=0, lty=2, col="grey")
  segments(LL90,0.5,UL90,0.5, lwd=3)
  segments(LL95,0.5,UL95,0.5, lwd=1)
  title(main=paste("Equivalence bounds ",round(low_eqbound,digits=3)," and ",round(high_eqbound,digits=3),"\nMean difference = ",round(dif,digits=3)," \n TOST: 90% CI [",round(LL90,digits=3),";",round(UL90,digits=3),"] ", TOSToutcome," \n NHST: 95% CI [",round(LL95,digits=3),";",round(UL95,digits=3),"] ", testoutcome,sep=""), cex.main=1)
  if(var.equal == TRUE) {
    message(cat("Using alpha = ",alpha," Student's t-test was ",testoutcome,", t(",degree_f,") = ",t,", p = ",pttest,sep=""))
    message(cat("Using alpha = ",alpha," the equivalence test based on Student's t-test was ",TOSToutcome,", t(",degree_f,") = ",ttost,", p = ",ptost,sep=""))
  } else {
    message(cat("Using alpha = ",alpha," Welch's t-test was ",testoutcome,", t(",degree_f,") = ",t,", p = ",pttest,sep=""))
    message(cat("Using alpha = ",alpha," the equivalence test based on Welch's t-test  was ",TOSToutcome,", t(",degree_f,") = ",ttost,", p = ",ptost,sep=""))
  }
  TOSTresults<-data.frame(t1,p1,t2,p2,degree_f)
  colnames(TOSTresults) <- c("t-value 1","p-value 1","t-value 2","p-value 2","df")
  bound_d_results<-data.frame(low_eqbound_d,high_eqbound_d)
  colnames(bound_d_results) <- c("low bound d","high bound d")
  bound_results<-data.frame(low_eqbound,high_eqbound)
  colnames(bound_results) <- c("low bound raw","high bound raw")
  CIresults<-data.frame(LL90,UL90)
  colnames(CIresults) <- c(paste("Lower Limit ",100*(1-alpha*2),"% CI raw",sep=""),paste("Upper Limit ",100*(1-alpha*2),"% CI raw",sep=""))
  cat("TOST results:\n")
  print(TOSTresults)  
  cat("\n")
  cat("Equivalence bounds (Cohen's d):\n")
  print(bound_d_results)  
  cat("\n")
  cat("Equivalence bounds (raw scores):\n")
  print(bound_results)  
  cat("\n")
  cat("TOST confidence interval:\n")
  print(CIresults)
  invisible(list(TOST_t1=t1,TOST_p1=p1,TOST_t2=t2,TOST_p2=p2, TOST_df=degree_f,alpha=alpha,low_eqbound=low_eqbound,high_eqbound=high_eqbound,low_eqbound_d=low_eqbound_d,high_eqbound_d=high_eqbound_d, LL_CI_TOST=LL90,UL_CI_TOST=UL90))
}