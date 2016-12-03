#' TOST function for a one-sample t-test (raw scores)
#' @param m mean
#' @param mu value to compare against
#' @param sd standard deviation
#' @param n sample size
#' @param low_eqbound_d lower equivalence bounds (e.g., -0.5) expressed in standardized mean difference (Cohen's d)
#' @param high_eqbound_d upper equivalence bounds (e.g., 0.5) expressed in standardized mean difference (Cohen's d)
#' @param alpha alpha level (default = 0.05)
#' @return Returns TOST t-value 1, TOST p-value 1, TOST t-value 2, TOST p-value 2, degrees of freedom, low equivalence bound, high equivalence bound, Lower limit confidence interval TOST, Upper limit confidence interval TOST
#' @examples
#' TOSTone.raw(m=0.52,mu=0.5,sd=0.5,n=300,low_eqbound=-0.1, high_eqbound=0.1, alpha=0.05)
#' @export
#' 

TOSTone.raw<-function(m,mu,sd,n,low_eqbound, high_eqbound, alpha){
  if(missing(alpha)) {
    alpha<-0.05
  }
  degree_f<-n-1
  t1<-(m-mu-low_eqbound)/(sd/sqrt(n))# t-test
  p1<-pt(t1, degree_f, lower=FALSE) 
  t2<-(m-mu-high_eqbound)/(sd/sqrt(n)) #t-test
  p2<-pt(t2, degree_f, lower=TRUE) 
  t<-(m-mu)/(sd/sqrt(n))
  pttest<-2*pt(-abs(t), df=degree_f)
  LL90<-(m-mu)-qt(1-alpha, degree_f)*(sd/sqrt(n))
  UL90<-(m-mu)+qt(1-alpha, degree_f)*(sd/sqrt(n))
  LL95<-(m-mu)-qt(1-(alpha/2), degree_f)*(sd/sqrt(n))
  UL95<-(m-mu)+qt(1-(alpha/2), degree_f)*(sd/sqrt(n))
  ptost<-max(p1,p2) #Get highest p-value for summary TOST result
  ttost<-ifelse(abs(t1) < abs(t2), t1, t2) #Get lowest t-value for summary TOST result
  results<-data.frame(t1,p1,t2,p2,degree_f,LL90,UL90)
  colnames(results) <- c("t-value 1","p-value 1","t-value 2","p-value 2","df", paste("Lower Limit ",100*(1-alpha*2),"% CI",sep=""),paste("Upper Limit ",100*(1-alpha*2),"% CI",sep=""))
  dif<-(m-mu)
  testoutcome<-ifelse(pttest<alpha,"significant","non-significant")
  TOSToutcome<-ifelse(ptost<alpha,"significant","non-significant")
  plot(NA, ylim=c(0,1), xlim=c(min(LL95,low_eqbound,dif)-max(UL95-LL95, high_eqbound-low_eqbound,dif)/10, max(UL95,high_eqbound,dif)+max(UL95-LL95, high_eqbound-low_eqbound, dif)/10), bty="l", yaxt="n", ylab="",xlab="Mean Difference")
  points(x=dif, y=0.5, pch=15, cex=2)
  abline(v=high_eqbound, lty=2)
  abline(v=low_eqbound, lty=2)
  abline(v=0, lty=2, col="grey")
  segments(LL90,0.5,UL90,0.5, lwd=3)
  segments(LL95,0.5,UL95,0.5, lwd=1)
  title(main=paste("Equivalence bounds ",round(low_eqbound,digits=3)," and ",round(high_eqbound,digits=3),"\nMean difference = ",round(dif,digits=3)," \n TOST: ", 100*(1-alpha*2),"% CI [",round(LL90,digits=3),";",round(UL90,digits=3),"] ", TOSToutcome," \n NHST: ", 100*(1-alpha),"% CI [",round(LL95,digits=3),";",round(UL95,digits=3),"] ", testoutcome,sep=""), cex.main=1)
  message(cat("Using alpha = ",alpha," the NHST one-sample t-test was ",testoutcome,", t(",degree_f,") = ",t,", p = ",pttest,sep=""))
  message(cat("Using alpha = ",alpha," the equivalence test was ",TOSToutcome,", t(",degree_f,") = ",ttost,", p = ",ptost,sep=""))
  TOSTresults<-data.frame(t1,p1,t2,p2,degree_f)
  colnames(TOSTresults) <- c("t-value 1","p-value 1","t-value 2","p-value 2","df")
  bound_results<-data.frame(low_eqbound,high_eqbound)
  colnames(bound_results) <- c("low bound raw","high bound raw")
  CIresults<-data.frame(LL90,UL90)
  colnames(CIresults) <- c(paste("Lower Limit ",100*(1-alpha*2),"% CI raw",sep=""),paste("Upper Limit ",100*(1-alpha*2),"% CI raw",sep=""))
  cat("TOST results:\n")
  print(TOSTresults)  
  cat("\n")
  cat("Equivalence bounds (raw scores):\n")
  print(bound_results)  
  cat("\n")
  cat("TOST confidence interval:\n")
  print(CIresults)
  invisible(list(TOST_t1=t1,TOST_p1=p1,TOST_t2=t2,TOST_p2=p2, TOST_df=degree_f,alpha=alpha,low_eqbound=low_eqbound,high_eqbound=high_eqbound, LL_CI_TOST=LL90,UL_CI_TOST=UL90))
}