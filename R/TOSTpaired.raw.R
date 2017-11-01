#' TOST function for a dependent t-test (raw scores)
#' @param m1 mean of group 1
#' @param m2 mean of group 2
#' @param sd1 standard deviation of group 1
#' @param sd2 standard deviation of group 2
#' @param n sample size (pairs)
#' @param r12 correlation of dependent variable between group 1 and  group 2
#' @param low_eqbound lower equivalence bounds (e.g., -0.5) expressed in raw scores
#' @param high_eqbound upper equivalence bounds (e.g., 0.5) expressed in raw scores
#' @param alpha alpha level (default = 0.05)
#' @param plot set whether results should be plotted (plot = TRUE) or not (plot = FALSE) - defaults to TRUE
#' @return Returns TOST t-value 1, TOST p-value 1, TOST t-value 2, TOST p-value 2, degrees of freedom, low equivalence bound, high equivalence bound, Lower limit confidence interval TOST, Upper limit confidence interval TOST
#' @examples
#' ## Test means of 5.83 and 5.75, standard deviations of 1.17 and 1.30 in sample of 65 pairs
#' ## with correlation between observations of 0.745 using equivalence bounds in raw units of
#' ## -0.34 and 0.34, (with default alpha setting of = 0.05).
#' TOSTpaired.raw(n=65,m1=5.83,m2=5.75,sd1=1.17,sd2=1.30,r12=0.745,low_eqbound=-0.34,high_eqbound=0.34)
#' @section References:
#' Mara, C. A., & Cribbie, R. A. (2012). Paired-Samples Tests of Equivalence. Communications in Statistics - Simulation and Computation, 41(10), 1928-1943. https://doi.org/10.1080/03610918.2011.626545, formula page 1932. Note there is a typo in the formula: n-1 should be n (personal communication, 31-8-2016)
#' @importFrom stats pnorm pt qnorm qt
#' @importFrom graphics abline plot points segments title
#' @export
#'

TOSTpaired.raw<-function(n,m1,m2,sd1,sd2,r12,low_eqbound, high_eqbound, alpha, plot = TRUE){
  if(missing(alpha)) {
    alpha <- 0.05
  }

  # Calculate TOST, t-test, 90% CIs and 95% CIs
  sdif<-sqrt(sd1^2+sd2^2-2*r12*sd1*sd2)
  se<-sdif/sqrt(n)
  t<-(m1-m2)/se
  degree_f<-n-1
  pttest<-2*pt(abs(t), degree_f, lower.tail=FALSE)
  t1<-((m1-m2)-(low_eqbound*sdif))/se
  p1<-pt(t1, degree_f, lower.tail=FALSE)
  t2<-((m1-m2)-(high_eqbound*sdif))/se
  p2<-pt(t2, degree_f, lower.tail=TRUE)
  ttost<-ifelse(abs(t1) < abs(t2), t1, t2)
  LL90<-((m1-m2)-qt(1-alpha, degree_f)*se)
  UL90<-((m1-m2)+qt(1-alpha, degree_f)*se)
  ptost<-max(p1,p2)
  dif<-(m1-m2)
  LL95<-((m1-m2)-qt(1-(alpha/2), degree_f)*se)
  UL95<-((m1-m2)+qt(1-(alpha/2), degree_f)*se)
  testoutcome<-ifelse(pttest<alpha,"significant","non-significant")
  TOSToutcome<-ifelse(ptost<alpha,"significant","non-significant")

  # Plot results
  if (plot == TRUE) {
  plot(NA, ylim=c(0,1), xlim=c(min(LL90,low_eqbound)-max(UL90-LL90, high_eqbound-low_eqbound)/10, max(UL90,high_eqbound)+max(UL90-LL90, high_eqbound-low_eqbound)/10), bty="l", yaxt="n", ylab="",xlab="Mean Difference")
  points(x=dif, y=0.5, pch=15, cex=2)
  abline(v=high_eqbound, lty=2)
  abline(v=low_eqbound, lty=2)
  abline(v=0, lty=2, col="grey")
  segments(LL90,0.5,UL90,0.5, lwd=3)
  segments(LL95,0.5,UL95,0.5, lwd=1)
  title(main=paste("Equivalence bounds ",round(low_eqbound,digits=3)," and ",round(high_eqbound,digits=3),"\nMean difference = ",round(dif,digits=3)," \n TOST: ", 100*(1-alpha*2),"% CI [",round(LL90,digits=3),";",round(UL90,digits=3),"] ", TOSToutcome," \n NHST: ", 100*(1-alpha),"% CI [",round(LL95,digits=3),";",round(UL95,digits=3),"] ", testoutcome,sep=""), cex.main=1)
  }

  # Print TOST and t-test results in message form
  message(cat("Using alpha = ",alpha," the NHST t-test was ",testoutcome,", t(",degree_f,") = ",t,", p = ",pttest,sep=""))
  cat("\n")
  message(cat("Using alpha = ",alpha," the equivalence test was ",TOSToutcome,", t(",degree_f,") = ",ttost,", p = ",ptost,sep=""))

  # Print TOST and t-test results in table form
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

  # Print TOST and t-test results in table form
  invisible(list(diff=dif,TOST_t1=t1,TOST_p1=p1,TOST_t2=t2,TOST_p2=p2, TOST_df=degree_f,alpha=alpha,low_eqbound=low_eqbound,high_eqbound=high_eqbound, LL_CI_TOST=LL90,UL_CI_TOST=UL90, LL_CI_TTEST=LL95, UL_CI_TTEST=UL95))
}
