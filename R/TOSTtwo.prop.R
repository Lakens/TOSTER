#' TOST function for two propotions (raw scores)
#' @param proprop1 proportion of group 1
#' @param proprop2 proportion of group 2
#' @param n1 sample size in group 1
#' @param n2 sample size in group 2
#' @param threshold equivalence bounds (e.g., 0.1) expressed in raw scale units (e.g., scalepoints)
#' @param alpha alpha level (default = 0.05)
#' @return Returns TOST z-value 1, TOST p-value 1, TOST z-value 2, TOST p-value 2, low equivalence bound, high equivalence bound, Lower limit confidence interval TOST, Upper limit confidence interval TOST
#' @importFrom stats pnorm pt qnorm qt
#' @importFrom graphics abline plot points segments title
#' @examples
#' ## Eskine (2013) showed that participants who had been exposed to organic
#' ## food were substantially harsher in their moral judgments relative to
#' ## those exposed to control (d = 0.81, 95% CI: [0.19, 1.45]). A
#' ## replication by Moery & Calin-Jageman (2016, Study 2) did not observe
#' ## a significant effect (Control: n = 95, M = 5.25, SD = 0.95, Organic
#' ## Food: n = 89, M = 5.22, SD = 0.83). Following Simonsohn's (2015)
#' ## recommendation the equivalence bound was set to the effect size the
#' ## original study had 33% power to detect (with n = 21 in each condition,
#' ## this means the equivalence bound is d = 0.48, which equals a
#' ## difference of 0.384 on a 7-point scale given the sample sizes and a
#' ## pooled standard deviation of 0.894). Using a TOST equivalence test
#' ## with alpha = 0.05, assuming equal variances, and equivalence
#' ## bounds of d = -0.43 and d = 0.43 is significant, t(182) = -2.69,
#' ## p = 0.004. We can reject effects larger than d = 0.43.
#'
#' TOSTtwo.raw(m1=5.25,m2=5.22,sd1=0.95,sd2=0.83,n1=95,n2=89,low_eqbound=-0.384,high_eqbound=0.384)
#' @section References:
#' Tunes da Silva, G., Logan, B. R., & Klein, J. P. (2008). Methods for Equivalence and Noninferiority Testing. Biology of Blood Marrow Transplant, 15(1 Suppl), 120-127.
#'
#' @export

TOSTtwo.prop <- function(prop1, prop2, n1, n2, threshold, alpha) {
  if(missing(alpha)) {
    alpha<-0.05
  }
  if(threshold < 0) {
    alpha<-abs(threshold)
  }
  
  prop_dif <- prop1 - prop2
  prop_se <- sqrt((prop1*(1-prop1))/n1 + (prop2*(1-prop2))/n2)
  
  ##defining z-statistic
  z_lower_numerator <- prop_dif + threshold
  z_upper_numerator <- prop_dif - threshold
  z_lower <- z_lower_numerator/prop_se
  z_upper <- z_upper_numerator/prop_se
  
  ##calculating p-value for both one-sided tests
  p_lower <- 1 - pnorm(z_lower)
  p_upper <- pnorm(z_upper)
  
  ##calculating CIs
  CI_lb <- prop_dif - (qnorm(1-alpha) * prop_se)
  CI_ub <- prop_dif + (qnorm(1-alpha) * prop_se)
  
  ##graph with 90% CIs
  plot(NA, ylim=c(0,1), xlim=c(min(CI_lb,(-threshold))-max(CI_ub-CI_lb, threshold-(-threshold))/10, max(CI_ub,threshold)+max(CI_ub-CI_lb, threshold-(-threshold))/10), bty="l", yaxt="n", ylab="",xlab="Proportion Difference")
  points(x=prop_dif, y=0.5, pch=15, cex=2)
  abline(v=threshold, lty=2)
  abline(v=-threshold, lty=2)
  abline(v=0, lty=2, col="grey")
  segments(CI_lb,0.5,CI_ub,0.5, lwd=3)
  title(main=paste("Equivalence bounds ",round(-threshold,digits=3)," and ",round(threshold,digits=3),"\nProportion Difference = "," \n TOST: ", 100*(1-alpha*2),"% CI [",round(CI_lb,digits=3),";",round(CI_ub,digits=3),"] ", cex.main=1))
        
  TOSTresults<-data.frame(z_lower,p_lower,z_upper,p_upper)
  colnames(TOSTresults) <- c("z-value 1","p-value 1","z-value 2","p-value 2")
  bound_results<-data.frame(-threshold,threshold)
  colnames(bound_results) <- c("low bound raw","high bound raw")
  CIresults<-data.frame(CI_lb,CI_ub)
  colnames(CIresults) <- c(paste("Lower Limit ",100*(1-alpha*2),"% CI raw",sep=""),paste("Upper Limit ",100*(1-alpha*2),"% CI raw",sep=""))
  cat("TOST results:\n")
  print(TOSTresults)
  cat("\n")
  cat("Equivalence bounds (raw scores):\n")
  print(bound_results)
  cat("\n")
  cat("TOST confidence interval:\n")
  print(CIresults)
  invisible(list(TOST_z1=z_lower,TOST_p1=p_lower,TOST_z2=z_upper,TOST_p2=p_upper, alpha=alpha,low_eqbound=-threshold,high_eqbound=threshold, LL_CI_TOST=CI_lb,UL_CI_TOST=CI_ub))
}
