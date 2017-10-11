#' TOST function for two proportions (raw scores)
#' @param prop1 proportion of group 1
#' @param prop2 proportion of group 2
#' @param n1 sample size in group 1
#' @param n2 sample size in group 2
#' @param low_eqbound lower equivalence bounds (e.g., -0.1) expressed in proportions
#' @param high_eqbound upper equivalence bounds (e.g., 0.1) expressed in proportions
#' @param alpha alpha level (default = 0.05)
#' @return Returns TOST z-value 1, TOST p-value 1, TOST z-value 2, TOST p-value 2, low equivalence bound, high equivalence bound, Lower limit confidence interval TOST, Upper limit confidence interval TOST
#' @importFrom stats pnorm pt qnorm qt
#' @importFrom graphics abline plot points segments title
#' @examples
#' ## Equivalence test for two independent proportions equal to .65 and .70, with 100 samples per group, a lower equivalence bound of -0.1,
#' ## a higher equivalence bound of 0.1, and an alpha of 0.05.
#' TOSTtwo.prop(prop1 = .65, prop2 = .70, n1 = 100, n2 = 100, low_eqbound = -0.1, high_eqbound = 0.1, alpha = .05)
#' @section References:
#' Tunes da Silva, G., Logan, B. R., & Klein, J. P. (2008). Methods for Equivalence and Noninferiority Testing. Biology of Blood Marrow Transplant, 15(1 Suppl), 120-127.
#' Yin, G. (2012). Clinical Trial Design: Bayesian and Frequentist Adaptive Methods. Hoboken, New Jersey: John Wiley & Sons, Inc.
#' @export

TOSTtwo.prop <- function(prop1, prop2, n1, n2, low_eqbound, high_eqbound, alpha) {
  if(missing(alpha)) {
    alpha <- 0.05
  }
  prop_dif <- prop1 - prop2
  prop_se <- sqrt((prop1*(1-prop1))/n1 + (prop2*(1-prop2))/n2)

  #calculating z-statistic
  z1 <- (prop_dif - low_eqbound)/prop_se
  z2 <- (prop_dif - high_eqbound)/prop_se

  #calculating p-value for both one-sided tests
  p1 <- 1 - pnorm(z1)
  p2 <- pnorm(z2)
  ptost <- max(p1,p2) #Get highest p-value for summary TOST result
  ztost <- ifelse(abs(z1) < abs(z2), z1, z2) #Get lowest z-value for summary TOST result
  TOSToutcome <- ifelse(ptost<alpha,"significant","non-significant")

  #calculating CIs
  CI_lb <- prop_dif - (qnorm(1-alpha) * prop_se)
  CI_ub <- prop_dif + (qnorm(1-alpha) * prop_se)

  #graph with 90% CIs
  plot(NA, ylim=c(0,1), xlim=c(min(CI_lb,(low_eqbound))-max(CI_ub-CI_lb, high_eqbound-(low_eqbound))/10, max(CI_ub,high_eqbound)+max(CI_ub-CI_lb, high_eqbound-(low_eqbound))/10), bty="l", yaxt="n", ylab="",xlab="Proportion Difference")
  points(x=prop_dif, y=0.5, pch=15, cex=2)
  abline(v=high_eqbound, lty=2)
  abline(v=low_eqbound, lty=2)
  abline(v=0, lty=2, col="grey")
  segments(CI_lb,0.5,CI_ub,0.5, lwd=3)
  title(main=paste("Equivalence bounds ",round(low_eqbound,digits=3)," and ",round(high_eqbound,digits=3),"\nProportion Difference = ",round(prop_dif,3)," \n TOST: ", 100*(1-alpha*2),"% CI [",round(CI_lb,digits=3),";",round(CI_ub,digits=3),"] ", TOSToutcome, sep=""), cex.main=1)

  TOSTresults<-data.frame(z1,p1,z2,p2)
  colnames(TOSTresults) <- c("z-value 1","p-value 1","z-value 2","p-value 2")
  bound_results<-data.frame(low_eqbound,high_eqbound)
  colnames(bound_results) <- c("low bound","high bound")
  CIresults<-data.frame(CI_lb,CI_ub)
  colnames(CIresults) <- c(paste("Lower Limit ",100*(1-alpha*2),"% CI",sep=""),paste("Upper Limit ",100*(1-alpha*2),"% CI",sep=""))
  cat("TOST results:\n")
  print(TOSTresults)
  cat("\n")
  cat("Equivalence bounds:\n")
  print(bound_results)
  cat("\n")
  cat("TOST confidence interval:\n")
  print(CIresults)
  invisible(list(TOST_z1=z1,TOST_p1=p1,TOST_z2=z2,TOST_p2=p2, alpha=alpha,low_eqbound=low_eqbound,high_eqbound=high_eqbound, LL_CI_TOST=CI_lb,UL_CI_TOST=CI_ub))
}
