#' TOST function for two propotions (raw scores)
#' @param proprop1 proportion of group 1
#' @param proprop2 proportion of group 2
#' @param n1 sample size in group 1
#' @param n2 sample size in group 2
#' @param margin equivalence bounds (e.g., 0.1) expressed in raw scale units (e.g., scalepoints)
#' @param alpha alpha level (default = 0.05)
#' @return Returns TOST z-value 1, TOST p-value 1, TOST z-value 2, TOST p-value 2, low equivalence bound, high equivalence bound, Lower limit confidence interval TOST, Upper limit confidence interval TOST
#' @importFrom stats pnorm pt qnorm qt
#' @importFrom graphics abline plot points segments title
#' @examples
#' ## Equivalence test for two independent proportions equal to .65 and .70, with 100 samples per group, an equivalence margin of .1 
#' ## and alpha = 0.05.
#' TOSTtwo.prop(prop1 = .65, prop2 = .70, n1 = 100, n2 = 100, margin = .1, alpha = .05)
#' @section References:
#' Tunes da Silva, G., Logan, B. R., & Klein, J. P. (2008). Methods for Equivalence and Noninferiority Testing. Biology of Blood Marrow Transplant, 15(1 Suppl), 120-127.
#' Yin, G. (2012). Clinical Trial Design: Bayesian and Frequentist Adaptive Methods. Hoboken, New Jersey: John Wiley & Sons, Inc.
#' @export

TOSTtwo.prop <- function(prop1, prop2, n1, n2, margin, alpha) {
  if(missing(alpha)) {
    alpha <- 0.05
  }
  if(margin < 0) {
    margin <- abs(margin)
  }
  
  prop_dif <- prop1 - prop2
  prop_se <- sqrt((prop1*(1-prop1))/n1 + (prop2*(1-prop2))/n2)
  
  ##defining z-statistic
  z_lower_numerator <- prop_dif + margin
  z_upper_numerator <- prop_dif - margin
  z_lower <- z_lower_numerator/prop_se
  z_upper <- z_upper_numerator/prop_se
  
  ##calculating p-value for both one-sided tests
  p_lower <- 1 - pnorm(z_lower)
  p_upper <- pnorm(z_upper)
  
  ##calculating CIs
  CI_lb <- prop_dif - (qnorm(1-alpha) * prop_se)
  CI_ub <- prop_dif + (qnorm(1-alpha) * prop_se)
  
  ##graph with 90% CIs
  plot(NA, ylim=c(0,1), xlim=c(min(CI_lb,(-margin))-max(CI_ub-CI_lb, margin-(-margin))/10, max(CI_ub,margin)+max(CI_ub-CI_lb, margin-(-margin))/10), bty="l", yaxt="n", ylab="",xlab="Proportion Difference")
  points(x=prop_dif, y=0.5, pch=15, cex=2)
  abline(v=margin, lty=2)
  abline(v=-margin, lty=2)
  abline(v=0, lty=2, col="grey")
  segments(CI_lb,0.5,CI_ub,0.5, lwd=3)
  title(main=paste("Equivalence margin ",round(-margin,digits=3)," and ",round(margin,digits=3),"\nProportion Difference = "," \n TOST: ", 100*(1-alpha*2),"% CI [",round(CI_lb,digits=3),";",round(CI_ub,digits=3),"] ", cex.main=1))
        
  TOSTresults<-data.frame(z_lower,p_lower,z_upper,p_upper)
  colnames(TOSTresults) <- c("z-value 1","p-value 1","z-value 2","p-value 2")
  bound_results<-data.frame(-margin,margin)
  colnames(bound_results) <- c("low bound raw","high bound raw")
  CIresults<-data.frame(CI_lb,CI_ub)
  colnames(CIresults) <- c(paste("Lower Limit ",100*(1-alpha*2),"% CI raw",sep=""),paste("Upper Limit ",100*(1-alpha*2),"% CI raw",sep=""))
  cat("TOST results:\n")
  print(TOSTresults)
  cat("\n")
  cat("Equivalence margin (raw scores):\n")
  print(bound_results)
  cat("\n")
  cat("TOST confidence interval:\n")
  print(CIresults)
  invisible(list(TOST_z1=z_lower,TOST_p1=p_lower,TOST_z2=z_upper,TOST_p2=p_upper, alpha=alpha,low_eqbound=-margin,high_eqbound=margin, LL_CI_TOST=CI_lb,UL_CI_TOST=CI_ub))
}
