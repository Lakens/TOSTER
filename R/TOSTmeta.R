#' TOST function for meta-analysis
#' @param ES meta-analytic effect size
#' @param var meta-analytic variance
#' @param se standard error
#' @param low_eqbound_d lower equivalence bounds (e.g., -0.5) expressed in standardized mean difference (Cohen's d)
#' @param high_eqbound_d upper equivalence bounds (e.g., 0.5) expressed in standardized mean difference (Cohen's d)
#' @param alpha alpha level (default = 0.05)
#' @return Returns TOST Z-value 1, TOST p-value 1, TOST Z-value 2, TOST p-value 2,  alpha, low equivalence bound d, high equivalence bound d, Lower limit confidence interval TOST, Upper limit confidence interval TOST
#' @examples
#' ## Run TOSTmeta by specifying the standard error
#' TOSTmeta(ES=0.12, se=0.09, low_eqbound_d=-0.2, high_eqbound_d=0.2, alpha=0.05)
#' ## Run TOSTmeta by specifying the variance
#' TOSTmeta(ES=0.12, var=0.0081, low_eqbound_d=-0.2, high_eqbound_d=0.2, alpha=0.05)
#' ## If both variance and se are specified, TOSTmeta will use standard error and ignore variance
#' TOSTmeta(ES=0.12, var=9999, se = 0.09, low_eqbound_d=-0.2, high_eqbound_d=0.2, alpha=0.05)
#' @section References:
#' Rogers, J. L., Howard, K. I., & Vessey, J. T. (1993). Using significance tests to evaluate equivalence between two experimental groups. Psychological Bulletin, 113(3), 553, formula page 557.
#' @export
#'

TOSTmeta<-function(ES,var,se,low_eqbound_d, high_eqbound_d, alpha, plot = TRUE){
  if(missing(alpha)) {
    alpha<-0.05
  }
  if(missing(se)) {
    if(missing(var)) {
      stop("Need to specify variance (var) or standard error (se).")
    }
    se<-sqrt(var)
  }
  if(missing(var)) {
    if(missing(se)) {
      stop("Need to specify variance (var) or standard error (se).")
    }
  }
  # Calculate TOST, z-test, 90% CIs and 95% CIs
  Z1<-(ES-low_eqbound_d)/se
  p1<-pnorm(Z1, lower.tail=FALSE)
  Z2<-(ES-high_eqbound_d)/se
  p2<-pnorm(Z2, lower.tail=TRUE)
  Z<-(ES/se)
  pttest<-2*pnorm(-abs(Z))
  LL90<-ES-qnorm(1-alpha)*(se)
  UL90<-ES+qnorm(1-alpha)*(se)
  LL95<-ES-qnorm(1-alpha/2)*(se)
  UL95<-ES+qnorm(1-alpha/2)*(se)
  ptost<-max(p1,p2) #Get highest p-value for summary TOST result
  Ztost<-ifelse(abs(Z1) < abs(Z2), Z1, Z2) #Get lowest t-value for summary TOST result
  results<-data.frame(Z1,p1,Z2,p2,LL90,UL90)
  colnames(results) <- c("Z-value 1","p-value 1","Z-value 2","p-value 2", paste("Lower Limit ",100*(1-alpha*2),"% CI",sep=""),paste("Upper Limit ",100*(1-alpha*2),"% CI",sep=""))
  testoutcome<-ifelse(pttest<alpha,"significant","non-significant")
  TOSToutcome<-ifelse(ptost<alpha,"significant","non-significant")

  # Plot results
  if (plot == TRUE) {
  plot(NA, ylim=c(0,1), xlim=c(min(LL95,low_eqbound_d,ES)-max(UL95-LL95, high_eqbound_d-low_eqbound_d,ES)/10, max(UL95,high_eqbound_d,ES)+max(UL95-LL95, high_eqbound_d-low_eqbound_d, ES)/10), bty="l", yaxt="n", ylab="",xlab="Effect size")
  points(x=ES, y=0.5, pch=15, cex=2)
  abline(v=high_eqbound_d, lty=2)
  abline(v=low_eqbound_d, lty=2)
  abline(v=0, lty=2, col="grey")
  segments(LL90,0.5,UL90,0.5, lwd=3)
  segments(LL95,0.5,UL95,0.5, lwd=1)
  title(main=paste("Equivalence bounds ",round(low_eqbound_d,digits=3)," and ",round(high_eqbound_d,digits=3),"\nEffect size = ",round(ES,digits=3)," \n TOST: ", 100*(1-alpha*2),"% CI [",round(LL90,digits=3),";",round(UL90,digits=3),"] ", TOSToutcome," \n NHST: ", 100*(1-alpha),"% CI [",round(LL95,digits=3),";",round(UL95,digits=3),"] ", testoutcome,sep=""), cex.main=1)
  }

  # Print TOST and z-test results in message form
  message(cat("Using alpha = ",alpha," the meta-analysis was ",testoutcome,", Z = ",Z,", p = ",pttest,sep=""))
  cat("\n")
  message(cat("Using alpha = ",alpha," the equivalence test was ",TOSToutcome,", Z = ",Ztost,", p = ",ptost,sep=""))

  # Print TOST and z test results in table form
  TOSTresults<-data.frame(Z1,p1,Z2,p2)
  colnames(TOSTresults) <- c("Z-value 1","p-value 1","Z-value 2","p-value 2")
  bound_d_results<-data.frame(low_eqbound_d,high_eqbound_d)
  colnames(bound_d_results) <- c("low bound d","high bound d")
  CIresults<-data.frame(LL90,UL90)
  colnames(CIresults) <- c(paste("Lower Limit ",100*(1-alpha*2),"% CI",sep=""),paste("Upper Limit ",100*(1-alpha*2),"% CI",sep=""))
  cat("TOST results:\n")
  print(TOSTresults)
  cat("\n")
  cat("Equivalence bounds (Cohen's d):\n")
  print(bound_d_results)
  cat("\n")
  cat("TOST confidence interval:\n")
  print(CIresults)

  # Return results in list()
  invisible(list(ES=ES,TOST_Z1=Z1,TOST_p1=p1,TOST_Z2=Z2,TOST_p2=p2,alpha=alpha,low_eqbound_d=low_eqbound_d,high_eqbound_d=high_eqbound_d, LL_CI_TOST=LL90,UL_CI_TOST=UL90,LL_CI_ZTEST=LL95,UL_CI_ZTEST=UL95))
}
