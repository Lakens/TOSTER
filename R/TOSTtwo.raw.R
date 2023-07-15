#' TOST function for an independent t-test (raw scores)
#' @description
#' `r lifecycle::badge('superseded')`
#'
#' Development on `TOSTtwo.raw` is complete, and for new code we recommend
#' switching to [tsum_TOST], which is easier to use, more featureful,
#' and still under active development.
#'
#' @param m1 mean of group 1
#' @param m2 mean of group 2
#' @param sd1 standard deviation of group 1
#' @param sd2 standard deviation of group 2
#' @param n1 sample size in group 1
#' @param n2 sample size in group 2
#' @param low_eqbound lower equivalence bounds (e.g., -0.5) expressed in raw scale units (e.g., scalepoints)
#' @param high_eqbound upper equivalence bounds (e.g., 0.5) expressed in raw scale units (e.g., scalepoints)
#' @param alpha alpha level (default = 0.05)
#' @param var.equal logical variable indicating whether equal variances assumption is assumed to be TRUE or FALSE. Defaults to FALSE.
#' @param plot set whether results should be plotted (plot = TRUE) or not (plot = FALSE) - defaults to TRUE
#' @param verbose logical variable indicating whether text output should be generated (verbose = TRUE) or not (verbose = FALSE) - default to TRUE
#' @return Returns TOST t-value 1, TOST p-value 1, TOST t-value 2, TOST p-value 2, degrees of freedom, low equivalence bound, high equivalence bound, Lower limit confidence interval TOST, Upper limit confidence interval TOST
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
#' Berger, R. L., & Hsu, J. C. (1996). Bioequivalence Trials, Intersection-Union Tests and Equivalence Confidence Sets. Statistical Science, 11(4), 283-302.
#'
#' Gruman, J. A., Cribbie, R. A., & Arpin-Cribbie, C. A. (2007). The effects of heteroscedasticity on tests of equivalence. Journal of Modern Applied Statistical Methods, 6(1), 133-140, formula for Welch's t-test on page 135
#' @export



TOSTtwo.raw<-function(m1,m2,sd1,sd2,n1,n2,low_eqbound, high_eqbound, alpha, var.equal, plot = TRUE, verbose = TRUE){
  message("Note: this function is defunct. Please use tsum_TOST instead")
  if(missing(alpha)) {
    alpha<-0.05
  }
  if(missing(var.equal)) {
    var.equal<-FALSE
  }
  if(low_eqbound >= high_eqbound) warning("The lower bound is equal to or larger than the upper bound. Check the plot and output to see if the bounds are specified as you intended.")
  if(n1 < 2 | n2 < 2) stop("The sample size should be larger than 1.")
  if(1 <= alpha | alpha <= 0) stop("The alpha level should be a positive value between 0 and 1.")
  if(sd1 <= 0 | sd2 <= 0) stop("The standard deviation should be a positive value.")

  # Calculate TOST, t-test, 90% CIs and 95% CIs
  if(var.equal==TRUE) {
    sdpooled<-sqrt((((n1 - 1)*(sd1^2)) + (n2 - 1)*(sd2^2))/((n1+n2)-2)) #calculate sd pooled
    t1<-((m1-m2)-low_eqbound)/(sdpooled*sqrt(1/n1 + 1/n2))
    degree_f<-n1+n2-2
    p1<-pt(t1, degree_f, lower.tail=FALSE)
    t2<-((m1-m2)-high_eqbound)/(sdpooled*sqrt(1/n1 + 1/n2))
    p2<-pt(t2, degree_f, lower.tail=TRUE)
    LL90<-(m1-m2)-qt(1-alpha, degree_f)*(sdpooled*sqrt(1/n1 + 1/n2))
    UL90<-(m1-m2)+qt(1-alpha, degree_f)*(sdpooled*sqrt(1/n1 + 1/n2))
    LL95<-(m1-m2)-qt(1-(alpha/2), degree_f)*(sdpooled*sqrt(1/n1 + 1/n2))
    UL95<-(m1-m2)+qt(1-(alpha/2), degree_f)*(sdpooled*sqrt(1/n1 + 1/n2))
    t<-(m1-m2)/(sdpooled*sqrt(1/n1 + 1/n2))
    pttest<-2*pt(-abs(t), df=degree_f)
  } else {
    sdpooled<-sqrt((sd1^2 + sd2^2)/2) #calculate sd root mean squared for Welch's t-test
    t1<-((m1-m2)-low_eqbound)/sqrt(sd1^2/n1 + sd2^2/n2) #welch's t-test lower bound
    degree_f<-(sd1^2/n1+sd2^2/n2)^2/(((sd1^2/n1)^2/(n1-1))+((sd2^2/n2)^2/(n2-1))) #degrees of freedom for Welch's t-test
    p1<-pt(t1, degree_f, lower.tail=FALSE) #p-value for Welch's TOST t-test
    t2<-((m1-m2)-high_eqbound)/sqrt(sd1^2/n1 + sd2^2/n2) #welch's t-test upper bound
    p2<-pt(t2, degree_f, lower.tail=TRUE) #p-value for Welch's TOST t-test
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

  if(missing(verbose)) {
    verbose <- TRUE
  }
  if(verbose == TRUE){
    cat("TOST results:\n")
    cat("t-value lower bound:",format(t1, digits = 3, nsmall = 2, scientific = FALSE),"\tp-value lower bound:",format(p1, digits = 1, nsmall = 3, scientific = FALSE))
    cat("\n")
    cat("t-value upper bound:",format(t2, digits = 3, nsmall = 2, scientific = FALSE),"\tp-value upper bound:",format(p2, digits = 1, nsmall = 3, scientific = FALSE))
    cat("\n")
    cat("degrees of freedom :",round(degree_f, digits = 2))
    cat("\n\n")
    cat("Equivalence bounds (raw scores):")
    cat("\n")
    cat("low eqbound:", paste0(round(low_eqbound, digits = 4)),"\nhigh eqbound:",paste0(round(high_eqbound, digits = 4)))
    cat("\n\n")
    cat("TOST confidence interval:")
    cat("\n")
    cat("lower bound ",100*(1-alpha*2),"% CI: ", paste0(round(LL90, digits = 3)),"\nupper bound ",100*(1-alpha*2),"% CI:  ",paste0(round(UL90,digits = 3)), sep = "")
    cat("\n\n")
    cat("NHST confidence interval:")
    cat("\n")
    cat("lower bound ",100*(1-alpha),"% CI: ", paste0(round(LL95, digits = 3)),"\nupper bound ",100*(1-alpha),"% CI:  ",paste0(round(UL95,digits = 3)), sep = "")
    cat("\n\n")
    cat("Equivalence Test Result:\n")
    message(cat("The equivalence test was ",TOSToutcome,", t(",round(degree_f, digits=2),") = ",format(ttost, digits = 3, nsmall = 3, scientific = FALSE),", p = ",format(ptost, digits = 3, nsmall = 3, scientific = FALSE),", given equivalence bounds of ",format(low_eqbound, digits = 3, nsmall = 3, scientific = FALSE)," and ",format(high_eqbound, digits = 3, nsmall = 3, scientific = FALSE)," (on a raw scale) and an alpha of ",alpha,".",sep=""))
    cat("\n")
    cat("Null Hypothesis Test Result:\n")
    message(cat("The null hypothesis test was ",testoutcome,", t(",round(degree_f, digits=2),") = ",format(t, digits = 3, nsmall = 3, scientific = FALSE),", p = ",format(pttest, digits = 3, nsmall = 3, scientific = FALSE),", given an alpha of ",alpha,".",sep=""))
    if(pttest <= alpha && ptost <= alpha){
      combined_outcome <- paste0("NHST: reject null significance hypothesis that the effect is equal to ", 0," \n",
                                 "TOST: reject null equivalence hypothesis")
    }
    if(pttest < alpha && ptost > alpha){
      combined_outcome <- paste0("NHST: reject null significance hypothesis that the effect is equal to ",0," \n",
                                 "TOST: don't reject null equivalence hypothesis")
    }
    if(pttest > alpha && ptost <= alpha){
      combined_outcome <- paste0("NHST: don't reject null significance hypothesis that the effect is equal to ",0," \n",
                                 "TOST: reject null equivalence hypothesis")
    }
    if(pttest > alpha && ptost > alpha){
      combined_outcome <-  paste0("NHST: don't reject null significance hypothesis that the effect is equal to ",0," \n",
                                  "TOST: don't reject null equivalence hypothesis")
    }
    cat("\n")
    message(combined_outcome)
  }
  # Return results in list()
  invisible(list(diff=dif,TOST_t1=t1,TOST_p1=p1,TOST_t2=t2,TOST_p2=p2, TOST_df=degree_f,alpha=alpha,low_eqbound=low_eqbound,high_eqbound=high_eqbound, LL_CI_TOST=LL90,UL_CI_TOST=UL90,LL_CI_TTEST=LL95, UL_CI_TTEST=UL95,NHST_t = t, NHST_p = pttest))
}
