#' TOST function for a one-sample t-test (raw scores)
#' @description
#' `r lifecycle::badge('superseded')`
#'
#' Development on this function is complete, and for new code we recommend
#' switching to [tsum_TOST], which is easier to use, more featureful,
#' and still under active development.
#'
#' @param m mean
#' @param mu value to compare against
#' @param sd standard deviation
#' @param n sample size
#' @param low_eqbound lower equivalence bounds (e.g., -0.5) expressed in raw units
#' @param high_eqbound upper equivalence bounds (e.g., 0.5) expressed in raw units
#' @param alpha alpha level (default = 0.05)
#' @param plot set whether results should be plotted (plot = TRUE) or not (plot = FALSE) - defaults to TRUE
#' @param verbose logical variable indicating whether text output should be generated (verbose = TRUE) or not (verbose = FALSE) - default to TRUE
#' @return Returns TOST t-value 1, TOST p-value 1, TOST t-value 2, TOST p-value 2, degrees of freedom, low equivalence bound, high equivalence bound, Lower limit confidence interval TOST, Upper limit confidence interval TOST
#' @examples
#' ## Test observed mean of 0.52 and standard deviation of 0.52 in sample of 300 participants
#' ## against 0.5 given equivalence bounds in raw units of -0.1 and 0.1, with an alpha = 0.05.
#' TOSTone.raw(m=0.52,mu=0.5,sd=0.5,n=300,low_eqbound=-0.1, high_eqbound=0.1, alpha=0.05)
#' @importFrom stats pnorm pt qnorm qt
#' @importFrom graphics abline plot points segments title
#' @export
#'

TOSTone.raw <-
  function(m,
           mu,
           sd,
           n,
           low_eqbound,
           high_eqbound,
           alpha,
           plot = TRUE,
           verbose = TRUE) {
    message("Note: this function is defunct. Please use tsum_TOST instead")
  if(missing(alpha)) {
    alpha<-0.05
  }
  if(low_eqbound >= high_eqbound) warning("The lower bound is equal to or larger than the upper bound. Check the plot and output to see if the bounds are specified as you intended.")
  if(n < 2) stop("The sample size should be larger than 1.")
  if(1 <= alpha | alpha <= 0) stop("The alpha level should be a positive value between 0 and 1.")
  if(sd <= 0) stop("The standard deviation should be a positive value.")
  # Calculate TOST, t-test, 90% CIs and 95% CIs
  degree_f<-n-1
  t1<-(m-mu-low_eqbound)/(sd/sqrt(n))# t-test
  p1<-pt(t1, degree_f, lower.tail=FALSE)
  t2<-(m-mu-high_eqbound)/(sd/sqrt(n)) #t-test
  p2<-pt(t2, degree_f, lower.tail=TRUE)
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

  # Plot results
  if (plot == TRUE) {
  plot(NA, ylim=c(0,1), xlim=c(min(LL95,low_eqbound,dif)-max(UL95-LL95, high_eqbound-low_eqbound,dif)/10, max(UL95,high_eqbound,dif)+max(UL95-LL95, high_eqbound-low_eqbound, dif)/10), bty="l", yaxt="n", ylab="",xlab="Mean Difference")
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
  # Return results in list(),NHST_t = t, NHST_p = pttest
  invisible(list(diff=dif,TOST_t1=t1,TOST_p1=p1,TOST_t2=t2,TOST_p2=p2, TOST_df=degree_f,alpha=alpha,low_eqbound=low_eqbound,high_eqbound=high_eqbound, LL_CI_TOST=LL90,UL_CI_TOST=UL90, LL_CI_TTEST=LL95, UL_CI_TTEST=UL95,NHST_t = t, NHST_p = pttest))
}
