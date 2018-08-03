#' TOST function for a dependent t-test (Cohen's dz)
#' @param m1 mean of group 1
#' @param m2 mean of group 2
#' @param sd1 standard deviation of group 1
#' @param sd2 standard deviation of group 2
#' @param n sample size (pairs)
#' @param r12 correlation of dependent variable between group 1 and  group 2
#' @param low_eqbound_dz lower equivalence bounds (e.g., -0.5) expressed in standardized mean difference (Cohen's dz)
#' @param high_eqbound_dz upper equivalence bounds (e.g., 0.5) expressed in standardized mean difference (Cohen's dz)
#' @param alpha alpha level (default = 0.05)
#' @param plot set whether results should be plotted (plot = TRUE) or not (plot = FALSE) - defaults to TRUE
#' @param verbose logical variable indicating whether text output should be generated (verbose = TRUE) or not (verbose = FALSE) - default to TRUE
#' @return Returns TOST t-value 1, TOST p-value 1, TOST t-value 2, TOST p-value 2, degrees of freedom, low equivalence bound, high equivalence bound, low equivalence bound in dz, high equivalence bound in dz, Lower limit confidence interval TOST, Upper limit confidence interval TOST
#' @examples
#' ## Test means of 5.83 and 5.75, standard deviations of 1.17 and 1.29 in sample of 65 pairs
#' ## with correlation between observations of 0.75 using equivalence bounds in Cohen's dz of
#' ## -0.4 and 0.4 (with default alpha setting of = 0.05).
#' TOSTpaired(n=65,m1=5.83,m2=5.75,sd1=1.17,sd2=1.29,r12=0.75,low_eqbound_dz=-0.4,high_eqbound_dz=0.4)
#' @section References:
#' Mara, C. A., & Cribbie, R. A. (2012). Paired-Samples Tests of Equivalence. Communications in Statistics - Simulation and Computation, 41(10), 1928-1943. https://doi.org/10.1080/03610918.2011.626545, formula page 1932. Note there is a typo in the formula: n-1 should be n (personal communication, 31-8-2016)
#' @importFrom stats pnorm pt qnorm qt
#' @importFrom graphics abline plot points segments title
#' @export
#'

TOSTpaired<-function(n,m1,m2,sd1,sd2,r12,low_eqbound_dz, high_eqbound_dz, alpha, plot = TRUE, verbose = TRUE){
  if(missing(alpha)) {
    alpha <- 0.05
  }
  if(low_eqbound_dz >= high_eqbound_dz) warning("The lower bound is equal to or larger than the upper bound. Check the plot and output to see if the bounds are specified as you intended.")
  if(n < 2) stop("The sample size should be larger than 1.")
  if(1 <= alpha | alpha <= 0) stop("The alpha level should be a positive value between 0 and 1.")
  if(sd1 <= 0 | sd2 <= 0) stop("The standard deviation should be a positive value.")
  # Calculate TOST, t-test, 90% CIs and 95% CIs
  sdif<-sqrt(sd1^2+sd2^2-2*r12*sd1*sd2)
  low_eqbound<-low_eqbound_dz*sdif
  high_eqbound<-high_eqbound_dz*sdif
  se<-sdif/sqrt(n)
  t<-(m1-m2)/se
  degree_f<-n-1
  pttest<-2*pt(abs(t), degree_f, lower.tail=FALSE)
  t1<-((m1-m2)-(low_eqbound_dz*sdif))/se
  p1<-pt(t1, degree_f, lower.tail=FALSE)
  t2<-((m1-m2)-(high_eqbound_dz*sdif))/se
  p2<-pt(t2, degree_f, lower.tail=TRUE)
  ttost<-ifelse(abs(t1) < abs(t2), t1, t2)
  LL90<-((m1-m2)-qt(1-alpha, degree_f)*se)
  UL90<-((m1-m2)+qt(1-alpha, degree_f)*se)
  ptost<-max(p1,p2)
  results<-data.frame(t1,p1,t2,p2,degree_f,LL90,UL90)
  colnames(results) <- c("t-value 1","p-value 1","t-value 2","p-value 2","df", paste("Lower Limit ",100*(1-alpha*2),"% CI",sep=""),paste("Upper Limit ",100*(1-alpha*2),"% CI",sep=""))
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
    cat("Equivalence bounds (Cohen's dz):")
    cat("\n")
    cat("low eqbound:", paste0(round(low_eqbound_dz, digits = 4)),"\nhigh eqbound:",paste0(round(high_eqbound_dz, digits = 4)))
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
      combined_outcome <- "statistically different from zero and statistically equivalent to zero"
    }
    if(pttest < alpha && ptost > alpha){
      combined_outcome <- "statistically different from zero and statistically not equivalent to zero"
    }
    if(pttest > alpha && ptost <= alpha){
      combined_outcome <- "statistically not different from zero and statistically equivalent to zero"
    }
    if(pttest > alpha && ptost > alpha){
      combined_outcome <- "statistically not different from zero and statistically not equivalent to zero"
    }
    cat("\n")
    message(cat("Based on the equivalence test and the null-hypothesis test combined, we can conclude that the observed effect is ",combined_outcome,".",sep=""))
  }
  # Return results in list()
  invisible(list(diff=dif,TOST_t1=t1,TOST_p1=p1,TOST_t2=t2,TOST_p2=p2, TOST_df=degree_f,alpha=alpha,low_eqbound=low_eqbound,high_eqbound=high_eqbound,low_eqbound_dz=low_eqbound_dz,high_eqbound_dz=high_eqbound_dz, LL_CI_TOST=LL90,UL_CI_TOST=UL90, LL_CI_TTEST=LL95, UL_CI_TTEST=UL95))
}
