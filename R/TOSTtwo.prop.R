#' TOST function for two proportions (raw scores)
#'
#' #' @description
#' `r lifecycle::badge('superseded')`
#'
#' Development on `TOSTtwo.prop` is complete, and for new code we recommend
#' switching to [tsum_TOST], which is easier to use, more featureful,
#' and still under active development.
#'
#' @param prop1 proportion of group 1
#' @param prop2 proportion of group 2
#' @param n1 sample size in group 1
#' @param n2 sample size in group 2
#' @param low_eqbound lower equivalence bounds (e.g., -0.1) expressed in proportions
#' @param high_eqbound upper equivalence bounds (e.g., 0.1) expressed in proportions
#' @param alpha alpha level (default = 0.05)
#' @param ci_type confidence interval type (default = "normal"). "wilson" produces Wilson score intervals with a Yates continuity correction while "normal" calculates the simple asymptotic method with no continuity correction.
#' @param plot set whether results should be plotted (plot = TRUE) or not (plot = FALSE) - defaults to TRUE
#' @param verbose logical variable indicating whether text output should be generated (verbose = TRUE) or not (verbose = FALSE) - default to TRUE
#' @return Returns TOST z-value 1, TOST p-value 1, TOST z-value 2, TOST p-value 2, low equivalence bound, high equivalence bound, Lower limit confidence interval TOST, Upper limit confidence interval TOST
#' @importFrom stats pnorm pt qnorm qt
#' @importFrom graphics abline plot points segments title
#' @examples
#' ## Equivalence test for two independent proportions equal to .65 and .70, with 100 samples
#' ## per group, lower equivalence bound of -0.1, higher equivalence bound of 0.1, and alpha of 0.05.
#'
#' TOSTtwo.prop(prop1 = .65, prop2 = .70, n1 = 100, n2 = 100,
#'    low_eqbound = -0.1, high_eqbound = 0.1, alpha = .05)
#'
#' @references
#' Tunes da Silva, G., Logan, B. R., & Klein, J. P. (2008). Methods for Equivalence and Noninferiority Testing. Biology of Blood Marrow Transplant, 15(1 Suppl), 120-127.
#'
#' Yin, G. (2012). Clinical Trial Design: Bayesian and Frequentist Adaptive Methods. Hoboken, New Jersey: John Wiley & Sons, Inc.
#' @export

TOSTtwo.prop <- function(prop1, prop2, n1, n2, low_eqbound, high_eqbound, alpha, ci_type = "normal", plot = TRUE, verbose = TRUE) {
  if(missing(alpha)) {
    alpha <- 0.05
  }
  if(low_eqbound >= high_eqbound) warning("The lower bound is equal to or larger than the upper bound. Check the plot and output to see if the bounds are specified as you intended.")
  if(n1 < 2 | n2 < 2) stop("The sample size should be larger than 1.")
  if(1 <= alpha | alpha <= 0) stop("The alpha level should be a positive value between 0 and 1.")
  if(1 < prop1 | prop1 < 0) stop("The proportion should be a positive value between 0 and 1.")
  if(1 < prop2 | prop2 < 0) stop("The proportion should be a positive value between 0 and 1.")
  if(ci_type != "normal" & ci_type != "wilson") stop("ci_type must be equal to wilson or normal")

    prop_dif <- prop1 - prop2
  # proportion se
  prop_se <- sqrt((prop1*(1-prop1))/n1 + (prop2*(1-prop2))/n2)

  #calculating z-statistic
  z1 <- (prop_dif - low_eqbound)/prop_se
  z2 <- (prop_dif - high_eqbound)/prop_se
  z  <- prop_dif / prop_se
  ztest <- 1 - pnorm(abs(z))

  #calculating p-value for both one-sided tests
  p1 <- 1 - pnorm(z1)
  p2 <- pnorm(z2)
  ptost <- max(p1,p2) #Get highest p-value for summary TOST result
  ztost <- ifelse(abs(z1) < abs(z2), z1, z2) #Get lowest z-value for summary TOST result
  TOSToutcome <- ifelse(ptost<alpha,"significant","non-significant")
  testoutcome <- ifelse(ztest<(alpha/2), "significant","non-significant")
  # Wilson CI method
  conf.level = 1 - alpha
  DELTA = prop_dif
  ESTIMATE = c(p1, p2)
  YATES <- abs(DELTA) / sum(1 / n1, 1 / n2)
  n = c(n1,n2)
  WIDTH90 <- qnorm(conf.level) * sqrt(sum(ESTIMATE * (1 - ESTIMATE) / n)) + YATES * sum(1 / n)
  WIDTH95 <- (qnorm((1+conf.level)/2))* sqrt(sum(ESTIMATE * (1 - ESTIMATE) / n)) + YATES * sum(1 / n)
  CINT90 <- c(max(DELTA - WIDTH90, -1), min(DELTA + WIDTH90, 1))
  CINT95 <- c(max(DELTA - WIDTH95, -1), min(DELTA + WIDTH95, 1))

  #calculating CIs
  if(ci_type == "normal"){
    LL90 <- prop_dif - (qnorm(1-alpha) * prop_se)
    UL90 <- prop_dif + (qnorm(1-alpha) * prop_se)
    LL95 <- prop_dif - (qnorm(1-(alpha/2)) * prop_se)
    UL95 <- prop_dif + (qnorm(1-(alpha/2)) * prop_se)
  }

  if(ci_type == "wilson"){
    LL90 <- CINT90[1]
    UL90 <- CINT90[2]
    LL95 <- CINT95[1]
    UL95 <- CINT95[2]
  }


  #plot results
  if (plot == TRUE) {
  plot(NA, ylim=c(0,1), xlim=c(min(LL90,(low_eqbound))-max(UL90-LL90, high_eqbound-(low_eqbound))/10, max(UL90,high_eqbound)+max(UL90-LL90, high_eqbound-(low_eqbound))/10), bty="l", yaxt="n", ylab="",xlab="Proportion Difference")
  points(x=prop_dif, y=0.5, pch=15, cex=2)
  abline(v=high_eqbound, lty=2)
  abline(v=low_eqbound, lty=2)
  abline(v=0, lty=2, col="grey")
  segments(LL90,0.5,UL90,0.5, lwd=3)
  segments(LL95,0.5,UL95,0.5, lwd=1)
  title(main=paste("Equivalence bounds ",round(low_eqbound,digits=3)," and ",round(high_eqbound,digits=3),"\nProportion Difference = ",round(prop_dif,3)," \n TOST: ", 100*(1-alpha*2),"% CI [",round(LL90,digits=3),";",round(UL90,digits=3),"] ", TOSToutcome, " \n NHST: ", 100*(1-alpha),"% CI [",round(LL95,digits=3),";",round(UL95,digits=3),"] ", testoutcome, sep=""), cex.main=1)
  }


  if(missing(verbose)) {
    verbose <- TRUE
  }
  if(verbose == TRUE){
    cat("TOST results:\n")
    cat("Z-value lower bound:",format(z1, digits = 3, nsmall = 2, scientific = FALSE),"\tp-value lower bound:",format(p1, digits = 1, nsmall = 3, scientific = FALSE))
    cat("\n")
    cat("Z-value upper bound:",format(z2, digits = 3, nsmall = 2, scientific = FALSE),"\tp-value upper bound:",format(p2, digits = 1, nsmall = 3, scientific = FALSE))
    cat("\n\n")
    cat("Equivalence bounds:")
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
    cat("Equivalence Test based on Fisher's exact z-test Result:\n")
    message(cat("The equivalence test was ",TOSToutcome,", Z = ",format(ztost, digits = 3, nsmall = 3, scientific = FALSE),", p = ",format(ptost, digits = 3, nsmall = 3, scientific = FALSE),", given equivalence bounds of ",format(low_eqbound, digits = 3, nsmall = 3, scientific = FALSE)," and ",format(high_eqbound, digits = 3, nsmall = 3, scientific = FALSE)," and an alpha of ",alpha,".",sep=""))
    cat("\n")
    cat("Null-Hypothesis Fisher's exact z-test Result:\n")
    message(cat("The null hypothesis test was ",testoutcome,", Z = ",format(z, digits = 3, nsmall = 3, scientific = FALSE),", p = ",format((ztest * 2), digits = 3, nsmall = 3, scientific = FALSE),", given an alpha of ",alpha,".",sep=""))
    if((ztest * 2) <= alpha && ptost <= alpha){
      combined_outcome <- paste0("NHST: reject null significance hypothesis\n",
                                 "TOST: reject null equivalence hypothesis")
      #"statistically different from zero and statistically equivalent to zero"
    }
    if((ztest * 2) < alpha && ptost > alpha){
      combined_outcome <- paste0("NHST: reject null significance hypothesis \n",
                                 "TOST: don't reject null equivalence hypothesis")
      #"statistically different from zero and statistically not equivalent to zero"
    }
    if((ztest * 2) > alpha && ptost <= alpha){
      combined_outcome <- paste0("NHST: don't reject null significance hypothesis that the effect is equal to ",0," \n",
                                 "TOST: reject null equivalence hypothesis")
      #"statistically not different from zero and statistically equivalent to zero"
    }
    if((ztest * 2) > alpha && ptost > alpha){
      combined_outcome <- paste0("NHST: don't reject null significance hypothesis that the effect is equal to ",0," \n",
                                 "TOST: don't reject null equivalence hypothesis")#"statistically not different from zero and statistically not equivalent to zero"
    }
    cat("\n")
    message(combined_outcome)
  }
  # Return results in list()
  invisible(list(dif = prop_dif, TOST_z1 = z1, TOST_p1 = p1,TOST_z2 = z2,TOST_p2 = p2, alpha = alpha, low_eqbound = low_eqbound, high_eqbound = high_eqbound, LL_CI_TOST = LL90, UL_CI_TOST = UL90, LL_CI_ZTEST = LL95, UL_CI_ZTEST = UL95, NHST_z = z, NHST_p = (ztest * 2)))
}
