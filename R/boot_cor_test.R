#' @title Bootstrapped correlation coefficients
#' @description
#' `r lifecycle::badge('stable')`
#'
#'  A function for a bootstrap, percentile, method for correlation coefficients.
#' @inheritParams boot_t_TOST
#' @inheritParams z_cor_test
#' @param method a character string indicating which correlation coefficient is to be used for the test. One of "winsorized", "bendpercent","pearson", "kendall", or "spearman", can be abbreviated.
#' @param boot_ci type of bootstrap confidence interval. Options include studentized (stud), empirical/basic (basic) and percentile (perc) confidence intervals.
#' @details This function uses a percentile bootstrap methods for the confidence intervals.
#' The returned p-values are calculated from a re-sampled null distribution (similar to [boot_t_TOST]).
#' See `vignette("correlations")` for more details.
#'
#' The bootstrap correlation methods in this package offer two other correlations: a Winsorized correlation and a percentage bend correlation (see Wilcox 2017).
#' These two can modified by adding the trim (Winsorized) or beta (percentage bend) arguments.
#' The default for both arguments is 0.2 but can be modified at the user's discretion.
#' These calculations are based on Rand Wilcox's R functions for his book (Wilcox, 2017), and adapted from their implementation in Guillaume Rousselet's R package "bootcorci".
#'
#' @return A list with class "htest" containing the following components:
#'
#'   - "p.value": the p-value of the test.
#'   - "estimate": the estimated measure of association, with name "pb", "wincor", "cor", "tau", or "rho" corresponding to the method employed.
#'   - "null.value": the value of the association measure under the null hypothesis.
#'   - "alternative": character string indicating the alternative hypothesis (the value of the input argument alternative).
#'   - "method": a character string indicating how the association was measured.
#'   - "data.name": a character string giving the names of the data.
#'   - "call": the matched call.
#'
#' @section References:
#'
#' Wilcox, R.R. (2009) Comparing Pearson Correlations: Dealing with Heteroscedasticity and Nonnormality.
#' Communications in Statistics - Simulation and Computation, 38, 2220–2234.
#'
#' Wilcox, R.R. (2017) Introduction to Robust Estimation and Hypothesis Testing, 4th edition. Academic Press.
#' @family Correlations
#' @export


boot_cor_test <- function(x,
                          y,
                          alternative = c("two.sided", "less", "greater",
                                          "equivalence", "minimal.effect"),
                          method = c("pearson", "kendall", "spearman",
                                     "winsorized", "bendpercent"),
                          alpha = 0.05,
                          null = 0,
                          boot_ci = c("basic","perc"),
                          R = 1999,
                          ...) {
  boot_ci = match.arg(boot_ci)
  DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
  alternative = match.arg(alternative)

  method = match.arg(method)
  nboot = R
  null.value = null
  if(!is.vector(x) || !is.vector(y)){
    stop("x and y must be vectors.")
  }
  if(length(x)!=length(y)){
    stop("the vectors do not have equal lengths.")
  }
  df <- cbind(x,y)
  df <- df[complete.cases(df), ]
  n <- nrow(df)
  x <- df[,1]
  y <- df[,2]

  if(alternative %in% c("equivalence", "minimal.effect")){
    if(length(null) == 1){
      null.value = c(null.value, -1*null.value)
    }
    TOST = TRUE
  } else {
    if(length(null) > 1){
      stop("null can only have 1 value for non-TOST procedures")
    }
    TOST = FALSE
  }

  #if(TOST && null <=0){
  #  stop("positive value for null must be supplied if using TOST.")
  #}
  #if(TOST){
  #  alternative = "less"
  #}

  if(alternative != "two.sided"){
    ci = 1 - alpha*2
    intmult = c(1,1)
  } else {
    ci = 1 - alpha
    if(TOST){
      intmult = c(1,1)
    } else if(alternative == "less"){
      intmult = c(1,NA)
    } else {
      intmult = c(NA,1)
    }
  }

  if(method %in% c("bendpercent","winsorized")){
    if(method == "bendpercent"){
      est <- pbcor(x, y, ...)
      data <- matrix(sample(n, size=n*nboot, replace=TRUE), nrow=nboot)
      bvec <- apply(data, 1, .corboot_pbcor, x, y, ...) # get bootstrap results corr
    }

    if(method == "winsorized"){
      est <- wincor(x, y, ...)
      data <- matrix(sample(n, size=n*nboot, replace=TRUE), nrow=nboot)
      bvec <- apply(data, 1, .corboot_wincor, x, y, ...) # get bootstrap results corr
    }


  } else {
    est <- cor(x, y, method = method)
    data <- matrix(sample(n, size=n*nboot, replace=TRUE), nrow=nboot)
    bvec <- apply(data, 1, .corboot, x, y, method = method, ...) # get bootstrap results corr
  }
  alpha2 = ifelse(alternative != "two.sided",
                  alpha*2,
                  alpha)
  boot.cint = switch(boot_ci,
                     "basic" = basic(m_vec, t0 = est, alpha2),
                     "perc" = perc(bvec, alpha2))
  #quantile(bvec, c((1 - ci) / 2, 1 - (1 - ci) / 2))
  attr(boot.cint, "conf.level") <- ci
  if(alternative == "two.sided"){
    phat <- (sum(bvec < null.value)+.5*sum(bvec==null.value))/nboot
    sig <- 2 * min(phat, 1 - phat)
  }
  if(alternative == "greater"){
    sig <- 1 - sum(bvec >= null.value)/nboot
  }
  if(alternative == "less"){
    sig <- 1 - sum(bvec <= null.value)/nboot
  }
  if(alternative == "equivalence"){
    #sig2 <- 1 - sum(bvec >= -1*null.value)/nboot
    #sig = max(sig,sig2)
    sig1 = 1 - sum(bvec >= min(null.value))/nboot
    sig2 = 1 - sum(bvec <= max(null.value))/nboot
    sig = max(sig1,sig2)
  }
  if(alternative == "minimal.effect"){
    #sig2 <- 1 - sum(bvec >= -1*null.value)/nboot
    #sig = max(sig,sig2)
    sig1 = 1 - sum(bvec >= max(null.value))/nboot
    sig2 = 1 - sum(bvec <= min(null.value))/nboot
    sig = min(sig1,sig2)
  }


  if (method == "pearson") {
    # Pearson # Fisher
    method2 <- "Bootstrapped Pearson's product-moment correlation"
    names(null.value) = rep("correlation",length(null.value))
    rfinal = c(cor = est)
  }
  if (method == "spearman") {
    method2 <- "Bootstrapped Spearman's rank correlation rho"
    #  # Fieller adjusted
    rfinal = c(rho = est)
    names(null.value) = rep("rho",length(null.value))

  }
  if (method == "kendall") {
    method2 <- "Bootstrapped Kendall's rank correlation tau"
    # # Fieller adjusted
    rfinal = c(tau = est)
    names(null.value) = rep("tau",length(null.value))

  }
  if (method == "bendpercent") {
    method2 <- "Bootstrapped percentage bend correlation pb"
    # # Fieller adjusted
    rfinal = c(pb = est)
    names(null.value) = rep("pb",length(null.value))

  }
  if (method == "winsorized") {
    method2 <- "Bootstrapped Winsorized correlation wincor"
    # # Fieller adjusted
    rfinal = c(wincor = est)
    names(null.value) = rep("wincor",length(null.value))

  }
  N = n
  names(N) = "N"
  # Store as htest
  rval <- list(p.value = sig,
               parameter = N,
               conf.int = boot.cint,
               estimate = rfinal,
               stderr = sd(bvec,na.rm=TRUE),
               null.value = null.value,
               alternative = alternative,
               method = method2,
               data.name = DNAME,
               boot_res = bvec,
               call = match.call())
  class(rval) <- "htest"
  return(rval)
}



