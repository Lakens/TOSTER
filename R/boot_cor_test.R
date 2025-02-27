#' @title Bootstrapped Correlation Coefficients
#' @description
#' `r lifecycle::badge('stable')`
#'
#' A function for bootstrap-based correlation tests using various correlation coefficients
#' including Pearson's, Kendall's, Spearman's, Winsorized, and percentage bend correlations.
#' This function supports standard, equivalence, and minimal effect testing with robust bootstrap methods.
#'
#' @inheritParams boot_t_TOST
#' @inheritParams z_cor_test
#' @param method a character string indicating which correlation coefficient to use:
#'   * "pearson": standard Pearson product-moment correlation
#'   * "kendall": Kendall's tau rank correlation
#'   * "spearman": Spearman's rho rank correlation
#'   * "winsorized": Winsorized correlation (robust to outliers)
#'   * "bendpercent": percentage bend correlation (robust to marginal outliers)
#'
#'   Can be abbreviated.
#' @param boot_ci type of bootstrap confidence interval:
#'   * "basic": basic/empirical bootstrap CI
#'   * "perc": percentile bootstrap CI (default)
#' @param R number of bootstrap replications (default = 1999).
#' @param ... additional arguments passed to correlation functions, such as:
#'   * trim: for Winsorized correlation (default = 0.2)
#'   * beta: for percentage bend correlation (default = 0.2)
#'
#' @details
#' This function uses bootstrap methods to calculate correlation coefficients and their
#' confidence intervals. P-values are calculated from a re-sampled null distribution.
#'
#' The bootstrap correlation methods in this package offer two robust correlations beyond
#' the standard methods:
#'
#' 1. **Winsorized correlation**: Replaces extreme values with less extreme values before
#'    calculating the correlation. The `trim` parameter (default = 0.2) determines the
#'    proportion of data to be Winsorized.
#'
#' 2. **Percentage bend correlation**: A robust correlation that downweights the influence
#'    of outliers. The `beta` parameter (default = 0.2) determines the bending constant.
#'
#' These calculations are based on Rand Wilcox's R functions for his book (Wilcox, 2017),
#' and adapted from their implementation in Guillaume Rousselet's R package "bootcorci".
#'
#' The function supports both standard hypothesis testing and equivalence/minimal effect testing:
#'
#' * For standard tests (two.sided, less, greater), the function tests whether the correlation
#'   differs from the null value (typically 0).
#'
#' * For equivalence testing ("equivalence"), it determines whether the correlation falls within
#'   the specified bounds, which can be set asymmetrically.
#'
#' * For minimal effect testing ("minimal.effect"), it determines whether the correlation falls
#'   outside the specified bounds.
#'
#' When performing equivalence or minimal effect testing:
#' * If a single value is provided for `null`, symmetric bounds ±value will be used
#' * If two values are provided for `null`, they will be used as the lower and upper bounds
#'
#' See `vignette("correlations")` for more details.
#'
#' @return A list with class "htest" containing the following components:
#'
#' * **p.value**: the bootstrap p-value of the test.
#' * **parameter**: the number of observations used in the test.
#' * **conf.int**: a bootstrap confidence interval for the correlation coefficient.
#' * **estimate**: the estimated correlation coefficient, with name "cor", "tau", "rho", "pb", or "wincor"
#'   corresponding to the method employed.
#' * **stderr**: the bootstrap standard error of the correlation coefficient.
#' * **null.value**: the value(s) of the correlation under the null hypothesis.
#' * **alternative**: character string indicating the alternative hypothesis.
#' * **method**: a character string indicating which bootstrapped correlation was measured.
#' * **data.name**: a character string giving the names of the data.
#' * **boot_res**: vector of bootstrap correlation estimates.
#' * **call**: the matched call.
#'
#' @examples
#' # Example 1: Standard bootstrap test with Pearson correlation
#' x <- c(44.4, 45.9, 41.9, 53.3, 44.7, 44.1, 50.7, 45.2, 60.1)
#' y <- c( 2.6,  3.1,  2.5,  5.0,  3.6,  4.0,  5.2,  2.8,  3.8)
#' boot_cor_test(x, y, method = "pearson", alternative = "two.sided",
#'               R = 999) # Fewer replicates for example
#'
#' # Example 2: Equivalence test with Spearman correlation
#' # Testing if correlation is equivalent to zero within ±0.3
#' boot_cor_test(x, y, method = "spearman", alternative = "equivalence",
#'              null = 0.3, R = 999)
#'
#' # Example 3: Using robust correlation methods
#' # Using Winsorized correlation with custom trim
#' boot_cor_test(x, y, method = "winsorized", trim = 0.1,
#'              R = 999)
#'
#' # Example 4: Using percentage bend correlation
#' boot_cor_test(x, y, method = "bendpercent", beta = 0.2,
#'              R = 999)
#'
#' # Example 5: Minimal effect test with asymmetric bounds
#' # Testing if correlation is outside bounds of -0.1 and 0.4
#' boot_cor_test(x, y, method = "pearson", alternative = "minimal.effect",
#'              null = c(-0.1, 0.4), R = 999)
#'
#' @section References:
#' Wilcox, R.R. (2009) Comparing Pearson Correlations: Dealing with Heteroscedasticity and Nonnormality.
#' Communications in Statistics - Simulation and Computation, 38, 2220–2234.
#'
#' Wilcox, R.R. (2017) Introduction to Robust Estimation and Hypothesis Testing, 4th edition. Academic Press.
#'
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
                     "basic" = basic(bvec, t0 = est, alpha2),
                     "perc" = perc(bvec, alpha2))
  #quantile(bvec, c((1 - ci) / 2, 1 - (1 - ci) / 2))
  attr(boot.cint, "conf.level") <- ci
  # pvalue
  ## note method different than Efron (e.g., t-test, SMDs, etc)
  ## Dervied from work of Wilcox
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



