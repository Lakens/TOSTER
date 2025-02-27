#' @title Comparing Correlations Between Independent Studies with Bootstrapping
#' @description
#' `r lifecycle::badge('stable')`
#'
#' A function to compare correlation coefficients between independent studies using bootstrap methods.
#' This function is intended to be used to compare the compatibility of original studies
#' with replication studies (lower p-values indicating lower compatibility).
#'
#' @param x1,y1 Numeric vectors of data values from study 1. x1 and y1 must have the same length.
#' @param x2,y2 Numeric vectors of data values from study 2. x2 and y2 must have the same length.
#' @param ... Additional arguments passed to the correlation functions.
#' @inheritParams boot_cor_test
#'
#' @details
#' This function tests for differences between correlation coefficients from independent studies
#' using bootstrap resampling methods. Unlike the `compare_cor` function, which uses Fisher's z
#' transformation or the Kraatz method with summary statistics, this function works with raw
#' data and uses bootstrapping to estimate confidence intervals and p-values.
#'
#' It is particularly useful for:
#'
#' * Comparing correlations when assumptions for parametric tests may not be met
#' * Obtaining robust confidence intervals for the difference between correlations
#' * Comparing an original study with its replication using raw data
#' * Testing if correlations from different samples are equivalent
#'
#' The function supports multiple correlation methods:
#'
#' * Standard correlation coefficients (Pearson, Kendall, Spearman)
#' * Robust correlation measures (Winsorized, percentage bend)
#'
#' The function also supports both standard hypothesis testing and equivalence/minimal effect testing:
#'
#' * For standard tests (two.sided, less, greater), the function tests whether the difference
#'   between correlations differs from the null value (typically 0).
#'
#' * For equivalence testing ("equivalence"), it determines whether the difference falls within
#'   the specified bounds, which can be set asymmetrically.
#'
#' * For minimal effect testing ("minimal.effect"), it determines whether the difference falls
#'   outside the specified bounds.
#'
#' When performing equivalence or minimal effect testing:
#' * If a single value is provided for `null`, symmetric bounds ±value will be used
#' * If two values are provided for `null`, they will be used as the lower and upper bounds
#'
#' @return A list with class "htest" containing the following components:
#'
#' * **p.value**: The p-value for the test under the null hypothesis.
#' * **parameter**: Sample sizes from each study.
#' * **conf.int**: Bootstrap confidence interval for the difference in correlations.
#' * **estimate**: Difference in correlations between studies.
#' * **stderr**: Standard error of the difference (estimated from bootstrap distribution).
#' * **null.value**: The specified hypothesized value(s) for the null hypothesis.
#' * **alternative**: Character string indicating the alternative hypothesis.
#' * **method**: Description of the correlation method used.
#' * **data.name**: Names of the input data vectors.
#' * **boot_res**: List containing the bootstrap samples for the difference and individual correlations.
#' * **call**: The matched call.
#'
#' @examples
#' # Example 1: Comparing Pearson correlations (standard test)
#' set.seed(123)
#' x1 <- rnorm(30)
#' y1 <- x1 * 0.6 + rnorm(30, 0, 0.8)
#' x2 <- rnorm(25)
#' y2 <- x2 * 0.3 + rnorm(25, 0, 0.9)
#'
#' # Two-sided test with Pearson correlation (use fewer bootstraps for example)
#' boot_compare_cor(x1, y1, x2, y2, method = "pearson",
#'                 alternative = "two.sided", R = 500)
#'
#' # Example 2: Testing for equivalence with Spearman correlation
#' # Testing if the difference in correlations is within ±0.2
#' boot_compare_cor(x1, y1, x2, y2, method = "spearman",
#'                 alternative = "equivalence", null = 0.2, R = 500)
#'
#' # Example 3: Testing with robust correlation measure
#' # Using percentage bend correlation for non-normal data
#' boot_compare_cor(x1, y1, x2, y2, method = "bendpercent",
#'                 alternative = "greater", R = 500)
#'
#' # Example 4: Using asymmetric bounds for equivalence testing
#' boot_compare_cor(x1, y1, x2, y2, method = "pearson",
#'                 alternative = "equivalence", null = c(-0.1, 0.3), R = 500)
#'
#' @name boot_compare_cor
#' @family compare studies
#' @importFrom stats median
#' @export boot_compare_cor
#'

boot_compare_cor <- function(
    x1,y1,
    x2,y2,
    alternative = c("two.sided","less", "greater",
                    "equivalence", "minimal.effect"),
    method = c("pearson", "kendall", "spearman",
               "winsorized", "bendpercent"),
    alpha = 0.05,
    null = 0,
    R = 1999,
    ...){
  DNAME <- paste(deparse(substitute(x1)), "and", deparse(substitute(y1)),
                 "vs.",
                 deparse(substitute(x2)), "and", deparse(substitute(y2)))
  nboot = R
  null.value = null
  if(!is.vector(x1) || !is.vector(x2) || !is.vector(y1) || !is.vector(y2)){
    stop("x1, x2, y1, y2  must be vectors.")
  }
  if(length(x1)!=length(y1)){
    stop("x1 and y1 do not have equal lengths.")
  }
  if(length(x2)!=length(y2)){
    stop("x2 and y2 do not have equal lengths.")
  }

  #if(TOST && null <=0){
  #  stop("positive value for null must be supplied if using TOST.")
  #}
  #if(TOST){
  #  alternative = "less"
  #}

  if(alternative %in% c("equivalence", "minimal.effect")){
    if(length(null) == 1){
      null = c(null, -1*null)
    }
    TOST = TRUE
  } else {
    if(length(null) > 1){
      stop("null can only have 1 value for non-TOST procedures")
    }
    TOST = FALSE
  }

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

  alternative = match.arg(alternative)
  method = match.arg(method)

  df <- cbind(x1,y1)
  df <- df[complete.cases(df), ]
  n1 <- nrow(df)
  x1 <- df[,1]
  y1 <- df[,2]
  df <- cbind(x2,y2)
  df <- df[complete.cases(df), ]
  n2 <- nrow(df)
  x2 <- df[,1]
  y2 <- df[,2]

  if(method %in% c("bendpercent","winsorized")){
    if(method == "bendpercent"){
      r1 <- pbcor(x1, y1, ...)
      r2 <- pbcor(x2, y2, ...)
      # bootstrap
      data1 <- matrix(sample(n1, size=n1*nboot, replace=TRUE), nrow=nboot)
      bvec1 <- apply(data1, 1, .corboot_pbcor, x1, y1, ...) # A 1 by nboot matrix.
      data2 <- matrix(sample(n2, size=n2*nboot, replace=TRUE), nrow=nboot)
      bvec2 <- apply(data2, 1, .corboot_pbcor, x2, y2, ...) # A 1 by nboot matrix.
    }

    if(method == "winsorized"){
      r1 <- wincor(x1, y1, ...)
      r2 <- wincor(x2, y2, ...)

      # bootstrap
      data1 <- matrix(sample(n1, size=n1*nboot, replace=TRUE), nrow=nboot)
      bvec1 <- apply(data1, 1, .corboot_wincor, x1, y1, ...) # A 1 by nboot matrix.
      data2 <- matrix(sample(n2, size=n2*nboot, replace=TRUE), nrow=nboot)
      bvec2 <- apply(data2, 1, .corboot_wincor, x2, y2, ...) # A 1 by nboot matrix
    }


  } else {
    # correlations
    r1 <- cor(x1,y1,method = method)
    r2 <- cor(x2,y2,method = method)
    # bootstrap
    data1 <- matrix(sample(n1, size=n1*nboot, replace=TRUE), nrow=nboot)
    bvec1 <- apply(data1, 1, .corboot, x1, y1,  method = method) # A 1 by nboot matrix.
    data2 <- matrix(sample(n2, size=n2*nboot, replace=TRUE), nrow=nboot)
    bvec2 <- apply(data2, 1, .corboot, x2, y2, method = method) # A 1 by nboot matrix.
  }


  bvec <- bvec1 - bvec2
  bsort <- sort(bvec)
  boot.cint = quantile(bvec, c((1 - ci) / 2, 1 - (1 - ci) / 2))
  attr(boot.cint, "conf.level") <- ci
  # p value
  ## note method different than Efron (e.g., t-test, SMDs, etc)
  ## Derived from work of Wilcox
  if(alternative == "two.sided"){
    phat <- (sum(bvec < null.value)+.5*sum(bvec==0))/nboot
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
    method2 <- "Bootstrapped difference in Pearson's correlation"
    names(null.value) = "difference in correlation"
    rfinal = c(cor = r1-r2)
  }
  if (method == "spearman") {
    method2 <- "Bootstrapped difference in Spearman's rho"
    #  # Fieller adjusted
    rfinal = c(rho = r1-r2)
    names(null.value) = "difference in rho"
  }
  if (method == "kendall") {
    method2 <- "Bootstrapped difference in Kendall's tau"
    # # Fieller adjusted
    rfinal = c(tau = r1-r2)
    names(null.value) = "difference in tau"
  }
  if (method == "bendpercent") {
    method2 <- "Bootstrapped difference in percentage bend correlation pb"
    # # Fieller adjusted
    rfinal = c(pb = r1-r2)
    names(null.value) = "difference in pb"
  }
  if (method == "winsorized") {
    method2 <- "Bootstrapped difference in Winsorized correlation wincor"
    # # Fieller adjusted
    rfinal = c(wincor = r1-r2)
    names(null.value) = "differnce in wincor"
  }
  N = c(n1 = n1, n2 = n2)
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
               boot_res = list(diff = bvec,
                           r1 = bvec1,
                           r2 = bvec2),
               call = match.call())
  class(rval) <- "htest"
  return(rval)
}
