#' @title Bootstrapped correlation coefficients
#' @description A function for a bootstrap, percentile, method for correlation coefficients.
#' @inheritParams boot_t_TOST
#' @inheritParams z_cor_test
#' @details This function uses a percentile bootstrap methods for the confidence intervals.
#' The returned p-values are calculated from a re-sampled null distribution (similar to boot_t_TOST).
#'  @return A list with class "htest" containing the following components:
#' \describe{
#'   \item{\code{"statistic"}}{z-score}
#'   \item{\code{"p.value"}}{the p-value of the test.}
#'   \item{\code{"estimate"}}{the estimated measure of association, with name "cor", "tau", or "rho" corresponding to the method employed.}
#'   \item{\code{"null.value"}}{the value of the association measure under the null hypothesis.}
#'   \item{\code{"alternative"}}{character string indicating the alternative hypothesis (the value of the input argument alternative). Possible values are "greater", "less", or "two-sided".}
#'   \item{\code{"method"}}{a character string indicating how the association was measured.}
#'   \item{\code{"data.name"}}{a character string giving the names of the data.}
#'   \item{\code{"call"}}{the matched call.}
#' }
#' @references
#' TBA
#' @section References:
#'
#' Efron, B., & Tibshirani, R. J. (1994). An introduction to the bootstrap. CRC press.
#' @export


boot_cor_test <- function(x,
                          y,
                          alternative = c("two.sided", "less", "greater"),
                          method = c("pearson", "kendall", "spearman"),
                          alpha = 0.05,
                          null = 0,
                          TOST = FALSE,
                          R = 1999,
                          ...) {
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
  n <- nrow(m)
  x <- df[,1]
  y <- df[,2]
  alternative = match.arg(alternative)
  method = match.arg(method)

  if(TOST && null <=0){
    stop("positive value for null must be supplied if using TOST.")
  }
  if(TOST){
    alternative = "less"
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

  est <- cor(x, y, method)
  data <- matrix(sample(n, size=n*nboot, replace=TRUE), nrow=nboot)
  bvec <- apply(data, 1, .corboot, x, y, method, ...) # Create a 1 by nboot matrix.

  corci = quantile(bvec, c((1 - ci) / 2, 1 - (1 - ci) / 2))

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
  if(saveboot){ # save bootstrap samples
    list(conf.int=corci, p.value=sig, estimate=est, bootsamples=bvec)
  } else {
    list(conf.int=corci, p.value=sig, estimate=est)
  }

  if (method == "pearson") {
    # Pearson # Fisher
    method2 <- "Pearson's product-moment correlation"
    names(null.value) = "correlation"
    rfinal = c(cor = est)
  }
  if (method == "spearman") {
    method2 <- "Spearman's rank correlation rho"
    #  # Fieller adjusted
    rfinal = c(rho = est)
    names(null.value) = "rho"

  }
  if (method == "kendall") {
    method2 <- "Kendall's rank correlation tau"
    # # Fieller adjusted
    rfinal = c(tau = est)
    names(null.value) = "tau"

  }
  DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
  # Store as htest
  rval <- list(p.value = sig,
               conf.int = boot.cint,
               estimate = rfinal,
               stderr = sd(bvec,na.rm=TRUE),
               null.value = null.value,
               alternative = alternative,
               method = method2,
               data.name = DNAME,
               boot = bvec,
               call = match.call())
  class(rval) <- "htest"
  return(rval)
}



