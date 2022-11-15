#' @title Equivalence Test using an F-test
#' @description Performs equivalence test on the partial eta-squared (pes) value for using an F-test.
#' @param Fstat The F-statistic from the F-test.
#' @param df1 Degrees of freedom for the numerator.
#' @param df2 Degrees of freedom for the denominator.
#' @param eqb Defunct argument for quivalence bound for the partial eta-squared.
#' @param eqbound Defunct argument for quivalence bound for the partial eta-squared. Default is NULL.
#' @param MET logical indicator to perform a minimal effect test rather than equivalence test (default is FALSE).
#' @param alpha alpha used for the test (e.g., 0.05).
#' @details For details on the calculations in this function see vignette("the_ftestTOSTER").
#' @return Object of class '"htest"
#' \describe{
#'   \item{\code{"statistic"}}{The value of the F-statistic.}
#'   \item{\code{"parameter"}}{The degrees of freedom for the F-statistic.}
#'   \item{\code{"p.value"}}{The he p-value for the test.}
#'   \item{\code{"conf.int"}}{A confidence interval for the partial eta-squared statistic.}
#'   \item{\code{"estimate"}}{Estimate of partial eta-squared.}
#'   \item{\code{"null.value"}}{The specified for the equivalence test.}
#'   \item{\code{"method"}}{A string indicating the type of F-test.}
#'   \item{\code{"data.name"}}{A required string indicating that this was calculated from summary statistics.}
#' }
#' @family f-test
#' @section References:
#' Campbell, H., & Lakens, D. (2021). Can we disregard the whole model? Omnibus non‐inferiority testing for R2 in multi‐variable linear regression and in ANOVA. British Journal of Mathematical and Statistical Psychology, 74(1), 64-89. doi: 10.1111/bmsp.12201
#' @importFrom stats pf qf
#' @export

equ_ftest <- function(Fstat,
                      df1,
                      df2,
                      eqbound = NULL,
                      eqb,
                      MET = FALSE,
                      alpha = 0.05){
  #message("Note: equ_ftest only validated for one-way ANOVA; use with caution")
  if(!is.null(eqbound)){
    eqb = eqbound
  }
  eqbound = eqb
  conf_level = 1 - alpha

  pes = Fstat * df1 / (Fstat*df1+df2)

  F_limits <- conf.limits.ncf(F.value = Fstat,
                              df.1 = df1,
                              df.2 = df2,
                              conf.level = conf_level)
  LL_lambda <- F_limits$Lower.Limit
  UL_lambda <- F_limits$Upper.Limit

  LL_partial_eta2 <- LL_lambda / (LL_lambda + df1 + df2 + 1)
  UL_partial_eta2 <- UL_lambda / (UL_lambda + df1 + df2 + 1)


  if (is.na(LL_partial_eta2)) {
    LL_partial_eta2 <- 0
  }

  if (is.na(UL_partial_eta2)) {
    UL_partial_eta2 <- 1
  }

  f2 = eqbound/(1 - eqbound)
  lambda = (f2 * (df1 + df2 + 1))

  l_tail = ifelse(MET,FALSE,TRUE)

  if(MET){
    method = "Minimal Effect Test from F-test"
  } else {
    method = "Equivalence Test from F-test"
  }

  pval = pf(Fstat,
            df1 = df1,
            df2 = df2,
            ncp = lambda,
            lower.tail = l_tail)
  cint = c(LL_partial_eta2 ,UL_partial_eta2)
  attr(cint,"conf.level") <- conf_level
  parameters = c(df1, df2)
  names(parameters) <- c("df1","df2")

  names(Fstat) = "F"

  alt_val = paste0("partial eta-squared less than ", eqbound)

  rval <- list(statistic = Fstat,
               parameter = parameters,
               p.value = pval,
               conf.int = cint,
               estimate = pes,
               null.value = eqbound,
               alternative = NULL,
               method = method,
               data.name = "Summary Statistics")
  class(rval) <- "htest"
  return(rval)
}

