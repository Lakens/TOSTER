#' @title Equivalence Test using an F-test
#' @description
#' `r lifecycle::badge('stable')`
#'
#' Performs equivalence or minimal effect testing on the partial eta-squared (pes) value using
#' an F-test. This function provides a low-level interface that works directly with F statistics
#' rather than ANOVA objects.
#'
#' @param Fstat The F-statistic from the F-test.
#' @param df1 Degrees of freedom for the numerator (effect degrees of freedom).
#' @param df2 Degrees of freedom for the denominator (error degrees of freedom).
#' @param eqbound Equivalence bound for the partial eta-squared. This value represents
#'   the smallest effect size considered meaningful or practically significant.
#' @param eqb Defunct argument for equivalence bound, use `eqbound` instead.
#' @param MET Logical indicator to perform a minimal effect test rather than equivalence
#'   test (default is FALSE). When TRUE, the alternative hypothesis becomes that the effect
#'   is larger than the equivalence bound.
#' @param alpha Alpha level used for the test (default = 0.05).
#'
#' @details
#' This function tests whether an effect is practically equivalent to zero (when
#' `MET = FALSE`) or meaningfully different from zero (when `MET = TRUE`) using the approach
#' described by Campbell & Lakens (2021).
#'
#' The function works by:
#'
#' 1. Converting the F-statistic to a partial eta-squared value
#' 2. Converting the equivalence bound for partial eta-squared to a non-centrality parameter
#' 3. Computing the confidence interval for the partial eta-squared
#' 4. Performing an equivalence test or minimal effect test based on the non-central F distribution
#'
#' For equivalence tests (`MET = FALSE`), a significant result (p < alpha) indicates that the
#' effect is statistically equivalent to zero (smaller than the equivalence bound).
#'
#' For minimal effect tests (`MET = TRUE`), a significant result (p < alpha) indicates that
#' the effect is meaningfully different from zero (larger than the equivalence bound).
#'
#' For details on the calculations in this function see `vignette("the_ftestTOSTER")`.
#'
#' @return
#' Object of class "htest" containing the following components:
#'
#' * **statistic**: The value of the F-statistic with name "F".
#' * **parameter**: The degrees of freedom for the F-statistic (df1 and df2).
#' * **p.value**: The p-value for the equivalence or minimal effect test.
#' * **conf.int**: A confidence interval for the partial eta-squared statistic.
#' * **estimate**: Estimate of partial eta-squared.
#' * **null.value**: The specified equivalence bound.
#' * **alternative**: NULL (not used in this test).
#' * **method**: A string indicating the type of test ("Equivalence Test from F-test" or "Minimal Effect Test from F-test").
#' * **data.name**: A string indicating that this was calculated from summary statistics.
#'
#' @examples
#' # Example 1: Equivalence test with a small effect
#' # F = 2.5, df1 = 2, df2 = 100, equivalence bound = 0.1
#' equ_ftest(Fstat = 2.5, df1 = 2, df2 = 100, eqbound = 0.1)
#'
#' # Example 2: Minimal effect test with a large effect
#' # F = 12, df1 = 3, df2 = 80, equivalence bound = 0.1
#' equ_ftest(Fstat = 12, df1 = 3, df2 = 80, eqbound = 0.1, MET = TRUE)
#'
#' # Example 3: Equivalence test with a very small effect
#' # F = 0.8, df1 = 1, df2 = 50, equivalence bound = 0.05
#' equ_ftest(Fstat = 0.8, df1 = 1, df2 = 50, eqbound = 0.05)
#'
#' @references
#' Campbell, H., & Lakens, D. (2021). Can we disregard the whole model? Omnibus non‐inferiority
#' testing for R2 in multi‐variable linear regression and in ANOVA. British Journal of
#' Mathematical and Statistical Psychology, 74(1), 64-89. doi: 10.1111/bmsp.12201
#'
#' @family f-test
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

