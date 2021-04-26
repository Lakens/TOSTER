#' Power analysis for TOST for an F-test
#' @param alpha alpha used for the test (e.g., 0.05)
#' @param df1 Degrees of freedom for the numerator
#' @param df2 Degrees of freedom for the denominator
#' @param eqbound Equivalence bound for the partial eta-squared
#' @return Object of class '"power.htest"
#' @examples
#' ## Statistical power for alpha = 0.05, 3 groups, n = 80 per group, equivalence bound of
#' ## partial eta squared = 0.01, assuming true effect = 0.
#' ## df1 = number of groups - 1 = 3 - 1 = 2.
#' ## df2 = Total N - number of groups = 240 - 3 = 237.
#' #' powerTOSToneway(alpha=0.05, df1=3, df2 = 237, eqbound = 0.01)
#' @section References:
#' Campbell, H., & Lakens, D. (2021). Can we disregard the whole model? Omnibus non‐inferiority testing for R2 in multi‐variable linear regression and in ANOVA. British Journal of Mathematical and Statistical Psychology, 74(1), 64-89. doi: 10.1111/bmsp.12201
#' @importFrom stats pf qf
#' @export

powerTOST_f <- function(alpha = 0.05,
                        df1,
                        df2,
                        eqbound){
  f2 = eqbound/(1 - eqbound)
  lambda = (f2 * (df1 + df2 + 1))


  Fstatstar = qf(
    alpha,
    df1 = df1,
    df2 = df2,
    ncp = lambda,
    lower.tail = TRUE
  )

  statistical_power = pf(Fstatstar,
                         df1 = df1,
                         df2 = df2,
                         lower.tail = TRUE)

  # return h-test

  structure(
    list(
      df1 = df1,
      df2 = df2,
      eqbound = eqbound,
      sig.level = alpha,
      power = statistical_power,
      method = "Equivalence Test for an F-test"
    ),
    class = "power.htest"
  )

}

