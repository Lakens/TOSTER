#' Power analysis for TOST for between subjects ANOVA.
#' @param alpha alpha used for the test (e.g., 0.05)
#' @param power desired power (e.g., 0.8)
#' @param df1 Degrees of freedom for the numerator
#' @param df2 Degrees of freedom for the denominator
#' @param eqbound Equivalence bound in partial eta-squared
#' @return Object of class '"power.htest"
#' @examples
#' ## Statistical power for alpha = 0.05, 3 groups, n = 80 per group, equivalence bound of
#' ## partial eta squared = 0.01, assuming true effect = 0.
#' ## df1 = number of groups - 1 = 3 - 1 = 2.
#' ## df2 = Total N - number of groups = 240 - 3 = 237.
#' #' powerTOSToneway(alpha=0.05, df1=3, df2 = 237, eqbound = 0.01)
#' @section References:
#' TOO BE ADDED
#' @importFrom stats pf qf
#' @export

powerTOSToneway <- function(alpha = 0.05, df1, df2, eqbound){
  N <- df2 + df1 + 1


  Fstatstar = qf(
    alpha,
    df1 = df1,
    df2 = df2,
    ncp = (N * eqbound) / (1 - eqbound),
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
      method = "Equivalence Test for one-way ANOVA"
    ),
    class = "power.htest"
  )

}

