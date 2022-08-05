#' @title Comparing two independent correlation coefficients
#' @description A function to compare correlations between studies. This function is intended to be used to compare the compatibility of original studies with replication studies (lower p-values indicating lower compatibility).
#' @param r1 Correlation study 1.
#' @param df1 Degrees of freedom from study 1 (if a simple correlation the df is N-2).
#' @param r2 Correlation study 2.
#' @param df2 Degrees of freedom from study 2 (if a simple correlation the df is N-2).
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less". You can specify just the initial letter.
#'  @return A list with class "htest" containing the following components:
#' \describe{
#'   \item{\code{"statistic"}}{z-score}
#'   \item{\code{"p.value"}}{numeric scalar containing the p-value for the test under the null hypothesis.}
#'   \item{\code{"estimate"}}{difference in SMD between studies}
#'   \item{\code{"null.value"}}{the specified hypothesized value for the null hypothesis.}
#'   \item{\code{"alternative"}}{character string indicating the alternative hypothesis (the value of the input argument alternative). Possible values are "greater", "less", or "two-sided".}
#'   \item{\code{"method"}}{Type of SMD}
#'   \item{\code{"data.name"}}{"Summary Statistics" to denote summary statistics were utilized to obtain results.}
#'   \item{\code{"cor"}}{Correlation input for the function.}
#'   \item{\code{"call"}}{the matched call.}
#' }
#' @name compare_cor
#' @export compare_cor
#'


compare_cor = function(r1,
                       df1,
                       r2,
                       df2,
                       alternative = c("two.sided", "less", "greater")){

  alternative <- match.arg(alternative)

  # z transform and SE

  z1 = (1/2) * log((1 + r1)/(1 - r1))
  z2 = (1/2) * log((1 + r2)/(1 - r2))
  z_diff = z1 - z2
  z_se = sqrt(1/(df1-1) + 1/(df2-1))
  z = z_diff/z_se


  pval = p_from_z(z, alternative = alternative)
  names(z) = "z"
  est2 = r1-r2
  names(est2) = "difference between correlations"
  null2 = 0
  names(null2) = "difference between correlations"
  meth = "Difference between two independent correlations"
  # Store as htest
  rval <- list(statistic = z, p.value = pval,
               #conf.int = cint,
               estimate = est2,
               null.value = null2,
               alternative = alternative,
               method = meth,
               cor = list(r1 = r1,
                          r2 = r2),
               data.name = "Summary Statistics",
               call = match.call())
  class(rval) <- "htest"
  return(rval)


}
