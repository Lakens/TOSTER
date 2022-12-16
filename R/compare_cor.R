#' @title Comparing two independent correlation coefficients
#' @description A function to compare correlations between studies. This function is intended to be used to compare the compatibility of original studies with replication studies (lower p-values indicating lower compatibility).
#' @param r1 Correlation study 1.
#' @param df1 Degrees of freedom from study 1 (if a simple correlation the df is N-2).
#' @param r2 Correlation study 2.
#' @param df2 Degrees of freedom from study 2 (if a simple correlation the df is N-2).
#' @param method Method for determining differences. Default, "z", will use Fisher's transformation, while "AH" will use the Anderson-Hauk approach.
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
#' @references
#' Counsell, A., & Cribbie, R. A. (2015). Equivalence tests for comparing correlation and regression coefficients. The British journal of mathematical and statistical psychology, 68(2), 292–309. https://doi.org/10.1111/bmsp.12045
#'
#' Anderson, S., & Hauck, W. W. (1983). A new procedure for testing equivalence in comparative bioavailability and other clinical trials. Communications in Statistics-Theory and Methods, 12(23), 2663-2692.
#'
#' @name compare_cor
#' @import R6
#' @export compare_cor
#'


compare_cor = function(r1,
                       df1,
                       r2,
                       df2,
                       method = c("fisher","ah","kraatz"),
                       alternative = c("two.sided", "less", "greater"),
                       null = 0,
                       TOST = FALSE){
  method = match.arg(method)
  alternative <- match.arg(alternative)

  if(TOST && null <=0){
    stop("positive value for null must be supplied if using TOST.")
  }
  if(TOST){
    alternative = "less"
  }

  if(!TOST && method == "ah"){
    stop("ah method can only be used")
  }

  if(method == "ah"){
    meth2 = "(Anderson-Hauk)"
  }
  if(method == "fisher"){
    meth2 = "(Fisher's z transform)"
  }
  if(method == "kraatz"){
    meth2 = "(Kraatz)"
  }
  # z transform and SE
  if(method == "fisher"){
    z1 = (1/2) * log((1 + r1)/(1 - r1))
    z2 = (1/2) * log((1 + r2)/(1 - r2))
    znull = (1/2) * log((1 + null)/(1 - null))
    if(TOST){
      z_diff = abs(z1 - z2) - znull
    } else{
      z_diff = z1 - z2 - znull
    }
    z_se = sqrt(1/(df1-1) + 1/(df2-1))
    z = z_diff/z_se

  }
  # Anderson-Hauk Method
  if(method == "ah" || method == "kraatz"){
    se = sqrt((1-r1^2)^2/(df1)+(1-r2^2)^2/(df2))
    if(TOST){
      diff = abs(r1 - r2) - null
    } else{
      diff = r1-r2-null
    }
    z = diff/se
  }
  if(method != "ah"){
    pval = p_from_z(z, alternative = alternative)
    names(z) = "z"
    est2 = r1-r2
    names(est2) = "difference between correlations"
    null2 = null
    names(null2) = "difference between correlations"
    meth = "Difference between two independent correlations"
    meth_final = paste0(meth," ",meth2)
    # Store as htest
    rval <- list(statistic = z, p.value = pval,
                 #conf.int = cint,
                 estimate = est2,
                 null.value = null2,
                 alternative = alternative,
                 method = meth_final,
                 cor = list(r1 = r1,
                            r2 = r2),
                 data.name = "Summary Statistics",
                 call = match.call())
  } else {
    diff = r1-r2-null
    pval <- pnorm((abs(diff) - null)/se) -
      pnorm((-abs(diff) - null)/se)
    est2 = r1-r2
    names(est2) = "difference between correlations"
    null2 = null
    names(null2) = "difference between correlations"
    meth = "Difference between two independent correlations"

    meth_final = paste0(meth," ",meth2)
    # Store as htest
    rval <- list(p.value = pval,
                 #conf.int = cint,
                 estimate = est2,
                 null.value = null2,
                 alternative = alternative,
                 method = meth_final,
                 cor = list(r1 = r1,
                            r2 = r2),
                 data.name = "Summary Statistics",
                 call = match.call())
  }


  class(rval) <- "htest"
  return(rval)


}
