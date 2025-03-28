#' @title Association/Correlation Test from Summary Statistics
#' @description
#' `r lifecycle::badge('stable')`
#'
#' Test for association between paired samples using only the correlation coefficient and sample size.
#' Supports Pearson's product moment correlation, Kendall's \eqn{\tau} (tau), or Spearman's \eqn{\rho} (rho).
#' This is the updated version of the `TOSTr` function.
#'
#' @param r correlation coefficient (the estimated value)
#' @param n sample size (number of pairs)
#' @inheritParams TOSTr
#' @inheritParams z_cor_test
#'
#' @details
#' This function uses Fisher's z transformation for the correlations,
#' but uses Fieller's correction of the standard error for Kendall's \eqn{\tau} or Spearman's \eqn{\rho}.
#'
#' Unlike `z_cor_test`, which requires raw data, this function only needs the correlation value
#' and sample size. This is particularly useful when:
#'
#' * You only have access to summary statistics (correlation coefficient and sample size)
#' * You want to reanalyze published results within an equivalence testing framework
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
#' @return A list with class "htest" containing the following components:
#'
#' * **statistic**: z-score with name "z".
#' * **p.value**: the p-value of the test.
#' * **parameter**: the sample size with name "N".
#' * **conf.int**: a confidence interval for the correlation appropriate to the specified alternative hypothesis.
#' * **estimate**: the estimated correlation coefficient, with name "cor", "tau", or "rho" corresponding to the method employed.
#' * **stderr**: the standard error of the test statistic.
#' * **null.value**: the value(s) of the correlation coefficient under the null hypothesis.
#' * **alternative**: character string indicating the alternative hypothesis.
#' * **method**: a character string indicating how the correlation was measured.
#' * **data.name**: a character string giving the names of the data.
#' * **call**: the matched call.
#'
#' @examples
#' # Example 1: Standard significance test for Pearson correlation
#' corsum_test(r = 0.45, n = 30, method = "pearson", alternative = "two.sided")
#'
#' # Example 2: Equivalence test for Spearman correlation
#' # Testing if correlation is equivalent to zero within ±0.3
#' corsum_test(r = 0.15, n = 40, method = "spearman",
#'             alternative = "equivalence", null = 0.3)
#'
#' # Example 3: Minimal effect test for Kendall's tau
#' # Testing if correlation is meaningfully different from ±0.25
#' corsum_test(r = 0.42, n = 50, method = "kendall",
#'             alternative = "minimal.effect", null = 0.25)
#'
#' # Example 4: One-sided test with non-zero null
#' # Testing if correlation is greater than 0.3
#' corsum_test(r = 0.45, n = 35, method = "pearson",
#'             alternative = "greater", null = 0.3)
#'
#' # Example 5: Using asymmetric bounds for equivalence testing
#' corsum_test(r = 0.1, n = 60, method = "pearson",
#'             alternative = "equivalence", null = c(-0.2, 0.3))
#'
#' @references
#' Goertzen, J. R., & Cribbie, R. A. (2010). Detecting a lack of association: An equivalence
#' testing approach. British Journal of Mathematical and Statistical Psychology, 63(3), 527-537.
#' https://doi.org/10.1348/000711009X475853, formula page 531.
#'
#' @family Correlations
#' @export


corsum_test = function(r,
                       n,
                       alternative = c("two.sided", "less", "greater",
                                       "equivalence", "minimal.effect"),
                       method = c("pearson", "kendall", "spearman"),
                       alpha = 0.05,
                       null = 0){
  alternative = match.arg(alternative)
  method = match.arg(method)

  #if(TOST && null <=0){
  #  stop("positive value for null must be supplied if using TOST.")
  #}
  #if(TOST){
  #  alternative = "less"
  #}
  #if(TOST && alternative %in% c("two.sided","greater","less")){
  #  alternative = "equivalence"
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
  r_xy = r
  n_obs = n

  z_xy = rho_to_z(r_xy)
  DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))

  # get absolute value if TOST
  #z_test = ifelse(TOST, abs(z_xy), z_xy)

  NVAL = null
  znull = rho_to_z(null)
  # get se ---
  if (method == "pearson") {
    # Pearson # Fisher
    method <- "Pearson's product-moment correlation"
    names(NVAL) = rep("correlation",length(NVAL))
    rfinal = c(cor = r_xy)
    z.se <- 1 / sqrt(n_obs - 3)
    cint = cor_to_ci(cor = r_xy, n = n_obs, ci = ci,
                     method = "pearson")
  }
  if (method == "spearman") {
    method <- "Spearman's rank correlation rho"
    #  # Fieller adjusted
    rfinal = c(rho = r_xy)
    names(NVAL) = rep("rho",length(NVAL))
    z.se <- (1.06 / (n_obs - 3)) ^ 0.5
    cint = cor_to_ci(cor = r_xy, n = n_obs, ci = ci,
                     method = "spearman",
                     correction = "fieller")
  }
  if (method == "kendall") {
    method <- "Kendall's rank correlation tau"
    # # Fieller adjusted
    rfinal = c(tau = r_xy)
    names(NVAL) = rep("tau",length(NVAL))
    z.se <- (0.437 / (n_obs - 4)) ^ 0.5

    cint = cor_to_ci(cor = r_xy, n = n_obs, ci = ci,
                     method = "kendall",
                     correction = "fieller")
  }
  if(alternative %in% c("equivalence", "minimal.effect")){
    if(alternative == "equivalence"){
      zlo = z_xy-min(znull)
      plo = p_from_z(zlo/z.se, alternative = 'greater')
      zhi = z_xy-max(znull)
      phi = p_from_z(zhi/z.se, alternative = 'less')
      if(phi >= plo){
        pvalue = phi
        z_test = zhi
      } else {
        pvalue = plo
        z_test = zlo
      }
    }
    if(alternative == "minimal.effect"){
      zlo = z_xy-min(znull)
      plo = p_from_z(zlo/z.se, alternative = 'less')
      zhi = z_xy-max(znull)
      phi = p_from_z(zhi/z.se, alternative = 'greater')
      if(phi <= plo){
        pvalue = phi
        z_test = zhi
      } else {
        pvalue = plo
        z_test = zlo
      }
    }

  } else{

    z_test = z_xy-znull
    pvalue = p_from_z(z_test/z.se, alternative = alternative)
  }

  z_test2 = z_test/z.se
  names(z_test2) = "z"
  attr(cint, "conf.level") <- ci
  N = n_obs
  names(N) = "N"

  # Store as htest
  rval <- list(statistic = z_test2, p.value = pvalue,
               parameter = N,
               conf.int = cint,
               estimate = rfinal,
               stderr = c(z.se = z.se),
               null.value = NVAL,
               alternative = alternative,
               method = method,
               data.name = DNAME,
               call = match.call())
  class(rval) <- "htest"
  return(rval)
}
