#' @title Power Calculations for Correlations
#' @description
#' `r lifecycle::badge('maturing')`
#'
#'  Calculates the approximate power for a z-test based on a correlation.
#'
#' @param n number of observations.
#' @param rho true correlation value (alternative hypothsis).
#' @param null the null hypothesis value.
#' @param alpha a priori alpha-level (i.e., significance level).
#' @param power statistical power (1-beta).
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater", "less", or "equivalence" (TOST). You can specify just the initial letter.
#'
#' @references
#' Lakens, D. (2017). Equivalence Tests: A Practical Primer for t Tests, Correlations, and Meta-Analyses. Social Psychological and Personality Science, 8(4), 355–362. https://doi.org/10.1177/1948550617697177
#'
#' Cohen, J. (1988). Statistical power analysis for the behavioral sciences (2nd Ed). Hillsdale, NJ: Lawrence Erlbaum Associates.
#'
#' Zhang, Z., & Yuan, K.-H. (2018). Practical Statistical Power Analysis Using Webpower and R (Eds). Granger, IN: ISDSA Press.
#'
#' @importFrom stats integrate uniroot qnorm pnorm dnorm
#' @family power
#' @export
#'

power_z_cor = function(n = NULL,
                       rho = NULL,
                       power = NULL,
                       null = 0,
                       alpha = NULL,
                       alternative = c("two.sided", "less", "greater", "equivalence")){

  alternative = match.arg(alternative)

  if(alternative == "equivalence"){
    if(length(null) == 1){
      if(null ==  0){
        stop("null cannot be zero if alternative is equivalence")
      }

      null = c(null,-1 * null)

    }
    pow_corr_tost(n = n, r = rho, power = power, null = null,
                  alpha = alpha)
  } else{
    pow_corr(n = null, r = r, power = power, null = null,
             alpha = alpha, alternative = alternative)
  }
}
