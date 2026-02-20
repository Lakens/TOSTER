#' @title Test for Association/Correlation Between Paired Samples
#' @description
#' `r lifecycle::badge('stable')`
#'
#' Test for association between paired samples, using one of Pearson's product moment correlation
#' coefficient, Kendall's \eqn{\tau} (tau) or Spearman's \eqn{\rho} (rho). Unlike the stats version of
#' cor.test, this function allows users to set the null to a value other than zero and perform
#' equivalence testing.
#'
#' @param x,y numeric vectors of data values. x and y must have the same length.
#' @param method a character string indicating which correlation coefficient is to be used for the test.
#'   One of "pearson", "kendall", or "spearman", can be abbreviated.
#' @param null a number or vector indicating the null hypothesis value(s):
#'   * For standard tests: a single value (default = 0)
#'   * For equivalence/minimal effect tests: either a single value (symmetric bounds ±value will be created)
#'     or a vector of two values representing the lower and upper bounds
#' @param alternative a character string specifying the alternative hypothesis:
#'   * "two.sided": correlation is not equal to null (default)
#'   * "greater": correlation is greater than null
#'   * "less": correlation is less than null
#'   * "equivalence": correlation is within the equivalence bounds (TOST)
#'   * "minimal.effect": correlation is outside the equivalence bounds (TOST)
#'
#'   You can specify just the initial letter.
#' @param alpha alpha level (default = 0.05)
#' @param se_method a character string indicating the method for computing the standard error.
#'   One of "analytic" (default) or "jackknife". The jackknife SE is computed on the Fisher z scale
#'   using leave-one-out resampling.
#'
#' @details
#' This function uses Fisher's z transformation for the correlations.
#' For Spearman's \eqn{\rho}, the Bonett-Wright \eqn{\rho}-dependent SE formula
#' \eqn{\sqrt{(1 + r^2/2) / (n - 3)}} is used rather than the fixed 1.06 constant.
#' For Kendall's \eqn{\tau}, Fieller's correction is applied.
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
#' * If a single value is provided for `null`, symmetric bounds ± value will be used
#' * If two values are provided for `null`, they will be used as the lower and upper bounds
#'
#' When `se_method = "jackknife"`, the standard error is computed via leave-one-out
#' resampling on the Fisher z scale, which can provide better calibration for small
#' samples or non-standard correlation methods. The jackknife SE is used for both
#' the test statistic and the confidence interval.
#'
#' See `vignette("correlations")` for more details.
#'
#' @return A list with class "htest" containing the following components:
#'
#' * **p.value**: the p-value of the test.
#' * **statistic**: the value of the test statistic with a name describing it.
#' * **parameter**: the degrees of freedom or number of observations.
#' * **conf.int**: a confidence interval for the measure of association appropriate to the specified alternative hypothesis.
#' * **estimate**: the estimated measure of association, with name "cor", "tau", or "rho" corresponding to the method employed.
#' * **stderr**: a named vector with `z.se` (standard error on the Fisher z scale, used for inference)
#'   and `cor.se` (delta method SE on the correlation scale, for descriptive purposes).
#'   Note that `cor.se` underestimates true sampling variability as |r| approaches 1.
#' * **null.value**: the value of the association measure under the null hypothesis.
#' * **alternative**: character string indicating the alternative hypothesis.
#' * **method**: a character string indicating how the association was measured.
#' * **data.name**: a character string giving the names of the data.
#' * **call**: the matched call.
#'
#' @examples
#' # Example 1: Standard significance test
#' x <- c(44.4, 45.9, 41.9, 53.3, 44.7, 44.1, 50.7, 45.2, 60.1)
#' y <- c( 2.6,  3.1,  2.5,  5.0,  3.6,  4.0,  5.2,  2.8,  3.8)
#' z_cor_test(x, y, method = "kendall", alternative = "t", null = 0)
#'
#' # Example 2: Minimal effect test
#' # Testing if correlation is meaningfully different from ±0.2
#' z_cor_test(x, y, method = "kendall", alternative = "min", null = 0.2)
#'
#' # Example 3: Equivalence test with Pearson correlation
#' # Testing if correlation is equivalent to zero within ±0.3
#' z_cor_test(x, y, method = "pearson", alternative = "equivalence", null = 0.3)
#'
#' # Example 4: Using asymmetric bounds
#' # Testing if correlation is within bounds of -0.1 and 0.4
#' z_cor_test(x, y, method = "spearman",
#'            alternative = "equivalence", null = c(-0.1, 0.4))
#'
#' @references
#' Goertzen, J. R., & Cribbie, R. A. (2010). Detecting a lack of association: An equivalence
#' testing approach. British Journal of Mathematical and Statistical Psychology, 63(3), 527-537.
#' https://doi.org/10.1348/000711009X475853, formula page 531.
#'
#' Bonett, D. G., & Wright, T. A. (2000). Sample size requirements for estimating
#' Pearson, Kendall and Spearman correlations. Psychometrika, 65(1), 23-28.
#'
#' @name z_cor_test
#' @family Correlations
#' @export z_cor_test

z_cor_test = function(x,
                      y,
                      alternative = c("two.sided", "less", "greater",
                                      "equivalence", "minimal.effect"),
                      method = c("pearson", "kendall", "spearman"),
                      alpha = 0.05,
                      null = 0,
                      se_method = c("analytic", "jackknife")){
  alternative = match.arg(alternative)
  method = match.arg(method)
  se_method = match.arg(se_method)

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
  r_xy = cor(x,y,
             method = method)
  df = data.frame(x=x,
                  y=y)
  df = na.omit(df)
  n_obs = nrow(df)

  z_xy = rho_to_z(r_xy)
  DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
  NVAL = null
  znull = rho_to_z(null)
  # capture raw method string before overwrite --------
  method_raw <- method
  # get se ---
  if (method == "pearson") {
    # Pearson # Fisher
    method <- "Pearson's product-moment correlation"
    names(NVAL) = rep("correlation",length(NVAL))
    rfinal = c(cor = r_xy)
    z.se <- 1 / sqrt(n_obs - 3)
    cor_correction <- "fieller"
  }
  if (method == "spearman") {
    method <- "Spearman's rank correlation rho"
    # Bonett-Wright
    rfinal = c(rho = r_xy)
    names(NVAL) = rep("rho",length(NVAL))
    z.se <- sqrt((1 + r_xy^2 / 2) / (n_obs - 3))
    cor_correction <- "bw"
  }
  if (method == "kendall") {
    method <- "Kendall's rank correlation tau"
    # Fieller adjusted
    rfinal = c(tau = r_xy)
    names(NVAL) = rep("tau",length(NVAL))
    z.se <- (0.437 / (n_obs - 4)) ^ 0.5
    cor_correction <- "fieller"
  }

  # jackknife SE --------
  if (se_method == "jackknife") {
    jk_cors <- vapply(seq_len(n_obs), function(i) {
      rho_to_z(cor(x[-i], y[-i], method = method_raw))
    }, numeric(1))
    z.se <- sqrt(((n_obs - 1) / n_obs) * sum((jk_cors - mean(jk_cors))^2))
  }

  # append SE label to method name --------
  se_label <- if (se_method == "jackknife") "with jackknifed SE" else "with approximate SE"
  method <- paste(method, se_label)

  # pass se to cor_to_ci only when jackknife --------
  se_override <- if (se_method == "jackknife") z.se else NULL

  cint = cor_to_ci(cor = r_xy, n = n_obs, ci = ci,
                   method = method_raw,
                   correction = cor_correction,
                   se = se_override)

  # delta method SE on correlation scale --------
  cor.se <- (1 - r_xy^2) * z.se

  # get absolute value if TOST
  #z_test = ifelse(TOST, abs(z_xy), z_xy)

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
               stderr = c(z.se = z.se, cor.se = cor.se),
               null.value = NVAL,
               alternative = alternative,
               method = method,
               data.name = DNAME,
               call = match.call())
  class(rval) <- "htest"
  return(rval)
}
