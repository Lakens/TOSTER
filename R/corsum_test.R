#' @title  Association/Correlation Test from Summary Statistics
#' @description   Test for association between paired samples, using one of Pearson's product moment correlation coefficient,Kendall's \eqn{\tau}{tau} or Spearman's \eqn{\rho}{rho}.
#' This is the updated version of the TOSTr function.
#' @inheritParams TOSTr
#' @inheritParams z_cor_test
#' @details
#' This function uses Fisher's z transformation for the correlations, but uses Fieller's correction of the standard error for Spearman's rho and Kendall's tau.
#'
#' @return A list with class "htest" containing the following components:
#'
#' \describe{
#'   \item{\code{"statistic"}}{z-score}
#'   \item{\code{"p.value"}}{the p-value of the test.}
#'   \item{\code{"estimate"}}{the estimated measure of association, with name "cor", "tau", or "rho" corresponding to the method employed.}
#'   \item{\code{"null.value"}}{the value of the association measure under the null hypothesis.}
#'   \item{\code{"alternative"}}{character string indicating the alternative hypothesis (the value of the input argument alternative). }
#'   \item{\code{"method"}}{a character string indicating how the association was measured.}
#'   \item{\code{"data.name"}}{a character string giving the names of the data.}
#'   \item{\code{"call"}}{the matched call.}
#' }
#' @references
#' Goertzen, J. R., & Cribbie, R. A. (2010). Detecting a lack of association: An equivalence testing approach. British Journal of Mathematical and Statistical Psychology, 63(3), 527-537. https://doi.org/10.1348/000711009X475853, formula page 531.
#' @export
#'


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
