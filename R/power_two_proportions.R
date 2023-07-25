#' @name power_twoprop
#' @aliases powerTOSTtwo.prop
#' @aliases power_twoprop
#' @title TOST Power for Tests of Two Proportions
#'
#' @description
#' `r lifecycle::badge("maturing")`
#'
#' Power analysis for TOST for difference between two proportions using Z-test (pooled)
#' @param n Sample size per group.
#' @inheritParams power_z_cor
#' @inheritParams twoprop_test
#' @param prop1 Deprecated. expected proportion in group 1.
#' @param prop2 Deprecated. expected proportion in group 2.
#' @param statistical_power Deprecated. desired power (e.g., 0.8)
#' @param N Deprecated. sample size (e.g., 108)
#' @param low_eqbound_prop Deprecated. lower equivalence bounds (e.g., -0.05) expressed in proportion
#' @param high_eqbound_prop Deprecated. upper equivalence bounds (e.g., 0.05) expressed in proportion
#' @return Calculate either achieved power, equivalence bounds, or required N, assuming a true effect size of 0.
#' Returns a string summarizing the power analysis, and a numeric variable for number of observations, equivalence bounds, or power.
#' @examples
#' ## Sample size for alpha = 0.05, 90% power, assuming true effect prop1 = prop 2 = 0.5,
#' ## equivalence bounds of 0.4 and 0.6 (so low_eqbound_prop = -0.1 and high_eqbound_prop = 0.1)
#'
#' powerTOSTtwo.prop(alpha = 0.05, statistical_power = 0.9, prop1 = 0.5, prop2 = 0.5,
#'    low_eqbound_prop = -0.1, high_eqbound_prop = 0.1)
#'
#' ## Power for alpha = 0.05, N 542 , assuming true effect prop1 = prop 2 = 0.5,
#' ## equivalence bounds of 0.4 and 0.6 (so low_eqbound_prop = -0.1 and high_eqbound_prop = 0.1)
#'
#' powerTOSTtwo.prop(alpha = 0.05, N = 542, prop1 = 0.5, prop2 = 0.5,
#'    low_eqbound_prop = -0.1, high_eqbound_prop = 0.1)
#'
#' ## Equivalence bounds for alpha = 0.05, N 542 , assuming true effect prop1 = prop 2 = 0.5,
#' ## and 90% power
#'
#' powerTOSTtwo.prop(alpha=0.05, statistical_power=0.9, N=542, prop1 = 0.5, prop2 = 0.5)
#'
#' #Example 4.2.4 from Chow, Wang, & Shao (2007, p. 93)
#' powerTOSTtwo.prop(alpha=0.05, statistical_power=0.8, prop1 = 0.75, prop2 = 0.8,
#'    low_eqbound_prop = -0.2, high_eqbound_prop = 0.2)
#'
#' # Example 5 from Julious & Campbell (2012, p. 2932)
#' powerTOSTtwo.prop(alpha=0.025, statistical_power=0.9, prop1 = 0.8, prop2 = 0.8,
#'    low_eqbound_prop=-0.1, high_eqbound_prop=0.1)
#' # From Machin, D. (Ed.). (2008). Sample size tables for clinical studies (3rd ed).
#'
#' # Example 9.4b equivalence of two proportions (p. 113) #
#' powerTOSTtwo.prop(alpha=0.010, statistical_power=0.8, prop1 = 0.5, prop2 = 0.5,
#'    low_eqbound_prop = -0.2, high_eqbound_prop = 0.2)/2
#' @references
#' Silva, G. T. da, Logan, B. R., & Klein, J. P. (2008). Methods for Equivalence and Noninferiority Testing. Biology of Blood and Marrow Transplantation: Journal of the American Society for Blood and Marrow Transplantation, 15(1 Suppl), 120-127. https://doi.org/10.1016/j.bbmt.2008.10.004
#'
#' Julious, S. A. & Campell, M. J. (2012). Tutorial in biostatistics: sample sizes for parallel group clinical trials with binary data. Statistics in Medicine, 31:2904-2936.
#'
#' Chow, S.-C., Wang, H., & Shao, J. (2007). Sample Size Calculations in Clinical Research, Second Edition (2 edition). Boca Raton: Chapman and Hall/CRC.
#' @importFrom stats pnorm pt qnorm qt
#' @export

powerTOSTtwo.prop <- function(alpha,
                              statistical_power,
                              prop1,
                              prop2,
                              N,
                              low_eqbound_prop,
                              high_eqbound_prop)
{
    if(missing(N)) {
    NT1<-(prop1*(1-prop1)+prop2*(1-prop2))*((qnorm(1-alpha)+qnorm(1-((1-statistical_power)/2)))/(abs(prop1-prop2)-abs(low_eqbound_prop)))^2
    NT2<-(prop1*(1-prop1)+prop2*(1-prop2))*((qnorm(1-alpha)+qnorm(1-((1-statistical_power)/2)))/(abs(prop1-prop2)-abs(high_eqbound_prop)))^2
    N<-max(NT1,NT2)
    message(cat("The required sample size to achieve",100*statistical_power,"% power with equivalence bounds of",low_eqbound_prop,"and",high_eqbound_prop,"is",ceiling(N)))
    return(N)
  }
  if(missing(statistical_power)) {
    statistical_power1<-2*(pnorm((abs(prop1-prop2)-low_eqbound_prop)/sqrt(prop1*(1-prop1)/N+prop2*(1-prop2)/N)-qnorm(1-alpha))+pnorm(-(abs(prop1-prop2)-low_eqbound_prop)/sqrt(prop1*(1-prop1)/N+prop2*(1-prop2)/N)-qnorm(1-alpha)))-1
    statistical_power2<-2*(pnorm((abs(prop1-prop2)-high_eqbound_prop)/sqrt(prop1*(1-prop1)/N+prop2*(1-prop2)/N)-qnorm(1-alpha))+pnorm(-(abs(prop1-prop2)-high_eqbound_prop)/sqrt(prop1*(1-prop1)/N+prop2*(1-prop2)/N)-qnorm(1-alpha)))-1
    statistical_power<-min(statistical_power1,statistical_power2)
    if(statistical_power<0) {statistical_power<-0}
    message(cat("The statistical power is",round(100*statistical_power,2),"% for equivalence bounds of",low_eqbound_prop,"and",high_eqbound_prop,"."))
    return(statistical_power)
  }
  if(missing(low_eqbound_prop) && missing(high_eqbound_prop)) {
    low_eqbound_prop<--sqrt((prop1*(1-prop1)+prop2*(1-prop2))*((qnorm(1-alpha)+qnorm(1-((1-statistical_power)/2))))^2/N)
    high_eqbound_prop<-sqrt((prop1*(1-prop1)+prop2*(1-prop2))*((qnorm(1-alpha)+qnorm(1-((1-statistical_power)/2))))^2/N)
    message(cat("The equivalence bounds to achieve",100*statistical_power,"% power with N =",N,"are",round(low_eqbound_prop,2),"and",round(high_eqbound_prop,2),"."))
    bounds<-c(low_eqbound_prop,high_eqbound_prop)
    return(bounds)
  }
}

#' @rdname power_twoprop
#' @export

power_twoprop = function(p1, p2,
                         n = NULL,
                         null = 0,
                         alpha = NULL,
                         power = NULL,
                         alternative = c("two.sided",
                                         "less",
                                         "greater",
                                         "equivalence")){
  if(missing(p1) || missing(p2)){
    stop("proportions (p1 & p2) must be supplied")
  }

  alternative = match.arg(alternative)

  if(alternative == "equivalence"){
    if(length(null) == 1){
      if(null ==  0){
        stop("null cannot be zero if alternative is equivalence")
      }

      null = c(null,-1 * null)

    }
    pow_prop_tost(p1, p2,
                  n,
                  null,
                  alpha,
                  power)
  } else{
    pow_prop(p1, p2,
             n,
             null,
             alpha,
             power,
             alternative)
  }



}


pow_prop = function (p1, p2,
                     n = NULL,
                     null = 0,
                     alpha = NULL,
                     power = NULL,
                     alternative = c("two.sided",
                                     "less",
                                     "greater"))
{
  if (sum(sapply(list(n, power, alpha), is.null)) != 1)
    stop("exactly one of n, power, and alpha must be NULL")
  if (!is.null(n) && min(n) < 1)
    stop("number of observations in each group must be at least 1")
  if (!is.null(alpha) && !is.numeric(alpha) || any(0 > alpha | alpha >
                                                   1))
    stop(sQuote("alpha"), " must be numeric in [0, 1]")
  if (!is.null(power) && !is.numeric(power) || any(0 > power | power >
                                                   1))
    stop(sQuote("power"), " must be numeric in [0, 1]")
  alternative <- match.arg(alternative)
  tside <- switch(alternative, less = 1, two.sided = 2, greater = 3)

  if (tside == 3) {
    p.body <- quote({
      prop_se <- sqrt((p1*(1-p1))/n + (p2*(1-p2))/n)
      prop_dif <- p1 - p2 - null
      zval = qnorm(alpha)
      pnorm(qnorm(alpha, lower = FALSE) - (zval*prop_se+prop_dif)/prop_se, lower = FALSE)
    })
  }
  if (tside == 2) {

    p.body <- quote({
      prop_se <- sqrt((p1*(1-p1))/n + (p2*(1-p2))/n)
      prop_dif <- p1 - p2 - null
      zval = qnorm(alpha/2)
      pnorm(qnorm(alpha/2, lower = FALSE) - (zval*prop_se+prop_dif)/prop_se, lower = FALSE) +
        pnorm(qnorm(alpha/2, lower = TRUE) - (zval*prop_se+prop_dif)/prop_se, lower = TRUE)
    })
  }
  if (tside == 1) {
    p.body <- quote({
      prop_se <- sqrt((p1*(1-p1))/n + (p2*(1-p2))/n)
      prop_dif <- p1 - p2 - null
      zval = qnorm(alpha)
      pnorm(qnorm(alpha, lower = TRUE) - (zval*prop_se+prop_dif)/prop_se, lower = TRUE)
    })
  }

  if (is.null(power))
    power <- eval(p.body)
  else if (is.null(n))
    n <- uniroot(function(n) eval(p.body) - power, c(4 + 1e-10, 1e+07))$root
  else if (is.null(alpha))
    alpha <- uniroot(function(alpha) eval(p.body) - power,
                     c(1e-10, 1 - 1e-10))$root
  else if (is.null(alpha))
    null <- uniroot(function(null) eval(p.body) - power,
                     c(1e-10, 1 - 1e-10))$root
   else stop("internal error")


  METHOD <- "Power for Test of Differences in Two Proportions (z-test)"
  NOTE = "Sample sizes for EACH group"
  structure(list(n = n,
                 proportions = c(p1,p2),
                 alpha = alpha,
                 beta = 1-power, power = power,
                 null = null, alternative = alternative,
                 method = METHOD,
                 NOTE = NOTE),
            class = "power.htest")
}

pow_prop_tost = function (p1, p2,
                          n = NULL,
                          null = NULL,
                          alpha = NULL,
                          power = NULL)
{

  if (sum(sapply(list(n, power, alpha), is.null)) != 1)
    stop("exactly one of n, power, and alpha must be NULL")
  if(is.null(r)){
    stop("r cannot be set to NULL at this time.")
  }
  if (!is.null(alpha) && !is.numeric(alpha) || any(0 > alpha |
                                                   alpha > 1))
    stop(sQuote("alpha"), " must be numeric in [0, 1]")
  if (!is.null(power) && !is.numeric(power) || any(0 > power |
                                                   power > 1))
    stop(sQuote("power"), " must be numeric in [0, 1]")
  if (!is.null(n) && min(n) < 4)
    stop("number of observations must be at least 4")

  alternative <- "equivalence"

  p.body =  quote({
    statistical_power1<-2*(pnorm((abs(prop1-prop2)-min(null))/sqrt(prop1*(1-prop1)/N+prop2*(1-prop2)/N)-qnorm(1-alpha))+pnorm(-(abs(prop1-prop2)-min(null))/sqrt(prop1*(1-prop1)/N+prop2*(1-prop2)/N)-qnorm(1-alpha)))-1
    statistical_power2<-2*(pnorm((abs(prop1-prop2)-max(null))/sqrt(prop1*(1-prop1)/N+prop2*(1-prop2)/N)-qnorm(1-alpha))+pnorm(-(abs(prop1-prop2)-max(null))/sqrt(prop1*(1-prop1)/N+prop2*(1-prop2)/N)-qnorm(1-alpha)))-1
    statistical_power<-min(statistical_power1,statistical_power2)
    if(statistical_power<0) {statistical_power<-0}
    statistical_power
  })

  if (is.null(power))
    power <- eval(p.body)
  else if (is.null(n))
    n <- uniroot(function(n) eval(p.body) - power, c(4 + 1e-10, 1e+07))$root
  else if (is.null(alpha))
    alpha <- uniroot(function(alpha) eval(p.body) - power,
                     c(1e-10, 1 - 1e-10))$root
  else stop("internal error")
  METHOD <- "Power for Test of Differences in Two Proportions (z-test)"
  NOTE = "Sample sizes for EACH group"
  structure(list(n = n, rho = r,
                 alpha = alpha, beta = 1-power, power = power,
                 null = null, alternative = alternative,
                 method = METHOD,
                 NOTE = NOTE),
            class = "power.htest")
}
