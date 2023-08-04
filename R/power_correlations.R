#' @name power_z_cor
#' @aliases powerTOSTr
#' @aliases power_z_cor
#' @title Power Calculations for Correlations
#'
#' @description
#' `r lifecycle::badge("maturing")`
#'
#' Calculates the approximate power for a z-test based on a Pearson product-moment correlation.
#'
#' @param n number of observations.
#' @param rho true correlation value (alternative hypothsis).
#' @param null the null hypothesis value.
#' @param alpha a priori alpha-level (i.e., significance level).
#' @param power statistical power (1-beta).
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater", "less", or "equivalence" (TOST). You can specify just the initial letter.
#' @param statistical_power Deprecated. desired power (e.g., 0.8)
#' @param N Deprecated. number of pairs (e.g., 96)
#' @param low_eqbound_r Deprecated. lower equivalence bounds (e.g., -0.3) expressed in a correlation effect size
#' @param high_eqbound_r Deprecated. upper equivalence bounds (e.g., 0.3) expressed in a correlation effect size
#'
#'
#' @return An object of the class power.htest.
#' This will include the sample size (n), power, beta (1-power), alpha (significance level), null value(s), alternative hypothesis, and a text string detailing the method.
#'
#' `powerTOSTr` has been replaced by the `power_z_cor` function. The function is only retained for historical purposes.
#'
#' @examples
#' ## Sample size for alpha = 0.05, 90% power, equivalence bounds of
#' ## r = -0.1 and r = 0.1, assuming true effect = 0
#' #powerTOSTr(alpha=0.05, statistical_power=0.9, low_eqbound_r=-0.1, high_eqbound_r=0.1)
#' power_z_cor(alternative = "equivalence", alpha = .05, null = .1, power = .9, rho = 0)
#'
#' ## Sample size for alpha = 0.05, N=536, equivalence bounds of
#' ## r = -0.1 and r = 0.1, assuming true effect = 0
#' #powerTOSTr(alpha=0.05, N=536, low_eqbound_r=-0.1, high_eqbound_r=0.1)
#' power_z_cor(alternative = "equivalence", alpha = .05, null = .1, n = 536, rho = 0)
#'
#' ## Equivalence bounds for alpha = 0.05, N=536, statistical power of
#' ## 0.9, assuming true effect = 0
#' #powerTOSTr(alpha=0.05, N=536, statistical_power=0.9)
#'
#' @importFrom stats pnorm pt qnorm qt integrate uniroot qnorm pnorm dnorm
#' @importFrom graphics abline plot points segments title
#' @family Correlations
#' @family power
#' @export

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
    pow_corr_tost(n = n,
                  rho = rho,
                  power = power,
                  null = null,
                  alpha = alpha)
  } else{
    pow_corr(n = n, rho = rho, power = power, null = null,
             alpha = alpha, alternative = alternative)
  }
}

#' @rdname power_z_cor
#' @export

powerTOSTr <- function(alpha,
                       statistical_power,
                       N,
                       low_eqbound_r,
                       high_eqbound_r) {
  lifecycle::deprecate_soft("0.4.0", "powerTOSTr()", "tsum_TOST()")
  if(missing(N)) {
    za <- qnorm(1-alpha) #one-sided, so alpha not alpha/2
    zb <- qnorm(1-((1-statistical_power)/2)) #Note beta/2
    C1 = 0.5 * log((1+low_eqbound_r)/(1-low_eqbound_r))
    C2 = 0.5 * log((1+high_eqbound_r)/(1-high_eqbound_r))
    NT1 <- ((za+zb)/C1)^2 + 3
    NT2 <- ((za+zb)/C2)^2 + 3
    N<-max(NT1,NT2)
    message(cat("The required sample size to achieve",100*statistical_power,"% power with equivalence bounds of",low_eqbound_r,"and",high_eqbound_r,"is",ceiling(N),"observations"))
    return(N)
  }
  if(missing(statistical_power)) {
    se <- 1/sqrt(N-3)
    C1 = 0.5 * log((1+low_eqbound_r)/(1-low_eqbound_r))
    C2 = 0.5 * log((1+high_eqbound_r)/(1-high_eqbound_r))
    statistical_power1 <- 2 * (pnorm((abs(C1)/se) - qnorm(1-alpha)) + pnorm(-(abs(C1)/se) - qnorm(1-alpha))) - 1
    statistical_power2 <- 2 * (pnorm((abs(C2)/se) - qnorm(1-alpha)) + pnorm(-(abs(C2)/se) - qnorm(1-alpha))) - 1
    statistical_power<-min(statistical_power1, statistical_power2)
    if(statistical_power<0) {statistical_power<-0}
    message(cat("The statistical power is",round(100*statistical_power,2),"% for equivalence bounds of",low_eqbound_r,"and",high_eqbound_r,"."))
    return(statistical_power)
  }
  if(missing(low_eqbound_r) && missing(high_eqbound_r)) {
    za <- qnorm(1-alpha) #one-sided, so alpha not alpha/2
    zb <- qnorm(1-((1-statistical_power)/2)) #Note beta/2
    C = (za+zb)/sqrt(N - 3)
    low_eqbound_r <- -(exp(1)^(2*C)-1)/(exp(1)^(2*C) + 1)
    high_eqbound_r <- (exp(1)^(2*C)-1)/(exp(1)^(2*C) + 1)
    message(cat("The equivalence bounds to achieve",100*statistical_power,"% power with N =",N,"are",round(low_eqbound_r,2),"and",round(high_eqbound_r,2),"."))
    bounds<-c(low_eqbound_r,high_eqbound_r)
    return(bounds)
  }
}

pow_corr = function (n = NULL, rho = NULL, power = NULL, null = 0,
                     alpha = NULL, alternative = c("two.sided", "less", "greater"))
{

  if (sum(sapply(list(n, rho, power, alpha), is.null)) != 1)
    stop("exactly one of n, rho, power, and alpha must be NULL")
  if (!is.null(alpha) && !is.numeric(alpha) || any(0 > alpha |
                                                   alpha > 1))
    stop(sQuote("alpha"), " must be numeric in [0, 1]")
  if (!is.null(power) && !is.numeric(power) || any(0 > power |
                                                   power > 1))
    stop(sQuote("power"), " must be numeric in [0, 1]")
  if (!is.null(n) && n < 4)
    stop("number of observations must be at least 4")
  p=0
  alternative <- match.arg(alternative)
  tside <- switch(alternative, less = 1, two.sided = 2, greater = 3)
  if (tside == 2 && !is.null(rho))
    rho <- abs(rho)
  if (tside == 3) {
    p.body <- quote({
      delta <- sqrt(n - 3 - p) * (log((1 + rho)/(1 - rho))/2 +
                                    rho/(n - 1 - p)/2 * (1 + (5 + rho^2)/(n - 1 - p)/4 +
                                                           (11 + 2 * rho^2 + 3 * rho^4)/(n - 1 - p)^2/8) -
                                    log((1 + null)/(1 - null))/2 - null/(n - 1 -
                                                                           p)/2)
      v <- (n - 3 - p)/(n - 1 - p) * (1 + (4 - rho^2)/(n -
                                                         1 - p)/2 + (22 - 6 * rho^2 - 3 * rho^4)/(n - 1 -
                                                                                                    p)^2/6)
      zalpha <- qnorm(1 - alpha)
      pnorm((delta - zalpha)/sqrt(v))
    })
  }
  if (tside == 1) {
    p.body <- quote({
      delta <- sqrt(n - 3 - p) * (log((1 + rho)/(1 - rho))/2 +
                                    rho/(n - 1 - p)/2 * (1 + (5 + rho^2)/(n - 1 - p)/4 +
                                                           (11 + 2 * rho^2 + 3 * rho^4)/(n - 1 - p)^2/8) -
                                    log((1 + null)/(1 - null))/2 - null/(n - 1 -
                                                                           p)/2)
      v <- (n - 3 - p)/(n - 1 - p) * (1 + (4 - rho^2)/(n -
                                                         1 - p)/2 + (22 - 6 * rho^2 - 3 * rho^4)/(n - 1 -
                                                                                                    p)^2/6)
      zalpha <- qnorm(1 - alpha)
      pnorm((-delta - zalpha)/sqrt(v))
    })
  }
  if (tside == 2) {
    p.body <- quote({
      delta <- sqrt(n - 3 - p) * (log((1 + rho)/(1 - rho))/2 +
                                    rho/(n - 1 - p)/2 * (1 + (5 + rho^2)/(n - 1 - p)/4 +
                                                           (11 + 2 * rho^2 + 3 * rho^4)/(n - 1 - p)^2/8) -
                                    log((1 + null)/(1 - null))/2 - null/(n - 1 -
                                                                           p)/2)
      v <- (n - 3 - p)/(n - 1 - p) * (1 + (4 - rho^2)/(n -
                                                         1 - p)/2 + (22 - 6 * rho^2 - 3 * rho^4)/(n - 1 -
                                                                                                    p)^2/6)
      zalpha <- qnorm(1 - alpha/2)
      pnorm((delta - zalpha)/sqrt(v)) + pnorm((-delta -
                                                 zalpha)/sqrt(v))
    })
  }
  if (is.null(power))
    power <- eval(p.body)
  else if (is.null(n))
    n <- uniroot(function(n) eval(p.body) - power, c(4 +
                                                       p + 1e-10, 1e+07))$root
  else if (is.null(rho)) {
    if (tside == 2) {
      rho <- uniroot(function(rho) eval(p.body) - power, c(1e-10,
                                                           1 - 1e-10))$root
    }
    else {
      rho <- uniroot(function(rho) eval(p.body) - power, c(-1 +
                                                             1e-10, 1 - 1e-10))$root
    }
  }
  else if (is.null(alpha))
    alpha <- uniroot(function(alpha) eval(p.body) - power,
                     c(1e-10, 1 - 1e-10))$root
  else stop("internal error")
  METHOD <- "Approximate Power for Pearson Product-Moment Correlation (z-test)"

  structure(list(n = n, rho = rho,
                 alpha = alpha, beta = 1-power, power = power,
                 null = null, alternative = alternative,
                 method = METHOD),
            class = "power.htest")
}

pow_corr_tost = function (n = NULL, rho = 0, power = NULL, null = NULL,
                          alpha = NULL)
{

  if (sum(sapply(list(n, rho, power, alpha), is.null)) != 1)
    stop("exactly one of n, rho, power, and alpha must be NULL")
  if(is.null(rho)){
    stop("rho cannot be set to NULL at this time.")
  }
  if(is.null(null)){
    stop("null cannot be set to NULL at this time.")
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

  p.body <- quote({
    se <- 1/sqrt(n-3)
    C1 = 0.5 * log((1+min(null))/(1-min(null)))
    C2 = 0.5 * log((1+max(null))/(1-max(null)))
    statistical_power1 <- 2 * (pnorm((abs(C1)/se) - qnorm(1-alpha)) + pnorm(-(abs(C1)/se) - qnorm(1-alpha))) - 1
    statistical_power2 <- 2 * (pnorm((abs(C2)/se) - qnorm(1-alpha)) + pnorm(-(abs(C2)/se) - qnorm(1-alpha))) - 1
    min(statistical_power1, statistical_power2)
  })

  if (is.null(power))
    power <- eval(p.body)
  else if (is.null(n))
    n <- uniroot(function(n) eval(p.body) - power, c(4 + 1e-10, 1e+07))$root
  else if (is.null(alpha))
    alpha <- uniroot(function(alpha) eval(p.body) - power,
                     c(1e-10, 1 - 1e-10))$root
  else stop("internal error")
  METHOD <- "Approximate Power for Pearson Product-Moment Correlation (z-test)"

  structure(list(n = n, rho = rho,
                 alpha = alpha, beta = 1-power, power = power,
                 null = null, alternative = alternative,
                 method = METHOD),
            class = "power.htest")
}

