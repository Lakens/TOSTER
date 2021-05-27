#' @title Power calculations for TOST for t-tests
#' @description  Calculates the exact power of two one sided t-tests (TOST) for one, two, and paired samples.
#' @param n number of observations per group. 2 sample sizes, in a vector, can be provided for the two sample case.
#' @param delta true difference in means (default is 0)
#' @param sd population standard deviation. Standard deviation of the differences for paired samples
#' @param low_eqbound The lower equivalance bound (raw units)
#' @param high_eqbound The upper equivalence bound (raw units)
#' @param alpha a priori alpha-level (i.e., significance level)
#' @param power power of the TOST procedure (1-beta)
#' @param type string specifying the type of t-test.
#' @details
#' The exact calculations of power are based on Owen’s Q-function or by direct integration of the bivariate non-central t-distribution (inspired by the PowerTOST package).
#' Approximate power is implemented via the non-central t-distribution or the ‘shifted’ central t-distribution.
#' @note
#' The power function in this package is limited. Please see the PowerTOST R package for more options.
#' @references
#' Phillips KF. Power of the Two One-Sided Tests Procedure in Bioequivalence. J Pharmacokin Biopharm. 1990;18(2):137–44. doi: 10.1007/BF01063556
#'
#' Diletti D, Hauschke D, Steinijans VW. Sample Size Determination for Bioequivalence Assessment by Means of Confidence Intervals. Int J Clin Pharmacol Ther Toxicol. 1991;29(1):1–8.
#' @importFrom stats integrate uniroot qnorm pnorm dnorm
#' @export
#'
power_t_TOST <- function(
  n = NULL,
  delta = 0,
  sd = 1,
  low_eqbound = NULL,
  high_eqbound = NULL,
  alpha = NULL,
  power = NULL,
  type = "two.sample"
){
  if(is.null(low_eqbound) || is.null(high_eqbound)){
    stop("Equivalence bounds must be provided (low_eqbound and high_eqbound)")
  }

  if(low_eqbound > delta || high_eqbound < delta){
    stop("True mean difference greater than bounds. TOST power calculation not possible.")
  }

  if (sum(sapply(list(n, power, alpha), is.null)) !=
      1){
    stop("exactly one of n, d, power, and sig.level must be NULL")
  }

  # Create quote to evalulate with uniroot
  p.body <- quote({
    pow_tTOST(alpha = alpha,
              theta1 = low_eqbound,
              theta2 = high_eqbound,
              theta0 = delta,
              sd = sd,
              n = n,
              type = type)
  })

  if (is.null(power)){
    power <- eval(p.body)
  }else if (is.null(n)){
    n <- uniroot(function(n) eval(p.body) - power, c(2 +
                                                       1e-10, 1e+09))$root
  }else if (is.null(alpha)){
    alpha <- uniroot(function(alpha) eval(p.body) - power,
                    c(1e-10, 1 - 1e-10))$root
  }else {stop("internal error")}

  beta = 1-power

  NOTE <- switch(type,
                 paired = "n is number of *pairs*",
                 two.sample = "n is number in *each* group",
                 NULL)
  METHOD <- paste(switch(type,
                         one.sample = "One-sample",
                         two.sample = "Two-sample",
                         paired = "Paired"),
                  "TOST power calculation")
  structure(list(power = power,
                 beta = beta,
                 alpha = alpha,
                 n = n,
                 delta = delta,
                 sd = sd,
                 bounds = c(low_eqbound, high_eqbound),
                 note = NOTE, method = METHOD),
            class = "power.htest")
}
