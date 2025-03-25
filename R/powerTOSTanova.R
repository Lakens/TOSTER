#' Power Analysis for F-test Equivalence Testing
#' @description
#' `r lifecycle::badge('stable')`
#'
#' Performs power analysis for equivalence testing with F-tests (ANOVA models). This function
#' calculates statistical power, sample size, equivalence bound, or alpha level when the
#' other parameters are specified.
#'
#' @param alpha Significance level (Type I error rate). Default is 0.05.
#' @param df1 Numerator degrees of freedom (e.g., groups - 1 for one-way ANOVA).
#' @param df2 Denominator degrees of freedom (e.g., N - groups for one-way ANOVA),
#'   where N is the total sample size.
#' @param eqbound Equivalence bound for partial eta-squared. This represents the threshold
#'   for what effect size would be considered practically insignificant.
#' @param power Desired statistical power (1 - Type II error rate). Default is NULL.
#'
#' @details
#' This function provides power analysis for the omnibus non-inferiority testing procedure
#' described by Campbell & Lakens (2021). Exactly one of the parameters `alpha`, `df1`,
#' `df2`, `eqbound`, or `power` must be NULL, and the function will solve for that
#' parameter.
#'
#' For one-way ANOVA:
#' * `df1` = number of groups - 1
#' * `df2` = total N - number of groups
#'
#' Common equivalence bounds (we do not recommend their use for choosing equivalence bounds) for partial eta-squared based on Cohen's benchmarks:
#' * Small effect: 0.01
#' * Medium effect: 0.06
#' * Large effect: 0.14
#'
#' Note that this function is primarily validated for one-way ANOVA designs; use with
#' caution for more complex designs.
#'
#' @return
#' An object of class "power.htest" containing the following components:
#'
#' * **df1**: Numerator degrees of freedom
#' * **df2**: Denominator degrees of freedom
#' * **eqbound**: Equivalence bound for partial eta-squared
#' * **sig.level**: Significance level (alpha)
#' * **power**: Statistical power
#' * **method**: Description of the test
#'
#' @examples
#' # Example 1: Calculate power given degrees of freedom and equivalence bound
#' # For a one-way ANOVA with 3 groups, 80 subjects per group, and equivalence bound of 0.01
#' power_eq_f(df1 = 2, df2 = 237, eqbound = 0.01)
#'
#' # Example 2: Calculate required denominator df (related to sample size)
#' # for 80% power with equivalence bound of 0.05
#' power_eq_f(df1 = 2, power = 0.8, eqbound = 0.05)
#'
#' # Example 3: Calculate detectable equivalence bound with 80% power
#' power_eq_f(df1 = 2, df2 = 100, power = 0.8)
#'
#' # Example 4: Calculate required alpha level for 90% power
#' power_eq_f(df1 = 2, df2 = 100, eqbound = 0.05, power = 0.9)
#'
#' @references
#' Campbell, H., & Lakens, D. (2021). Can we disregard the whole model? Omnibus non‐inferiority
#' testing for R2 in multi‐variable linear regression and in ANOVA. British Journal of
#' Mathematical and Statistical Psychology, 74(1), 64-89. doi: 10.1111/bmsp.12201
#'
#' @importFrom stats pf qf uniroot optimize
#' @family power
#' @export

power_eq_f <- function(alpha = 0.05,
                       df1 = NULL,
                       df2 = NULL,
                       eqbound = NULL,
                       power = NULL){

  # Check that exactly one parameter is NULL
  if (sum(c(is.null(alpha), is.null(df1), is.null(df2),
            is.null(eqbound), is.null(power))) != 1) {
    stop("Exactly one of alpha, df1, df2, eqbound, or power must be NULL")
  }

  # Parameter validation
  if (!is.null(alpha) && (alpha <= 0 || alpha >= 1)) {
    stop("alpha must be between 0 and 1")
  }
  if (!is.null(power) && (power <= 0 || power >= 1)) {
    stop("power must be between 0 and 1")
  }
  if (!is.null(df1) && df1 <= 0) {
    stop("df1 must be positive")
  }
  if (!is.null(df2) && df2 <= 0) {
    stop("df2 must be positive")
  }
  if (!is.null(eqbound) && (eqbound <= 0 || eqbound >= 1)) {
    stop("eqbound must be between 0 and 1")
  }

  # Function to calculate power
  calc_power <- function(alpha, df1, df2, eqbound) {
    f2 = eqbound/(1 - eqbound)
    lambda = (f2 * (df1 + df2 + 1))

    Fstatstar = qf(
      alpha,
      df1 = df1,
      df2 = df2,
      ncp = lambda,
      lower.tail = TRUE
    )

    statistical_power = pf(Fstatstar,
                           df1 = df1,
                           df2 = df2,
                           lower.tail = TRUE)
    return(statistical_power)
  }

  # Solve for the missing parameter
  if (is.null(power)) {
    power <- calc_power(alpha, df1, df2, eqbound)
    message(paste("Power =", round(power, 4)))
  } else if (is.null(df2)) {
    # Find df2 that gives desired power
    f1 <- function(df2) calc_power(alpha, df1, df2, eqbound) - power
    df2 <- uniroot(f1, c(4, 1e5))$root
    message(paste("Required df2 =", round(df2, 2),
                  "(approximately", ceiling(df2 + df1 + 1),
                  "total observations for a one-way ANOVA with", df1 + 1, "groups)"))
  } else if (is.null(eqbound)) {
    # Find eqbound that gives desired power
    f2 <- function(eqbound) calc_power(alpha, df1, df2, eqbound) - power
    eqbound <- uniroot(f2, c(0.001, 0.999))$root
    message(paste("Detectable equivalence bound =", round(eqbound, 4)))
  } else if (is.null(alpha)) {
    # Find alpha that gives desired power
    f3 <- function(alpha) calc_power(alpha, df1, df2, eqbound) - power
    alpha <- uniroot(f3, c(1e-10, 1-1e-10))$root
    message(paste("Required alpha =", round(alpha, 4)))
  } else if (is.null(df1)) {
    # Find df1 that gives desired power
    f4 <- function(df1) calc_power(alpha, df1, df2, eqbound) - power
    df1 <- uniroot(f4, c(1, 100))$root
    message(paste("Required df1 =", round(df1, 2)))
  }

  # Warning message about validation
  message("Note: This function is primarily validated for one-way ANOVA; use with caution for more complex designs")

  # Return results
  structure(
    list(
      df1 = df1,
      df2 = df2,
      eqbound = eqbound,
      sig.level = alpha,
      power = power,
      method = "Power Analysis for F-test Equivalence Testing"
    ),
    class = "power.htest"
  )
}
