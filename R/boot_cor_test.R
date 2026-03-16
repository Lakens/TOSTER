#' @title Bootstrapped Correlation Coefficients
#' @description
#' `r lifecycle::badge('stable')`
#'
#' A function for bootstrap-based correlation tests using various correlation coefficients
#' including Pearson's, Kendall's, Spearman's, Winsorized, and percentage bend correlations.
#' This function supports standard, equivalence, and minimal effect testing with robust bootstrap methods.
#'
#' @inheritParams boot_t_TOST
#' @inheritParams z_cor_test
#' @param method a character string indicating which correlation coefficient to use:
#'   * "pearson": standard Pearson product-moment correlation
#'   * "kendall": Kendall's tau rank correlation
#'   * "spearman": Spearman's rho rank correlation
#'   * "winsorized": Winsorized correlation (robust to outliers)
#'   * "bendpercent": percentage bend correlation (robust to marginal outliers)
#'
#'   Can be abbreviated.
#' @param boot_ci type of bootstrap confidence interval:
#'   * "basic": basic/empirical bootstrap CI (default)
#'   * "perc": percentile bootstrap CI
#'   * "bca": bias-corrected and accelerated bootstrap CI. Provides second-order
#'     accuracy by correcting for bias and skewness, but requires additional
#'     computation via the jackknife (n extra evaluations of the statistic).
#'   * "stud": studentized (bootstrap-t) CI. Uses pivot statistics on the Fisher z
#'     scale with analytical SEs. Only available for `method = "pearson"`,
#'     `"kendall"`, or `"spearman"`.
#' @param R number of bootstrap replications (default = 1999).
#' @param ... additional arguments passed to correlation functions, such as:
#'   * tr: trim for Winsorized correlation (default = 0.2)
#'   * beta: for percentage bend correlation (default = 0.2)
#'
#' @details
#' This function uses bootstrap methods to calculate correlation coefficients and their
#' confidence intervals. P-values are computed by inverting the selected CI method,
#' which guarantees that `p < alpha` if and only if the corresponding CI excludes the
#' null value.
#'
#' **P-value computation by CI method:**
#'
#' * `boot_ci = "perc"`: p-values are computed from the raw bootstrap distribution
#'   (proportion of replicates beyond the null). This is the original approach from
#'   Wilcox (2017).
#'
#' * `boot_ci = "basic"`: p-values use the reflected bootstrap distribution
#'   (`2 * est - bvec`), which is the exact inversion of the basic CI.
#'
#' * `boot_ci = "bca"`: p-values are derived from the BCa probability transformation,
#'   using the same bias correction and acceleration parameters as the BCa CI.
#'
#' * `boot_ci = "stud"`: p-values are derived from the bootstrap pivot distribution
#'   on the Fisher z scale. Each replicate's pivot is `(z_star - z_obs) / se_star`,
#'   where `se_star` is the analytical SE. This method is only available for
#'   Pearson, Kendall, and Spearman correlations, since robust methods lack
#'   analytical SEs on the Fisher z scale.
#'
#' The bootstrap correlation methods in this package offer two robust correlations beyond
#' the standard methods:
#'
#' 1. **Winsorized correlation**: Replaces extreme values with less extreme values before
#'    calculating the correlation. The `trim` parameter (default: `tr = 0.2`) determines the
#'    proportion of data to be Winsorized.
#'
#' 2. **Percentage bend correlation**: A robust correlation that downweights the influence
#'    of outliers. The `beta` parameter (default = 0.2) determines the bending constant.
#'
#' These calculations are based on Rand Wilcox's R functions for his book (Wilcox, 2017),
#' and adapted from their implementation in Guillaume Rousselet's R package "bootcorci".
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
#' See `vignette("correlations")` for more details.
#'
#' @return A list with class "htest" containing the following components:
#'
#' * **p.value**: the bootstrap p-value of the test.
#' * **parameter**: the number of observations used in the test.
#' * **conf.int**: a bootstrap confidence interval for the correlation coefficient.
#' * **estimate**: the estimated correlation coefficient, with name "r", "tau", "rho", "pb", or "wincor"
#'   corresponding to the method employed.
#' * **stderr**: the bootstrap standard error of the correlation coefficient.
#' * **null.value**: the value(s) of the correlation under the null hypothesis.
#' * **alternative**: character string indicating the alternative hypothesis.
#' * **method**: a character string indicating which bootstrapped correlation was measured.
#' * **data.name**: a character string giving the names of the data.
#' * **boot_res**: vector of bootstrap correlation estimates.
#' * **boot_ci**: character string indicating which bootstrap CI method was used.
#' * **call**: the matched call.
#'
#' @examples
#' # Example 1: Standard bootstrap test with Pearson correlation
#' x <- c(44.4, 45.9, 41.9, 53.3, 44.7, 44.1, 50.7, 45.2, 60.1)
#' y <- c( 2.6,  3.1,  2.5,  5.0,  3.6,  4.0,  5.2,  2.8,  3.8)
#' boot_cor_test(x, y, method = "pearson", alternative = "two.sided",
#'               R = 999) # Fewer replicates for example
#'
#' # Example 2: Equivalence test with Spearman correlation
#' # Testing if correlation is equivalent to zero within ±0.3
#' boot_cor_test(x, y, method = "spearman", alternative = "equivalence",
#'              null = 0.3, R = 999)
#'
#' # Example 3: Using robust correlation methods
#' # Using Winsorized correlation with custom trim
#' boot_cor_test(x, y, method = "winsorized", tr = 0.1,
#'              R = 999)
#'
#' # Example 4: Using percentage bend correlation
#' boot_cor_test(x, y, method = "bendpercent", beta = 0.2,
#'              R = 999)
#'
#' # Example 5: Minimal effect test with asymmetric bounds
#' # Testing if correlation is outside bounds of -0.1 and 0.4
#' boot_cor_test(x, y, method = "pearson", alternative = "minimal.effect",
#'              null = c(-0.1, 0.4), R = 999)
#'
#' @section References:
#' Wilcox, R.R. (2009) Comparing Pearson Correlations: Dealing with Heteroscedasticity and Nonnormality.
#' Communications in Statistics - Simulation and Computation, 38, 2220–2234.
#'
#' Wilcox, R.R. (2017) Introduction to Robust Estimation and Hypothesis Testing, 4th edition. Academic Press.
#'
#' @family Correlations
#' @export

boot_cor_test <- function(x,
                          y,
                          alternative = c("two.sided", "less", "greater",
                                          "equivalence", "minimal.effect"),
                          method = c("pearson", "kendall", "spearman",
                                     "winsorized", "bendpercent"),
                          alpha = 0.05,
                          null = 0,
                          boot_ci = c("bca","stud", "basic", "perc"),
                          R = 1999,
                          ...) {
  boot_ci = match.arg(boot_ci)
  DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
  alternative = match.arg(alternative)

  method = match.arg(method)

  if (boot_ci == "stud" && method %in% c("winsorized", "bendpercent")) {
    stop(
      "Studentized bootstrap requires an analytical SE and is only available ",
      "for method = 'pearson', 'kendall', or 'spearman'.",
      call. = FALSE
    )
  }
  nboot = R
  null.value = null
  if(!is.vector(x) || !is.vector(y)){
    stop("x and y must be vectors.")
  }
  if(length(x)!=length(y)){
    stop("the vectors do not have equal lengths.")
  }
  df <- cbind(x,y)
  df <- df[complete.cases(df), ]
  n <- nrow(df)
  x <- df[,1]
  y <- df[,2]

  if(alternative %in% c("equivalence", "minimal.effect")){
    if(length(null) == 1){
      null.value = c(null.value, -1*null.value)
    }
    TOST = TRUE
  } else {
    if(length(null) > 1){
      stop("null can only have 1 value for non-TOST procedures")
    }
    TOST = FALSE
  }

  #if(TOST && null <=0){
  #  stop("positive value for null must be supplied if using TOST.")
  #}
  #if(TOST){
  #  alternative = "less"
  #}

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

  if(method %in% c("bendpercent","winsorized")){
    if(method == "bendpercent"){
      est <- pbcor(x, y, ...)
      data <- matrix(sample(n, size=n*nboot, replace=TRUE), nrow=nboot)
      bvec <- apply(data, 1, .corboot_pbcor, x, y, ...) # get bootstrap results corr
    }

    if(method == "winsorized"){
      est <- wincor(x, y, ...)
      data <- matrix(sample(n, size=n*nboot, replace=TRUE), nrow=nboot)
      bvec <- apply(data, 1, .corboot_wincor, x, y, ...) # get bootstrap results corr
    }

  } else {
    est <- cor(x, y, method = method)
    data <- matrix(sample(n, size=n*nboot, replace=TRUE), nrow=nboot)
    bvec <- apply(data, 1, .corboot, x, y, method = method, ...) # get bootstrap results corr
  }

  alpha2 = ifelse(alternative != "two.sided",
                  alpha*2,
                  alpha)

  # Compute pivots for studentized bootstrap
  tvec <- NULL
  se_obs <- NULL
  if (boot_ci == "stud") {
    se_obs <- .fisher_z_se(est, n, method)
    se_star <- .fisher_z_se(bvec, n, method)
    z_star <- atanh(bvec)
    z_obs <- atanh(est)
    tvec <- (z_star - z_obs) / se_star
  }

  # Jackknife for BCa (if needed)
  z0 <- NULL
  acc <- NULL
  if (boot_ci == "bca") {
    jack_est <- numeric(n)
    for (j in seq_len(n)) {
      if (method == "bendpercent") {
        jack_est[j] <- pbcor(x[-j], y[-j], ...)
      } else if (method == "winsorized") {
        jack_est[j] <- wincor(x[-j], y[-j], ...)
      } else {
        jack_est[j] <- cor(x[-j], y[-j], method = method)
      }
    }
    # Pre-compute BCa parameters for p-value use
    bca_par <- bca_params(bvec, est, jack_est)
    z0 <- bca_par$z0
    acc <- bca_par$acc
  }

  # CI computation
  boot.cint = switch(boot_ci,
                     "basic" = basic(bvec, t0 = est, alpha2),
                     "perc" = perc(bvec, alpha2),
                     "bca" = bca_ci(boots_est = bvec, t0 = est,
                                    jack_est = jack_est, alpha = alpha2),
                     "stud" = stud_ci(tvec, t0_z = atanh(est),
                                      se_obs = se_obs, alpha = alpha2))
  attr(boot.cint, "conf.level") <- ci

  # P-value computation (method-consistent)
  if (alternative %in% c("two.sided", "greater", "less")) {
    sig <- boot_pvalue(bvec = bvec, est = est, null = null.value,
                       alternative = alternative, boot_ci = boot_ci,
                       tvec = tvec, se_obs = se_obs,
                       z0 = z0, acc = acc, nboot = nboot,
                       z_transform= TRUE)
  } else if (alternative == "equivalence") {
    sig1 <- boot_pvalue(bvec = bvec, est = est, null = min(null.value),
                        alternative = "greater", boot_ci = boot_ci,
                        tvec = tvec, se_obs = se_obs,
                        z0 = z0, acc = acc, nboot = nboot,
                        z_transform= TRUE)
    sig2 <- boot_pvalue(bvec = bvec, est = est, null = max(null.value),
                        alternative = "less", boot_ci = boot_ci,
                        tvec = tvec, se_obs = se_obs,
                        z0 = z0, acc = acc, nboot = nboot,
                        z_transform= TRUE)
    sig <- max(sig1, sig2)
  } else if (alternative == "minimal.effect") {
    sig1 <- boot_pvalue(bvec = bvec, est = est, null = max(null.value),
                        alternative = "greater", boot_ci = boot_ci,
                        tvec = tvec, se_obs = se_obs,
                        z0 = z0, acc = acc, nboot = nboot,
                        z_transform= TRUE)
    sig2 <- boot_pvalue(bvec = bvec, est = est, null = min(null.value),
                        alternative = "less", boot_ci = boot_ci,
                        tvec = tvec, se_obs = se_obs,
                        z0 = z0, acc = acc, nboot = nboot,
                        z_transform= TRUE)
    sig <- min(sig1, sig2)
  }

  # CI method label for method string
  ci_label <- switch(boot_ci,
                     "basic" = "(basic)",
                     "perc" = " (percentile)",
                     "bca" = " (BCa)",
                     "stud" = " (studentized)")

  if (method == "pearson") {
    method2 <- paste0("Bootstrapped Pearson's product-moment correlation", ci_label)
    names(null.value) = rep("correlation",length(null.value))
    rfinal = c(r = est)
  }
  if (method == "spearman") {
    method2 <- paste0("Bootstrapped Spearman's rank correlation rho", ci_label)
    rfinal = c(rho = est)
    names(null.value) = rep("rho",length(null.value))
  }
  if (method == "kendall") {
    method2 <- paste0("Bootstrapped Kendall's rank correlation tau", ci_label)
    rfinal = c(tau = est)
    names(null.value) = rep("tau",length(null.value))
  }
  if (method == "bendpercent") {
    method2 <- paste0("Bootstrapped percentage bend correlation pb", ci_label)
    rfinal = c(pb = est)
    names(null.value) = rep("pb",length(null.value))
  }
  if (method == "winsorized") {
    method2 <- paste0("Bootstrapped Winsorized correlation wincor", ci_label)
    rfinal = c(wincor = est)
    names(null.value) = rep("wincor",length(null.value))
  }
  N = n
  names(N) = "N"

  # SE: for stud, report both bootstrap SE and analytical z-scale SE
  if (boot_ci == "stud") {
    se_out <- c(boot.se = sd(bvec, na.rm = TRUE), z.se = se_obs)
  } else {
    se_out <- sd(bvec, na.rm = TRUE)
  }

  # Store as htest
  rval <- list(p.value = sig,
               parameter = N,
               conf.int = boot.cint,
               estimate = rfinal,
               stderr = se_out,
               null.value = null.value,
               alternative = alternative,
               method = method2,
               data.name = DNAME,
               boot_res = bvec,
               boot_ci = boot_ci,
               call = match.call())
  class(rval) <- "htest"
  return(rval)
}



