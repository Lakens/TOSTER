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
#'   * "auto" (default): `"stud"` for `method = "pearson"` and `"bca"` for all
#'     other methods. In simulations, the studentized interval kept Pearson's r
#'     near nominal coverage for skewed, heteroscedastic, and heavy-tailed data,
#'     where the BCa interval under-covered; for the rank and robust
#'     correlations BCa performed as well or better (see the "Studentized
#'     bootstrap" section). The method actually used is returned in `boot_ci`.
#'   * "bca": bias-corrected and accelerated bootstrap CI. Provides
#'     second-order accuracy by correcting for bias and skewness, but requires
#'     additional computation via the jackknife (n extra evaluations of the
#'     statistic).
#'   * "stud": studentized (bootstrap-t) CI. Each bootstrap replicate is
#'     standardized by a standard error estimated from that replicate's data
#'     (see the "Studentized bootstrap" section below). Only available for
#'     `method = "pearson"`, `"kendall"`, or `"spearman"`.
#'   * "basic": basic/empirical bootstrap CI
#'   * "perc": percentile bootstrap CI
#' @param R number of bootstrap replications (default = 1999).
#' @param boot_scale scale on which the `"basic"` and `"stud"` intervals and
#'   p-values are computed before being back-transformed to the correlation
#'   scale:
#'   * "z": Fisher z scale, \eqn{z = \mathrm{atanh}(r)} (default). Keeps
#'     the interval within \[-1, 1\] and is closer to variance-stabilizing.
#'   * "r": correlation scale. Studentized limits are truncated to \[-1, 1\];
#'     basic limits are not.
#'
#'   Ignored for `"perc"` and `"bca"`, which give the same result on either
#'   scale.
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
#'   (`2 * est - bvec`), which is the exact inversion of the basic CI. Both are
#'   computed on the scale chosen by `boot_scale`.
#'
#' * `boot_ci = "bca"`: p-values are derived from the BCa probability transformation,
#'   using the same bias correction and acceleration parameters as the BCa CI.
#'
#' * `boot_ci = "stud"`: p-values are derived from the bootstrap pivot
#'   distribution on the scale chosen by `boot_scale`. Each replicate's pivot is
#'   `(t_star - t_obs) / se_star`, where `se_star` is the standard error
#'   estimated from that replicate.
#'
#' @section Studentized bootstrap:
#'
#' A studentized (bootstrap-t) interval only improves on the basic interval if
#' the standard error is re-estimated from each bootstrap resample. The
#' normal-theory standard errors on the Fisher z scale (e.g.,
#' \eqn{1/\sqrt{n-3}} for Pearson's r) depend only on \eqn{n} (and, for
#' Spearman, on the estimate itself), so they cannot serve this purpose.
#' Instead, `boot_cor_test()` uses an influence-function (sandwich) standard
#' error that is estimated from the data and does not assume bivariate
#' normality:
#'
#' \deqn{\widehat{SE}(\hat\rho) = \frac{1}{n}\sqrt{\sum_{i=1}^{n} \psi_i^2},}
#'
#' where \eqn{\psi_i} is the empirical influence of observation \eqn{i} on the
#' coefficient. The same formula is applied to the observed data and to every
#' bootstrap resample. With `boot_scale = "z"`, the standard error is moved to
#' the Fisher z scale by the delta method,
#' \eqn{\widehat{SE}(\hat z) = \widehat{SE}(\hat\rho) / (1 - \hat\rho^2)}.
#'
#' **Pearson's r** uses the asymptotic distribution-free (fourth-moment)
#' influence function (Steiger & Hakstian, 1982). With
#' \eqn{z_{x,i} = (x_i - \bar x)/s_x} and \eqn{z_{y,i} = (y_i - \bar y)/s_y}
#' (standardized with denominator \eqn{n}),
#' \deqn{\psi^{0}_i = z_{x,i} z_{y,i} - \frac{r}{2}\left(z_{x,i}^2 + z_{y,i}^2\right).}
#' Under bivariate normality this gives the familiar
#' \eqn{\mathrm{Var}(r) \approx (1 - \rho^2)^2 / n}. On its own, this standard
#' error is too small in small samples and with heavy-tailed data, because a
#' sample under-represents the extreme points that dominate the variance of
#' \eqn{r}. Following the HC4 heteroscedasticity-consistent estimator for
#' regression (Cribari-Neto, 2004), which Wilcox (2017) also uses for testing
#' correlations, each influence value is therefore inflated according to its
#' leverage:
#' \deqn{\psi_i = \psi^{0}_i\,(1 - h_i)^{-\delta_i/2}, \qquad h_i = \frac{1}{n} + \frac{D_i^2}{n-1}, \qquad \delta_i = \min\left(4, \frac{n h_i}{3}\right),}
#' where \eqn{D_i} is the Mahalanobis distance of \eqn{(x_i, y_i)} from the
#' bivariate mean (so the \eqn{h_i} sum to 3), and \eqn{h_i} is capped at 0.99.
#' This leverage correction of the ADF influence function is an adaptation
#' made for this package, chosen because in simulations (normal,
#' heteroscedastic, \eqn{t_5}, \eqn{t_3}, and discrete data; n = 30 and 80) it
#' kept coverage of the studentized interval near nominal, including for
#' \eqn{t_3} data, where the uncorrected standard error under-covered badly. It
#' is slightly conservative for strongly heteroscedastic data. For data with
#' infinite fourth moments, no standard error of \eqn{r} is consistent, and a
#' robust correlation (`"winsorized"` or `"bendpercent"`) is still the better
#' choice.
#'
#' **Spearman's rho** is Pearson's r computed on the mid-distribution
#' transforms \eqn{u_i = (R(x_i) - 1/2)/n} and \eqn{v_i = (R(y_i) - 1/2)/n},
#' where \eqn{R} denotes midranks. Its influence function adds, to the Pearson
#' influence function for \eqn{(u, v)}, terms for estimating the transforms.
#' Writing \eqn{\omega_{ij}^x = \mathbf{1}(x_j > x_i) + \frac{1}{2}\mathbf{1}(x_j = x_i)}
#' (and likewise \eqn{\omega_{ij}^y}), the influence of observation \eqn{i} on
#' each moment is
#' \deqn{d_i(\overline{uv}) = u_i v_i + \frac{1}{n}\sum_j \omega_{ij}^x v_j + \frac{1}{n}\sum_j \omega_{ij}^y u_j - 3\,\overline{uv},}
#' \deqn{d_i(\overline{u}) = u_i + \frac{1}{n}\sum_j \omega_{ij}^x - 2\bar u, \qquad d_i(\overline{u^2}) = u_i^2 + \frac{2}{n}\sum_j \omega_{ij}^x u_j - 3\,\overline{u^2},}
#' (similarly for \eqn{v}), and \eqn{\psi_i} follows from the delta method for
#' \eqn{r_s = C / \sqrt{S_u S_v}}, with \eqn{C = \overline{uv} - \bar u \bar v}
#' and \eqn{S_u = \overline{u^2} - \bar u^2}:
#' \deqn{\psi_i = \frac{d_i(C)}{\sqrt{S_u S_v}} - \frac{r_s}{2}\left(\frac{d_i(S_u)}{S_u} + \frac{d_i(S_v)}{S_v}\right),}
#' where \eqn{d_i(C) = d_i(\overline{uv}) - \bar v\, d_i(\bar u) - \bar u\, d_i(\bar v)}
#' and \eqn{d_i(S_u) = d_i(\overline{u^2}) - 2 \bar u\, d_i(\bar u)}.
#' Without ties this is the influence function of Spearman's rho given by
#' Croux and Dehon (2010); the tie terms keep the standard error accurate for
#' discrete data (and for bootstrap resamples, which always contain ties).
#'
#' **Kendall's tau** (tau-b, as returned by `cor(method = "kendall")`) is a
#' ratio of U-statistics,
#' \eqn{\tau_b = \bar a / \sqrt{\bar b^x \bar b^y}}, with
#' \deqn{a_i = \frac{1}{n-1}\sum_{j \ne i} \mathrm{sign}(x_i - x_j)\,\mathrm{sign}(y_i - y_j), \qquad b_i^x = \frac{1}{n-1}\sum_{j \ne i} \mathbf{1}(x_i \ne x_j),}
#' and \eqn{b_i^y} defined likewise. Using Hoeffding's (1948) projection and
#' the delta method,
#' \deqn{\psi_i = 2\left[\frac{a_i - \bar a}{\sqrt{\bar b^x \bar b^y}} - \frac{\tau_b}{2}\left(\frac{b_i^x - \bar b^x}{\bar b^x} + \frac{b_i^y - \bar b^y}{\bar b^y}\right)\right].}
#' Without ties this reduces to the usual U-statistic variance
#' \eqn{\frac{4}{n}\widehat{\mathrm{Var}}(a_i)}. Its computation is
#' \eqn{O(n^2)} per resample, so the studentized Kendall interval is slower
#' for large samples.
#'
#' Under independence, these standard errors reduce to the familiar
#' \eqn{\mathrm{Var}(r) \approx \mathrm{Var}(r_s) \approx 1/n} and
#' \eqn{\mathrm{Var}(\tau) \approx 4/(9n)}.
#'
#' Because the pivots divide by an estimated standard error, studentized
#' intervals can be wide or erratic in small samples (roughly n < 30), where
#' the standard error of each resample is itself imprecise. This matters most
#' for the rank correlations. Kendall's standard error, for example, is nearly
#' unbiased even at n = 20, but it varies a lot between samples (especially
#' with ties) and is strongly correlated with the estimate. In simulations at
#' n = 30, the studentized Spearman and Kendall intervals under-covered by up
#' to about 2.5 percentage points for heavy-tailed or tied data, while
#' `boot_ci = "bca"` stayed closer to nominal; at n = 80 both were accurate.
#' For Spearman's rho and Kendall's tau with fewer than about 50 observations,
#' `boot_ci = "bca"` (the `"auto"` choice for these methods) is therefore
#' preferable.
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
#' * **stderr**: the bootstrap standard error of the correlation coefficient. For
#'   `boot_ci = "stud"`, a second element gives the influence-function standard
#'   error of the observed estimate on the working scale (`z.se` or `r.se`).
#' * **null.value**: the value(s) of the correlation under the null hypothesis.
#' * **alternative**: character string indicating the alternative hypothesis.
#' * **method**: a character string indicating which bootstrapped correlation was measured.
#' * **data.name**: a character string giving the names of the data.
#' * **boot_res**: vector of bootstrap correlation estimates.
#' * **boot_ci**: character string indicating which bootstrap CI method was used
#'   (with `boot_ci = "auto"`, the method it selected).
#' * **boot_scale**: the scale (`"z"` or `"r"`) on which the CI and p-value were
#'   computed. Always `"r"` for `"perc"` and `"bca"`.
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
#'
#' Wilcox, R. R. (1994). The Percentage Bend Correlation Coefficient. Psychometrika, 59(4), 601–616. https://doi.org/10.1007/bf02294395
#'
#' Wilcox, R. R. (1993). Some results on a Winsorized correlation coefficient. British Journal of Mathematical and Statistical Psychology, 46(2), 339–349. https://doi.org/10.1111/j.2044-8317.1993.tb01020.x
#'
#' Wilcox, R.R. (2017) Introduction to Robust Estimation and Hypothesis Testing, 4th edition. Academic Press.
#'
#' Cribari-Neto, F. (2004). Asymptotic inference under heteroskedasticity of unknown form. Computational Statistics & Data Analysis, 45(2), 215–233. https://doi.org/10.1016/S0167-9473(02)00366-3
#'
#' Croux, C., & Dehon, C. (2010). Influence functions of the Spearman and Kendall correlation measures. Statistical Methods & Applications, 19(4), 497–515. https://doi.org/10.1007/s10260-010-0142-z
#'
#' Efron, B., & Tibshirani, R. J. (1993). An Introduction to the Bootstrap. Chapman & Hall/CRC.
#'
#' Hoeffding, W. (1948). A class of statistics with asymptotically normal distribution. The Annals of Mathematical Statistics, 19(3), 293–325. https://doi.org/10.1214/aoms/1177730196
#'
#' Steiger, J. H., & Hakstian, A. R. (1982). The asymptotic distribution of elements of a correlation matrix: Theory and application. British Journal of Mathematical and Statistical Psychology, 35(2), 208–215. https://doi.org/10.1111/j.2044-8317.1982.tb00653.x
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
                          boot_ci = c("auto", "bca", "stud", "basic", "perc"),
                          R = 1999,
                          boot_scale = c("z", "r"),
                          ...) {
  boot_ci = match.arg(boot_ci)
  boot_scale = match.arg(boot_scale)
  DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
  alternative = match.arg(alternative)

  method = match.arg(method)

  # Default CI method by correlation type (see simulation notes in docs)
  if (boot_ci == "auto") {
    boot_ci <- if (method == "pearson") "stud" else "bca"
  }

  if (boot_ci == "stud" && method %in% c("winsorized", "bendpercent")) {
    stop(
      "Studentized bootstrap requires a closed-form SE and is only available ",
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

  # Working scale for scale-dependent intervals (basic, stud) -----
  # Percentile and BCa intervals are invariant to monotone transformations,
  # so they are always computed on the correlation scale.
  use_z <- boot_scale == "z" && boot_ci %in% c("basic", "stud")
  to_work <- if (use_z) atanh else identity
  from_work <- if (use_z) tanh else identity
  b_work <- to_work(bvec)
  est_work <- to_work(est)

  # Pivots for studentized bootstrap -----
  # The SE is re-estimated from every resample (influence-function SE)
  tvec <- NULL
  se_obs <- NULL
  if (boot_ci == "stud") {
    se_obs <- .cor_se(x, y, method)
    se_star <- apply(data, 1, function(i) .cor_se(x[i], y[i], method))
    if (use_z) {
      se_obs <- se_obs / (1 - est^2)
      se_star <- se_star / (1 - bvec^2)
    }
    if (!is.finite(se_obs) || se_obs <= 0) {
      stop(
        "Studentized bootstrap is undefined when the observed correlation ",
        "is +/-1 or its standard error is zero.",
        call. = FALSE
      )
    }
    tvec <- (b_work - est_work) / se_star
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
                     "basic" = from_work(basic(b_work, t0 = est_work, alpha2)),
                     "perc" = perc(bvec, alpha2),
                     "bca" = bca_ci(boots_est = bvec, t0 = est,
                                    jack_est = jack_est, alpha = alpha2),
                     "stud" = stud_ci(tvec, t0 = est_work, se_obs = se_obs,
                                      alpha = alpha2, back = from_work))
  attr(boot.cint, "conf.level") <- ci

  # P-value computation (method-consistent, on the same working scale as the CI)
  pval <- function(null, alternative) {
    boot_pvalue(bvec = b_work, est = est_work, null = to_work(null),
                alternative = alternative, boot_ci = boot_ci,
                tvec = tvec, se_obs = se_obs,
                z0 = z0, acc = acc, nboot = nboot)
  }
  if (alternative %in% c("two.sided", "greater", "less")) {
    sig <- pval(null.value, alternative)
  } else if (alternative == "equivalence") {
    sig <- max(pval(min(null.value), "greater"),
               pval(max(null.value), "less"))
  } else if (alternative == "minimal.effect") {
    sig <- min(pval(max(null.value), "greater"),
               pval(min(null.value), "less"))
  }

  # CI method label for method string
  ci_label <- switch(boot_ci,
                     "basic" = " (basic)",
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

  # SE: for stud, report both the bootstrap SE and the observed
  # influence-function SE on the working scale
  if (boot_ci == "stud") {
    se_out <- c(boot.se = sd(bvec, na.rm = TRUE), se_obs)
    names(se_out)[2] <- if (use_z) "z.se" else "r.se"
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
               boot_scale = if (use_z) "z" else "r",
               call = match.call())
  class(rval) <- "htest"
  return(rval)
}



