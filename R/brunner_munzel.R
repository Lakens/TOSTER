#' @title Brunner-Munzel Test
#' @description
#' `r lifecycle::badge("maturing")`
#'
#' This is a generic function that performs a generalized asymptotic Brunner-Munzel test in a fashion similar to [t.test].
#' @param paired a logical indicating whether you want a paired test.
#' @param mu a number or vector specifying the null hypothesis value(s):
#'
#'   * For standard alternatives ("two.sided", "less", "greater"): a single value
#'     representing the hypothesized relative effect (default = 0.5, i.e., stochastic equality).
#'   * For "equivalence" or "minimal.effect": two values representing the lower and upper bounds
#'     for the relative effect. Values must be between 0 and 1.
#'
#' @param test_method a character string specifying the test method to use:
#'
#'   * "t" (default): approximate t-distribution with Satterthwaite-Welch degrees of freedom
#'   * "logit": logit transformation for range-preserving confidence intervals
#'   * "perm": studentized permutation test (recommended when sample size per condition is less than 15)
#'
#' @param R the number of permutations for the permutation test (default is 10000). Only used when `test_method = "perm"`.
#' @param p_method a character string specifying the method for computing permutation p-values. Only used when `test_method = "perm"`:
#'
#'   * "plusone" (default): Uses the (b+1)/(R+1) formula from Phipson & Smyth (2010), which provides
#'     exact p-values when permutations are sampled without replacement.
#'   * "original": Uses the traditional b/R formula with a minimum of 1/R.
#' @param perm `r lifecycle::badge("deprecated")` Use `test_method = "perm"` instead.
#' @param max_n_perm `r lifecycle::badge("deprecated")` Use `R` instead.
#' @param alternative a character string specifying the alternative hypothesis, must be one of:
#'
#'   * "two.sided" (default): true relative effect is not equal to `mu`
#'   * "less": true relative effect is less than `mu`
#'   * "greater": true relative effect is greater than `mu`
#'   * "equivalence": true relative effect is between the lower and upper bounds specified in `mu`
#'   * "minimal.effect": true relative effect is less than the lower bound or greater than the upper bound specified in `mu`
#'
#' @inheritParams t_TOST
#' @param ...  further arguments to be passed to or from methods.
#' @details
#'
#' This function is made to provide a test of stochastic equality between two samples (paired or independent), and is referred to as the Brunner-Munzel test.
#'
#' This tests the hypothesis that the relative effect, discussed below, is equal to the null value (default is `mu = 0.5`).
#'
#' The estimate of the relative effect, which can be considered as value similar to the probability of superiority, refers to the following:
#'
#'  \deqn{\hat p = P(X>Y) + \frac{1}{2} \cdot P(X=Y)}
#'
#'  Note, for paired samples, this does *not* refer to the probability of an increase/decrease in paired sample but rather the probability that a randomly sampled value of X is greater than a randomly sampled value of Y.
#'  This is also referred to as the "relative" effect in the literature. Therefore, the results will differ from the concordance probability provided by the ses_calc function.
#'
#'  The brunner_munzel function is based on the `npar.t.test` and `npar.t.test.paired` functions within the `nparcomp` package (Konietschke et al. 2015).
#'
#' ## Test Methods
#'
#' Three test methods are available:
#'
#' * "t": The default method uses a t-distribution approximation with Satterthwaite-Welch
#'   degrees of freedom. This is appropriate for moderate to large sample sizes.
#'
#' * "logit": Uses a logit transformation to produce range-preserving confidence intervals
#'   that are guaranteed to stay within `[0, 1]`. This method is recommended when the estimated
#'   relative effect is close to 0 or 1.
#'
#' * "perm": A studentized permutation test following Neubert & Brunner (2007). This method
#'   is highly recommended when sample sizes are small (< 15 per group) as it provides better
#'   control of Type I error rates in these situations.
#'
#' ## Hypothesis Testing
#'
#' For the standard alternatives, the null hypothesis is that the relative effect equals `mu`:
#'
#' * "two.sided": H0: p = mu vs H1: p ≠ mu
#' * "less": H0: p ≥ mu vs H1: p < mu
#' * "greater": H0: p ≤ mu vs H1: p > mu
#'
#' For equivalence and minimal effect testing using the two one-sided tests (TOST) procedure:
#'
#' * "equivalence": H0: p ≤ mu\[1\] OR p ≥ mu\[2\] vs H1: mu\[1\] < p < mu\[2\]
#'
#'   Tests whether the relative effect falls within the specified bounds.
#'   The p-value is the maximum of the two one-sided p-values.
#'
#' * "minimal.effect": H0: mu\[1\] < p < mu\[2\] vs H1: p ≤ mu\[1\] OR p ≥ mu\[2\]
#'
#'   Tests whether the relative effect falls outside the specified bounds.
#'   The p-value is the minimum of the two one-sided p-values.
#'
#' ## Test Statistic and P-value Calculation
#'
#' The test statistic is calculated as:
#'
#' \deqn{t = \sqrt{N} \cdot \frac{\hat{p} - p_0}{s}}
#'
#' where \eqn{N} is the total sample size (or \eqn{n} for paired samples), \eqn{\hat{p}} is the
#' estimated relative effect, \eqn{p_0} is the null hypothesis value, and \eqn{s} is the
#' rank-based standard error.
#'
#' For equivalence testing, two test statistics are computed:
#'
#' \deqn{t_{low} = \sqrt{N} \cdot \frac{\hat{p} - p_{low}}{s}}
#' \deqn{t_{high} = \sqrt{N} \cdot \frac{\hat{p} - p_{high}}{s}}
#'
#' where \eqn{p_{low}} and \eqn{p_{high}} are the lower and upper equivalence bounds.
#' The one-sided p-values are:
#'
#' * \eqn{p_1}: p-value for H1: p > \eqn{p_{low}} (from the lower bound test)
#' * \eqn{p_2}: p-value for H1: p < \eqn{p_{high}} (from the upper bound test)
#'
#' For equivalence: \eqn{p_{TOST} = \max(p_1, p_2)}
#'
#' For minimal effect: \eqn{p_{MET} = \min(1 - p_1, 1 - p_2)}
#'
#' ## Confidence Intervals for Equivalence Testing
#'
#' When `alternative = "equivalence"` or `alternative = "minimal.effect"`, the confidence
#' interval is computed at the \eqn{1 - 2\alpha} level (default: 90% CI when \eqn{\alpha = 0.05}).
#' This follows the standard TOST procedure where the \eqn{(1 - 2\alpha) \times 100\%} CI
#' corresponds to two one-sided tests at level \eqn{\alpha}.
#'
#' ## Permutation Tests with Non-0.5 Null Values
#'
#' When `test_method = "perm"` and `mu != 0.5`, the permutation distribution is constructed by centering
#' the permuted test statistics at 0.5 (the value implied by exchangeability), while the observed
#' test statistic is centered at the hypothesized null value. This approach is valid because
#' the studentized permutation distribution converges to the same limit regardless of the
#' centering, following the asymptotic theory of Janssen (1997) and Neubert & Brunner (2007).
#'
#' @return A list with class `"htest"` containing the following components:
#'
#'   - "statistic": the value of the test statistic.
#'   - "parameter": the degrees of freedom for the test statistic.
#'   - "p.value": the p-value for the test.
#'   - "conf.int": a confidence interval for the relative effect appropriate to the specified alternative hypothesis.
#'   - "estimate": the estimated relative effect.
#'   - "null.value": the specified hypothesized value of the relative effect.
#'   - "stderr": the standard error of the relative effect.
#'   - "alternative": a character string describing the alternative hypothesis.
#'   - "method": a character string indicating what type of test was performed.
#'   - "data.name": a character string giving the name(s) of the data.
#'
#' @examples
#' data(mtcars)
#' # Standard test of stochastic equality
#' brunner_munzel(mpg ~ am, data = mtcars)
#'
#' # Test using logit transformation for range-preserving CIs
#' brunner_munzel(mpg ~ am, data = mtcars, test_method = "logit")
#'
#' # Test against a specific null value
#' brunner_munzel(mpg ~ am, data = mtcars, mu = 0.3)
#'
#' # Equivalence test: is the relative effect between 0.35 and 0.65?
#' brunner_munzel(mpg ~ am, data = mtcars,
#'                alternative = "equivalence",
#'                mu = c(0.35, 0.65))
#'
#' # Minimal effect test: is the relative effect outside 0.4 to 0.6?
#' brunner_munzel(mpg ~ am, data = mtcars,
#'                alternative = "minimal.effect",
#'                mu = c(0.4, 0.6))
#'
#' # Permutation-based equivalence test
#' brunner_munzel(mpg ~ am, data = mtcars,
#'                alternative = "equivalence",
#'                mu = c(0.35, 0.65),
#'                test_method = "perm")
#'
#' @references
#' Brunner, E., Munzel, U. (2000). The Nonparametric Behrens-Fisher Problem: Asymptotic Theory and a Small Sample Approximation. Biometrical Journal 42, 17 -25.
#'
#' Neubert, K., Brunner, E., (2006). A Studentized Permutation Test for the Nonparametric Behrens-Fisher Problem. Computational Statistics and Data Analysis.
#'
#' Munzel, U., Brunner, E. (2002). An Exact Paired Rank Test. Biometrical Journal 44, 584-593.
#'
#' Munzel, U., Hauschke, D. (2003). A nonparametric test for proving noninferiority in clinical trials with ordered categorical data. Pharmaceutical Statistics, 2, 31-37.
#'
#' Konietschke, F., Placzek, M., Schaarschmidt, F., & Hothorn, L. A. (2015). nparcomp: an R software package for nonparametric multiple comparisons and simultaneous confidence intervals. Journal of Statistical Software 64 (2015), Nr. 9, 64(9), 1-17. http://www.jstatsoft.org/v64/i09/
#'
#' Janssen, A. (1997). Studentized permutation tests for non-i.i.d. hypotheses and the generalized Behrens-Fisher problem. Statistics & Probability Letters, 36(1), 9-21.
#'
#' Phipson, B., & Smyth, G. K. (2010). Permutation P-values should never be zero: calculating exact P-values when permutations are randomly drawn. Statistical Applications in Genetics and Molecular Biology, 9(1), Article 39.
#' @name brunner_munzel
#' @importFrom stats var quantile
#' @family Robust tests
#' @export brunner_munzel
NULL

#' @noRd
bm_perm_indices <- function(N, n.x, R) {
  # Internal helper function to generate permutation indices without replacement
  # Returns a matrix where each column contains indices for group x
  # Sample permutations without replacement to avoid duplicates
  perms <- matrix(NA, nrow = R, ncol = N)
  for (i in 1:R) {
    perms[i, ] <- sample(1:N, N)
  }
  # Remove duplicate permutations
  perms <- unique(perms)
  iter <- 0
  max_iter <- 10
  while (nrow(perms) < R && iter < max_iter) {
    needed <- R - nrow(perms)
    new_perms <- matrix(NA, nrow = needed * 2, ncol = N)
    for (i in 1:(needed * 2)) {
      new_perms[i, ] <- sample(1:N, N)
    }
    perms <- unique(rbind(perms, new_perms))
    iter <- iter + 1
  }
  if (nrow(perms) > R) {
    perms <- perms[1:R, , drop = FALSE]
  }
  return(t(perms))  # Return as N x R matrix for compatibility with original code
}

#' @noRd
bm_compute_perm_pval <- function(b, R, p_method) {
  # Internal helper function to compute permutation p-value
  # Based on Phipson & Smyth (2010) for sampling without replacement
  if (p_method == "plusone") {
    # Phipson & Smyth (2010) formula: (b+1)/(R+1)
    # Provides exact p-values when permutations are sampled without replacement
    return((b + 1) / (R + 1))
  } else {
    # Original formula: b/R with minimum 1/R
    pval <- b / R
    return(max(pval, 1 / R))
  }
}

#brunner_munzel <- setClass("brunner_munzel")
brunner_munzel <- function(x,
                           ...,
                           paired = FALSE,
                           alternative = c("two.sided",
                                           "less",
                                           "greater",
                                           "equivalence",
                                           "minimal.effect"),
                           mu = 0.5,
                           alpha = 0.05,
                           test_method = c("t", "logit", "perm"),
                           R = 10000,
                           p_method = c("plusone", "original"),
                           perm = "deprecated",
                           max_n_perm = "deprecated") {

  UseMethod("brunner_munzel")
}

#' @rdname brunner_munzel
#' @importFrom stats sd cor na.omit setNames t.test terms nlm optim optimize
#' @method brunner_munzel default
#' @export

# @method brunner_munzel default
brunner_munzel.default = function(x,
                                  y,
                                  paired = FALSE,
                                  alternative = c("two.sided",
                                                  "less",
                                                  "greater",
                                                  "equivalence",
                                                  "minimal.effect"),
                                  mu = 0.5,
                                  alpha = 0.05,
                                  test_method = c("t", "logit", "perm"),
                                  R = 10000,
                                  p_method = c("plusone", "original"),
                                  perm = "deprecated",
                                  max_n_perm = "deprecated",
                                  ...) {

  # Handle deprecated arguments
  if (!identical(perm, "deprecated")) {
    lifecycle::deprecate_warn(
      when = "0.9.0",
      what = "brunner_munzel(perm)",
      with = "brunner_munzel(test_method)",
      details = "Use `test_method = 'perm'` instead of `perm = TRUE`."
    )
    if (isTRUE(perm)) {
      test_method <- "perm"
    }
  }

  if (!identical(max_n_perm, "deprecated")) {
    lifecycle::deprecate_warn(
      when = "0.9.0",
      what = "brunner_munzel(max_n_perm)",
      with = "brunner_munzel(R)",
      details = "Use `R` to specify the number of permutations."
    )
    R <- max_n_perm
  }

  alternative = match.arg(alternative)
  test_method = match.arg(test_method)
  p_method = match.arg(p_method)

  # Validate mu based on alternative

  if(alternative %in% c("equivalence", "minimal.effect")) {
    if(length(mu) != 2) {
      stop("'mu' must be a vector of length 2 for equivalence or minimal.effect alternatives")
    }
    if(any(!is.finite(mu)) || any(mu < 0) || any(mu > 1)) {
      stop("'mu' values must be finite numbers between 0 and 1")
    }
    if(mu[1] >= mu[2]) {
      stop("'mu[1]' (lower bound) must be less than 'mu[2]' (upper bound)")
    }
    low_eqbound <- min(mu)
    high_eqbound <- max(mu)
  } else {
    if(length(mu) != 1 || !is.finite(mu) || mu < 0 || mu > 1) {
      stop("'mu' must be a single number between 0 and 1")
    }
    if(alternative == "two.sided") {
      if(mu == 0 || mu == 1) {
        stop("'mu' cannot be 0 or 1 when alternative is 'two.sided'")
      }
    }
  }

  if(!is.numeric(x)) stop("'x' must be numeric")
  if(!missing(y)) {
    if(!is.numeric(y)) stop("'y' must be numeric")
    DNAME <- paste(deparse(substitute(x)), "and",
                   deparse(substitute(y)))
    if(paired) {
      if(length(x) != length(y))
        stop("'x' and 'y' must have the same length")
      ok <- complete.cases(x,y)
      x = x[ok]
      y = y[ok]
    }
    else {
      x <- x[is.finite(x)]
      y <- y[is.finite(y)]
    }

    if(min(length(x),length(y)) < 15 && test_method != "perm"){
      message("Sample size in at least one group is small. Permutation test (test_method = 'perm') is highly recommended.")
    }
    if(min(length(x),length(y)) > 250 &&
       sum(length(x),length(y)) > 500 &&
       test_method == "perm"){
      message("Sample size is fairly large. Use of a permutation test is probably unnecessary.")
    }
  } else {

    stop("'y' is missing. One sample tests currently not supported.")

  }

  # Set confidence level based on alternative
  if(alternative %in% c("equivalence", "minimal.effect")) {
    conf.level <- 1 - alpha * 2
  } else {
    conf.level <- 1 - alpha
  }

  # Paired -----
  if(paired){

    n = length(x)
    n1 = n+1
    df.sw = n - 1
    all_data <- c(y, x)
    N = length(all_data)

    xinverse <- c(x, y)
    x1 <- y
    x2 <- x
    rx <- rank(all_data)
    rxinverse <- rank(xinverse)
    rx1 <- rx[1:n]
    rx2 <- rx[(n+1):N]
    rix1 <- rank(x1)
    rix2 <- rank(x2)
    BM1 <- 1 / n * (rx1 - rix1)
    BM2 <- 1 / n * (rx2 - rix2)
    BM3 <- BM1 - BM2
    BM4 <- 1 / (2 * n) * (rx1 - rx2)
    pd <- mean(BM2)

    m <- mean(BM3)
    v <- (sum(BM3 ^ 2) - n * m ^ 2) / (n - 1)
    v0 <- (v == 0)
    v[v0] <- 1 / n
    std_err = sqrt(v/n)

    if(test_method == "perm"){
      METHOD = "Paired Brunner-Munzel permutation test"
      if(alternative %in% c("equivalence", "minimal.effect")) {
        message("NOTE: Permutation-based TOST for equivalence/minimal.effect testing.")
      }

      # Directly from nparcomp
      if(n<=13){
        n_perm_actual <- 2^n
        p<-0
        for (i in 1:n){
          a<-rep(c(rep(c(i,i+n),n_perm_actual/(2^i)),rep(c(i+n,i),n_perm_actual/(2^i))),2^(i-1))
          p<-rbind(p,a)
        }
        p<-p[2:(n+1),]
        P<-matrix(p,ncol=n_perm_actual)

        xperm<-matrix(all_data[P],nrow=N,ncol=n_perm_actual)
        rxperm<-matrix(rx[P],nrow=N,ncol=n_perm_actual)
      }
      else{
        n_perm_actual <- R
        P<-matrix(nrow=n,ncol=n_perm_actual)
        permu<-function(all_data){
          n<-length(all_data)
          result<-sample(c(0,1),size=n,replace=TRUE)
          return(result)
        }
        P1<-apply(P,2,permu)
        P2<-rbind(P1,P1)
        xperm<-all_data*P2+xinverse*(1-P2)
        rxperm<-rx*P2+rxinverse*(1-P2)
      }
      xperm1<-xperm[1:n,]
      xperm2<-xperm[n1:N,]
      rperm1<-rxperm[1:n,]
      rperm2<-rxperm[n1:N,]
      riperm1<-apply(xperm1,2,rank)
      riperm2<-apply(xperm2,2,rank)
      BMperm2<-1/n*(rperm2-riperm2)
      BMperm3<-1/n*(rperm1-riperm1)-BMperm2
      pdperm<-colMeans(BMperm2)
      mperm3<-colMeans(BMperm3)
      vperm3<-(colSums(BMperm3^2)-n*mperm3^2)/(n-1)
      vperm30<-(vperm3==0)
      vperm3[vperm30]<-1/n

      # FIXED: Permuted statistics always centered at 0.5 (exchangeability assumption)
      Tperm <- sqrt(n) * (pdperm - 0.5) / sqrt(vperm3)

      # Calculate p-values based on alternative
      if(alternative %in% c("equivalence", "minimal.effect")) {
        # Observed test statistics for each bound
        test_stat_low <- sqrt(n) * (pd - low_eqbound) / sqrt(v)
        test_stat_high <- sqrt(n) * (pd - high_eqbound) / sqrt(v)

        # One-sided p-values
        # For lower bound test: H0: p <= low vs H1: p > low
        # Large test_stat_low supports H1, so p-value = P(T >= t_obs)
        p_greater_low <- mean(Tperm >= test_stat_low)

        # For upper bound test: H0: p >= high vs H1: p < high
        # Small (negative) test_stat_high supports H1, so p-value = P(T <= t_obs)
        p_less_high <- mean(Tperm <= test_stat_high)

        if(alternative == "equivalence") {
          # Both conditions must be met: p > low AND p < high
          # p-value is the maximum of the two one-sided tests
          p.value <- max(p_greater_low, p_less_high)

          # Determine which bound is "binding" for reporting
          if(p_greater_low >= p_less_high) {
            test_stat <- test_stat_low
          } else {
            test_stat <- test_stat_high
          }
        } else { # minimal.effect
          # At least one condition must be met: p <= low OR p >= high
          # p-value is the minimum of the two one-sided tests
          p.value <- min(1 - p_greater_low, 1 - p_less_high)

          # Determine which bound is "binding" for reporting
          if((1 - p_greater_low) <= (1 - p_less_high)) {
            test_stat <- test_stat_low
          } else {
            test_stat <- test_stat_high
          }
        }

        # Confidence interval quantiles from permutation distribution
        # Use actual number of permutations for indexing
        pq1 <- sort(Tperm)[(floor((1-alpha)*n_perm_actual)+1)]

        pd.lower <- pd - pq1*sqrt(v/n)
        pd.upper <- pd + pq1*sqrt(v/n)

      } else {
        # Standard alternatives (two.sided, less, greater)
        # Observed test statistic centered at mu
        test_stat <- sqrt(n) * (pd - mu) / sqrt(v)

        p1perm <- mean(Tperm <= test_stat)
        pq1 <- sort(Tperm)[(floor((1-alpha/2)*n_perm_actual)+1)]
        pq2 <- sort(Tperm)[(floor((1-alpha)*n_perm_actual)+1)]

        p.value = switch(alternative,
                         "two.sided" = min(2*p1perm, 2*(1-p1perm)),
                         "less" = p1perm,
                         "greater" = 1-p1perm)

        pd.lower = switch(alternative,
                          "two.sided" = pd - pq1*sqrt(v/n),
                          "less" = 0,
                          "greater" = pd - pq2*sqrt(v/n))

        pd.upper = switch(alternative,
                          "two.sided" = pd + pq1*sqrt(v/n),
                          "less" = pd + pq2*sqrt(v/n),
                          "greater" = 1)
      }

      # Clamp bounds to [0, 1]
      pd.lower <- ifelse(pd.lower < 0, 0, pd.lower)
      pd.upper <- ifelse(pd.upper > 1, 1, pd.upper)

    } else if(test_method == "logit") {
      # Logit transformation method for paired samples
      METHOD = "Exact paired Brunner-Munzel test (logit)"

      # Logit transformation for range-preserving CIs
      pd_logit <- log(pd / (1 - pd))
      se_logit <- sqrt(v/n) / (pd * (1 - pd))

      if(alternative %in% c("equivalence", "minimal.effect")) {
        # Test statistics for each bound on logit scale
        low_logit <- log(low_eqbound / (1 - low_eqbound))
        high_logit <- log(high_eqbound / (1 - high_eqbound))

        test_stat_low <- (pd_logit - low_logit) / se_logit
        test_stat_high <- (pd_logit - high_logit) / se_logit

        # One-sided p-values
        p_low_greater <- 1 - pt(test_stat_low, df.sw)
        p_high_less <- pt(test_stat_high, df.sw)

        if(alternative == "equivalence") {
          p.value <- max(p_low_greater, p_high_less)

          if(p_low_greater >= p_high_less) {
            test_stat <- test_stat_low
          } else {
            test_stat <- test_stat_high
          }
        } else { # minimal.effect
          p.value <- min(1 - p_low_greater, 1 - p_high_less)

          if((1 - p_low_greater) <= (1 - p_high_less)) {
            test_stat <- test_stat_low
          } else {
            test_stat <- test_stat_high
          }
        }

        # CI on logit scale, then back-transform
        ci_logit_lower <- pd_logit - qt(1-alpha, df.sw) * se_logit
        ci_logit_upper <- pd_logit + qt(1-alpha, df.sw) * se_logit

        pd.lower <- exp(ci_logit_lower) / (1 + exp(ci_logit_lower))
        pd.upper <- exp(ci_logit_upper) / (1 + exp(ci_logit_upper))

      } else {
        # Standard alternatives
        mu_logit <- log(mu / (1 - mu))
        test_stat <- (pd_logit - mu_logit) / se_logit

        p.value = switch(alternative,
                         "two.sided" = 2*min(pt(test_stat, df.sw),
                                             1-pt(test_stat, df.sw)),
                         "less" = pt(test_stat, df=df.sw),
                         "greater" = 1-pt(test_stat, df=df.sw))

        ci_logit_lower = switch(alternative,
                                "two.sided" = pd_logit - qt(1-alpha/2, df.sw) * se_logit,
                                "less" = -Inf,
                                "greater" = pd_logit - qt(1-alpha, df.sw) * se_logit)

        ci_logit_upper = switch(alternative,
                                "two.sided" = pd_logit + qt(1-alpha/2, df.sw) * se_logit,
                                "less" = pd_logit + qt(1-alpha, df.sw) * se_logit,
                                "greater" = Inf)

        pd.lower <- exp(ci_logit_lower) / (1 + exp(ci_logit_lower))
        pd.upper <- exp(ci_logit_upper) / (1 + exp(ci_logit_upper))
      }

      # Handle edge cases
      pd.lower <- ifelse(is.nan(pd.lower) | pd.lower < 0, 0, pd.lower)
      pd.upper <- ifelse(is.nan(pd.upper) | pd.upper > 1, 1, pd.upper)

    } else {
      # Asymptotic (t-distribution) approach
      METHOD = "Exact paired Brunner-Munzel test"

      if(alternative %in% c("equivalence", "minimal.effect")) {
        # Test statistics for each bound
        test_stat_low <- sqrt(n) * (pd - low_eqbound) / sqrt(v)
        test_stat_high <- sqrt(n) * (pd - high_eqbound) / sqrt(v)

        # One-sided p-values
        p_low_greater <- 1 - pt(test_stat_low, df.sw)  # p > low_eqbound
        p_high_less <- pt(test_stat_high, df.sw)       # p < high_eqbound

        if(alternative == "equivalence") {
          p.value <- max(p_low_greater, p_high_less)

          # Report the binding statistic
          if(p_low_greater >= p_high_less) {
            test_stat <- test_stat_low
          } else {
            test_stat <- test_stat_high
          }
        } else { # minimal.effect
          p.value <- min(1 - p_low_greater, 1 - p_high_less)

          # Report the binding statistic
          if((1 - p_low_greater) <= (1 - p_high_less)) {
            test_stat <- test_stat_low
          } else {
            test_stat <- test_stat_high
          }
        }

        # 1-2*alpha CI for TOST
        pd.lower <- pd - qt(1-alpha, df.sw)*sqrt(v/n)
        pd.upper <- pd + qt(1-alpha, df.sw)*sqrt(v/n)

      } else {
        # Standard alternatives
        test_stat <- sqrt(n) * (pd - mu) / sqrt(v)

        p.value = switch(alternative,
                         "two.sided" = 2*min(pt(test_stat, df.sw),
                                             1-pt(test_stat, df.sw)),
                         "less" = pt(test_stat, df=df.sw),
                         "greater" = 1-pt(test_stat, df=df.sw))

        pd.lower = switch(alternative,
                          "two.sided" = pd - qt(1-alpha/2, df.sw)*sqrt(v/n),
                          "less" = 0,
                          "greater" = pd - qt(1-alpha, df.sw)*sqrt(v/n))

        pd.upper = switch(alternative,
                          "two.sided" = pd + qt(1-alpha/2, df.sw)*sqrt(v/n),
                          "less" = pd + qt(1-alpha, df.sw)*sqrt(v/n),
                          "greater" = 1)
      }
    }

  } else {

    # Two-sample ------
    rxy <- rank(c(x, y))
    rx <- rank(x)
    ry <- rank(y)
    n.x <- as.double(length(x))
    n.y <- as.double(length(y))
    N = n.x + n.y

    pl2 <- 1/n.y*(rxy[1:n.x]-rx)
    pl1 <- 1/n.x*(rxy[(n.x+1):N]-ry)
    pd <- mean(pl2)
    pd1 <- (pd == 1)
    pd0 <- (pd == 0)
    pd[pd1] <- 0.9999
    pd[pd0] <- 0.0001
    s1 <- var(pl2)/n.x
    s2 <- var(pl1)/n.y

    V <- N*(s1+s2)
    singular.bf <- (V == 0)
    V[singular.bf] <- N/(2 * n.x * n.y)
    std_err = sqrt(V/N)

    df.sw <- (s1 + s2)^2/(s1^2/(n.x - 1) + s2^2/(n.y - 1))
    df.sw[is.nan(df.sw)] <- 1000

    if(test_method == "perm"){

      ## permutation -----
      METHOD = "Two-sample Brunner-Munzel permutation test"
      if(alternative %in% c("equivalence", "minimal.effect")) {
        message("NOTE: Permutation-based TOST for equivalence/minimal.effect testing.")
      }

      Tprob<-qnorm(pd)*exp(-0.5*qnorm(pd)^2)*sqrt(N/(V*2*pi))
      # Use sampling without replacement to avoid duplicate permutations
      P <- bm_perm_indices(N, n.x, R)
      R_actual <- ncol(P)  # Actual number of unique permutations obtained
      Px<-matrix(c(x,y)[P],ncol=R_actual)

      # perm_loop already centers at 0.5 (see res1[1,]<-(pdP-1/2)/sqrt(vP))
      Tperm<-t(apply(perm_loop(x=Px[1:n.x,],y=Px[(n.x+1):N,],
                               n.x=n.x,n.y=n.y,R=R_actual),1,sort))

      if(alternative %in% c("equivalence", "minimal.effect")) {
        # Observed test statistics for each bound
        test_stat_low <- sqrt(N) * (pd - low_eqbound) / sqrt(V)
        test_stat_high <- sqrt(N) * (pd - high_eqbound) / sqrt(V)

        # Count extreme values for p-value calculation
        # For lower bound test (H1: p > low): count permutations with T >= t_obs_low
        # For upper bound test (H1: p < high): count permutations with T <= t_obs_high
        b_greater_low <- sum(Tperm[1,] >= test_stat_low)
        b_less_high <- sum(Tperm[1,] <= test_stat_high)

        # One-sided p-values using selected method
        # p_greater_low = P(T >= t_obs | H0) for testing H1: p > low
        # p_less_high = P(T <= t_obs | H0) for testing H1: p < high
        p_greater_low <- bm_compute_perm_pval(b_greater_low, R_actual, p_method)
        p_less_high <- bm_compute_perm_pval(b_less_high, R_actual, p_method)

        if(alternative == "equivalence") {
          # Both conditions must be met: p > low AND p < high
          # p-value is the maximum of the two one-sided tests
          p.value <- max(p_greater_low, p_less_high)

          # Determine which bound is "binding" for reporting
          if(p_greater_low >= p_less_high) {
            test_stat <- test_stat_low
          } else {
            test_stat <- test_stat_high
          }
        } else { # minimal.effect
          # At least one condition must be met: p <= low OR p >= high
          # p-value is the minimum of the two one-sided tests
          p.value <- min(1 - p_greater_low, 1 - p_less_high)

          if((1 - p_greater_low) <= (1 - p_less_high)) {
            test_stat <- test_stat_low
          } else {
            test_stat <- test_stat_high
          }
        }

        # CI quantiles for 1-2*alpha level
        c1 <- 0.5*(Tperm[1, floor((1-alpha)*R_actual)] +
                     Tperm[1, ceiling((1-alpha)*R_actual)])

        pd.lower <- pd - sqrt(V/N)*c1
        pd.upper <- pd + sqrt(V/N)*c1

      } else {
        # Standard alternatives
        # Observed test statistic centered at mu
        test_stat <- sqrt(N) * (pd - mu) / sqrt(V)

        # Count extreme values for p-value calculation
        b_less <- sum(test_stat <= Tperm[1,])
        b_greater <- sum(test_stat >= Tperm[1,])

        if(alternative == "two.sided"){
          c1<-0.5*(Tperm[1,floor((1-alpha/2)*R_actual)]+Tperm[1,ceiling((1-alpha/2)*R_actual)])
          c2<-0.5*(Tperm[1,floor(alpha/2*R_actual)]+Tperm[1,ceiling(alpha/2*R_actual)])
        } else {
          c1<-0.5*(Tperm[1, floor((1-alpha)*R_actual)]+Tperm[1, ceiling((1-alpha)*R_actual)])
          c2<-0.5*(Tperm[1, floor(alpha*R_actual)]+Tperm[1, ceiling(alpha*R_actual)])
        }

        lower_ci = pd - sqrt(V/N)*c1
        upper_ci = pd - sqrt(V/N)*c2

        # Compute p-values using selected method
        p_less <- bm_compute_perm_pval(b_less, R_actual, p_method)
        p_greater <- bm_compute_perm_pval(b_greater, R_actual, p_method)

        p.value = switch(alternative,
                         "two.sided" = min(2 * p_less, 2 * p_greater),
                         "less" = p_less,
                         "greater" = p_greater)

        pd.lower = switch(alternative,
                          "two.sided" = lower_ci,
                          "less" = 0,
                          "greater" = lower_ci)

        pd.upper = switch(alternative,
                          "two.sided" = upper_ci,
                          "less" = upper_ci,
                          "greater" = 1)
      }

      pd.lower = ifelse(pd.lower < 0, 0, pd.lower)
      pd.upper = ifelse(pd.upper > 1, 1, pd.upper)

    } else if(test_method == "logit") {

      ## logit transformation ----
      METHOD = "Two-sample Brunner-Munzel test (logit)"

      # Logit transformation for range-preserving CIs
      pd_logit <- log(pd / (1 - pd))
      se_logit <- sqrt(V/N) / (pd * (1 - pd))

      if(alternative %in% c("equivalence", "minimal.effect")) {
        # Test statistics for each bound on logit scale
        low_logit <- log(low_eqbound / (1 - low_eqbound))
        high_logit <- log(high_eqbound / (1 - high_eqbound))

        test_stat_low <- (pd_logit - low_logit) / se_logit
        test_stat_high <- (pd_logit - high_logit) / se_logit

        # One-sided p-values
        p_low_greater <- 1 - pt(test_stat_low, df=df.sw)
        p_high_less <- pt(test_stat_high, df=df.sw)

        if(alternative == "equivalence") {
          p.value <- max(p_low_greater, p_high_less)

          if(p_low_greater >= p_high_less) {
            test_stat <- test_stat_low
          } else {
            test_stat <- test_stat_high
          }
        } else { # minimal.effect
          p.value <- min(1 - p_low_greater, 1 - p_high_less)

          if((1 - p_low_greater) <= (1 - p_high_less)) {
            test_stat <- test_stat_low
          } else {
            test_stat <- test_stat_high
          }
        }

        # CI on logit scale, then back-transform
        ci_logit_lower <- pd_logit - qt(1-alpha, df=df.sw) * se_logit
        ci_logit_upper <- pd_logit + qt(1-alpha, df=df.sw) * se_logit

        pd.lower <- exp(ci_logit_lower) / (1 + exp(ci_logit_lower))
        pd.upper <- exp(ci_logit_upper) / (1 + exp(ci_logit_upper))

      } else {
        # Standard alternatives
        mu_logit <- log(mu / (1 - mu))
        test_stat <- (pd_logit - mu_logit) / se_logit

        p.value = switch(alternative,
                         "two.sided" = min(c(2 - 2 * pt(test_stat, df=df.sw),
                                             2 * pt(test_stat, df=df.sw))),
                         "less" = pt(test_stat, df=df.sw),
                         "greater" = 1-pt(test_stat, df=df.sw))

        ci_logit_lower = switch(alternative,
                                "two.sided" = pd_logit - qt(1-alpha/2, df=df.sw) * se_logit,
                                "less" = -Inf,
                                "greater" = pd_logit - qt(1-alpha, df=df.sw) * se_logit)

        ci_logit_upper = switch(alternative,
                                "two.sided" = pd_logit + qt(1-alpha/2, df=df.sw) * se_logit,
                                "less" = pd_logit + qt(1-alpha, df=df.sw) * se_logit,
                                "greater" = Inf)

        pd.lower <- exp(ci_logit_lower) / (1 + exp(ci_logit_lower))
        pd.upper <- exp(ci_logit_upper) / (1 + exp(ci_logit_upper))
      }

      # Handle edge cases
      pd.lower <- ifelse(is.nan(pd.lower) | pd.lower < 0, 0, pd.lower)
      pd.upper <- ifelse(is.nan(pd.upper) | pd.upper > 1, 1, pd.upper)

    } else{

      ## t-approx ----
      METHOD = "Two-sample Brunner-Munzel test"

      if(alternative %in% c("equivalence", "minimal.effect")) {
        # Test statistics for each bound
        test_stat_low <- sqrt(N) * (pd - low_eqbound) / sqrt(V)
        test_stat_high <- sqrt(N) * (pd - high_eqbound) / sqrt(V)

        # One-sided p-values
        p_low_greater <- 1 - pt(test_stat_low, df=df.sw)  # p > low_eqbound
        p_high_less <- pt(test_stat_high, df=df.sw)       # p < high_eqbound

        if(alternative == "equivalence") {
          p.value <- max(p_low_greater, p_high_less)

          if(p_low_greater >= p_high_less) {
            test_stat <- test_stat_low
          } else {
            test_stat <- test_stat_high
          }
        } else { # minimal.effect
          p.value <- min(1 - p_low_greater, 1 - p_high_less)

          if((1 - p_low_greater) <= (1 - p_high_less)) {
            test_stat <- test_stat_low
          } else {
            test_stat <- test_stat_high
          }
        }

        # 1-2*alpha CI for TOST
        pd.lower <- pd - qt(1-alpha, df=df.sw)/sqrt(N)*sqrt(V)
        pd.upper <- pd + qt(1-alpha, df=df.sw)/sqrt(N)*sqrt(V)
        pd.lower = ifelse(pd.lower < 0, 0, pd.lower)
        pd.upper = ifelse(pd.upper > 1, 1, pd.upper)

      } else {
        # Standard alternatives
        test_stat <- sqrt(N) * (pd - mu) / sqrt(V)

        p.value = switch(alternative,
                         "two.sided" = min(c(2 - 2 * pt(test_stat, df=df.sw),
                                             2 * pt(test_stat, df=df.sw))),
                         "less" = pt(test_stat, df=df.sw),
                         "greater" = 1-pt(test_stat, df=df.sw))

        pd.lower = switch(alternative,
                          "two.sided" = pd - qt(1-alpha/2, df=df.sw)/sqrt(N)*sqrt(V),
                          "less" = 0,
                          "greater" = pd - qt(1-alpha, df=df.sw)/sqrt(N)*sqrt(V))
        pd.lower = ifelse(pd.lower < 0, 0, pd.lower)

        pd.upper = switch(alternative,
                          "two.sided" = pd + qt(1-alpha/2, df=df.sw)/sqrt(N)*sqrt(V),
                          "less" = pd + qt(1-alpha, df=df.sw)/sqrt(N)*sqrt(V),
                          "greater" = 1)
        pd.upper = ifelse(pd.upper > 1, 1, pd.upper)
      }
    }
  }

  # Prepare output
  if(alternative %in% c("equivalence", "minimal.effect")) {
    names(mu) <- c("lower bound", "upper bound")
  } else {
    names(mu) <- "relative effect"
  }

  if(test_method == "perm"){
    names(test_stat) = "t-observed"
  } else {
    names(test_stat) = "t"
  }

  names(df.sw) = "df"
  cint = c(pd.lower, pd.upper)
  attr(cint, "conf.level") = conf.level
  estimate = pd
  names(estimate) = "P(X>Y) + .5*P(X=Y)"

  rval <- list(statistic = test_stat,
               parameter = df.sw,
               p.value = as.numeric(p.value),
               estimate = estimate,
               stderr = std_err,
               conf.int = cint,
               null.value = mu,
               alternative = alternative,
               method = METHOD,
               data.name = DNAME)
  class(rval) <- "htest"
  return(rval)

}

#' @rdname brunner_munzel
#' @method brunner_munzel formula
#' @export

brunner_munzel.formula = function(formula,
                                  data,
                                  subset,
                                  na.action, ...) {

  if(missing(formula)
     || (length(formula) != 3L)
     || (length(attr(terms(formula[-2L]), "term.labels")) != 1L))
    stop("'formula' missing or incorrect")
  
  # Check for paired argument in ... and warn user
  dots <- list(...)
  if("paired" %in% names(dots)){
    if(isTRUE(dots$paired)){
      message("Using 'paired = TRUE' with the formula interface is not recommended. Please ensure your data is sorted appropriately to make the correct paired comparison.")
    }
  }

  m <- match.call(expand.dots = FALSE)
  if(is.matrix(eval(m$data, parent.frame())))
    m$data <- as.data.frame(data)
  ## need stats:: for non-standard evaluation
  m[[1L]] <- quote(stats::model.frame)
  m$... <- NULL
  mf <- eval(m, parent.frame())
  DNAME <- paste(names(mf), collapse = " by ")
  names(mf) <- NULL
  response <- attr(attr(mf, "terms"), "response")
  g <- factor(mf[[-response]])
  if(nlevels(g) != 2L)
    stop("grouping factor must have exactly 2 levels")
  DATA <- setNames(split(mf[[response]], g), c("x", "y"))
  y <- do.call("brunner_munzel", c(DATA, list(...)))
  y$data.name <- DNAME
  y

}


perm_loop <-function(x,y,n.x,n.y,R){

  pl1P<-matrix(0,nrow=n.x,ncol=R)
  pl2P<-matrix(0,nrow=n.y,ncol=R)

  for(h1 in 1:n.x){
    help1<-matrix(t(x[h1,]),
                  ncol=R,
                  nrow=n.y,byrow=TRUE)
    pl1P[h1,]<-1/n.y*(colSums((y<help1)+1/2*(y==help1)))
  }
  for(h2 in 1:n.y){
    help2<-matrix(t(y[h2,]),
                  ncol=R,
                  nrow=n.x,byrow=TRUE)
    pl2P[h2,]<-1/n.x*(colSums((x<help2)+1/2*(x==help2)))
  }

  pdP<-colMeans(pl2P)
  pd2P<-colMeans(pl1P)

  v1P<-(colSums(pl1P^2)-n.x*pd2P^2)/(n.x-1)
  v2P<-(colSums(pl2P^2)-n.y*pdP^2)/(n.y-1)
  vP<-v1P/n.x + v2P/n.y

  v0P<-(vP==0)
  vP[v0P]<-0.5/(n.x*n.y)^2

  res1<-matrix(rep(0,R*3),nrow=3)

  # Note: This centers at 0.5, which is correct for the permutation distribution
  res1[1,]<-(pdP-1/2)/sqrt(vP)

  pdP0<-(pdP==0)
  pdP1<-(pdP==1)
  pdP[pdP0]<-0.01
  pdP[pdP1]<-0.99

  res1[2,]<-log(pdP/(1-pdP))*pdP*(1-pdP)/sqrt(vP)
  res1[3,]<-qnorm(pdP)*exp(-0.5*qnorm(pdP)^2)/(sqrt(2*pi*vP))

  res1
}
