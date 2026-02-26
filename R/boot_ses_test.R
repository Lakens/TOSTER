#' @title Parametric Bootstrap Test for Rank-Based Effect Sizes
#' @description
#' `r lifecycle::badge('experimental')`
#'
#' Performs hypothesis testing for rank-based effect sizes using a parametric
#' bootstrap. This function is designed primarily for equivalence testing (TOST)
#' and minimal effect testing with non-zero null hypotheses, where
#' permutation-based approaches are not valid for rank-based effect sizes.
#'
#' @section Warning:
#' This function is experimental. Important caveats:
#' - Validity depends on parametric assumptions (Lehmann alternative for
#'   two-sample, sign-randomization for paired). Unlike permutation tests,
#'   this is not assumption-free.
#' - The procedure is most reliable for continuous data, moderate n (>= 20),
#'   and equivalence bounds not too close to +/-1.
#' - Results should be interpreted with caution and ideally cross-checked
#'   against asymptotic methods from [ses_calc()].
#' - This function exists because no assumption-free alternative for non-zero
#'   null TOST is available for rank-based effect sizes. The choice is between
#'   this and no test at all, not between this and a better test.
#'
#' @inheritParams ses_calc
#' @param x numeric vector of data values (first group or pre-treatment).
#' @param y numeric vector of data values (second group or post-treatment).
#' @param ses a character string specifying the effect size measure:
#'     - "rb": rank-biserial correlation (default)
#'     - "cstat": concordance statistic (C-statistic/AUC)
#'     - "odds": Wilcoxon-Mann-Whitney odds
#'     - "logodds": Wilcoxon-Mann-Whitney log-odds
#' @param alpha significance level (default = 0.05).
#' @param mu the null hypothesis value(s) on the scale of the chosen effect size.
#'     - For standard alternatives: a single value (default = 0 for rb/logodds,
#'       0.5 for cstat, 1 for odds)
#'     - For equivalence/minimal.effect: two values representing the lower and
#'       upper bounds, or a single value for symmetric bounds
#' @param alternative a character string specifying the alternative hypothesis:
#'     - "two.sided": Test whether effect differs from mu
#'     - "less": Test whether effect is less than mu
#'     - "greater": Test whether effect is greater than mu
#'     - "equivalence": TOST equivalence test (effect inside bounds)
#'     - "minimal.effect": Minimal effect test (effect outside bounds)
#' @param B integer; the number of bootstrap replicates (default = 2000).
#'     Increase to 5000+ for publication-quality results.
#' @param keep_boot logical; if `TRUE` (default), return the bootstrap
#'     distributions in the output.
#' @param ... further arguments to be passed to or from methods.
#'
#' @details
#' This function calculates p-values for rank-based effect sizes using a
#' parametric bootstrap. It generates data under the null hypothesis and
#' compares the observed effect size to the resulting null reference
#' distribution.
#'
#' ## Why Not Permutation?
#'
#' Permutation tests are exact and assumption-free when testing the null
#' \eqn{\rho = 0}. However, for non-zero nulls — as required by equivalence
#' (TOST) and minimal effect testing — the permutation distribution cannot be
#' shifted to the correct null by arithmetic operations. The rank-biserial
#' correlation is a nonlinear, bounded function of the data, so there is no
#' data transformation that shifts rb by a fixed \eqn{\Delta}. Studentization
#' (as used in [perm_t_test()] and [brunner_munzel()]) cannot rescue this
#' because rb is not a studentized statistic.
#'
#' This function uses parametric models to generate data under the non-zero
#' null. The tradeoff is that validity now depends on how well the model
#' approximates the true data-generating process.
#'
#' Users who need TOST for means should use [boot_t_TOST()], which handles
#' non-zero nulls correctly via studentization without any parametric
#' assumption.
#'
#' ## Why No Confidence Intervals?
#'
#' This function intentionally omits confidence intervals. The parametric
#' bootstrap here generates data under a specific null to produce a reference
#' distribution for p-value computation. This is fundamentally different from
#' the nonparametric bootstrap in [boot_ses_calc()], which resamples from the
#' observed data to characterize the sampling distribution of the estimator.
#' Users who need confidence intervals should use [boot_ses_calc()] or the
#' asymptotic intervals from [ses_calc()].
#'
#' ## Algorithm
#'
#' **Two-sample (independent)**: Uses the Lehmann alternative (proportional
#' hazards) model. Data are pooled and treated as a common empirical CDF
#' \eqn{F}. Under the null \eqn{P(X > Y) = \text{target}}, group Y is
#' resampled from \eqn{F} and group X is resampled from a transformed CDF
#' \eqn{G(t) = 1 - (1 - F(t))^{1/\gamma}} where
#' \eqn{\gamma = \text{target\_cstat} / (1 - \text{target\_cstat})}. This
#' produces data with the correct population rank-biserial under the
#' proportional hazards assumption.
#'
#' **Paired samples**: Uses a sign-randomization model. The absolute
#' differences \eqn{|d_i| = |x_i - y_i|} are resampled with replacement, and
#' signs are assigned independently with probability
#' \eqn{P(\text{sign} = +) = (1 + \Delta) / 2}. This produces a bootstrap
#' distribution with the correct expected signed-rank rb under the assumption
#' of rank-independent sign probabilities.
#'
#' For equivalence testing, two null distributions are generated (one per bound)
#' and the TOST p-value is the maximum of the two one-sided p-values.
#'
#' @section Future Work:
#' \itemize{
#'   \item **Confidence intervals**: Not currently provided. Adding CIs would
#'     require test inversion across a grid of null values (Berger & Boos, 1994),
#'     which is computationally expensive. Use [boot_ses_calc()] or [ses_calc()]
#'     for interval estimation.
#'   \item **Rank-dependent sign model**: The current paired bootstrap assigns
#'     signs independently of rank. A rank-weighted sign model could improve
#'     accuracy by allowing the sign probability to depend on the magnitude of
#'     the difference.
#'   \item **Copula-based generation**: A normal copula model could provide an
#'     alternative data generation mechanism, especially for paired data where
#'     the joint distribution matters beyond marginals. This would require
#'     numerical calibration of the copula parameter to the target rb.
#' }
#'
#' @return A list with class `"htest"` containing:
#'   \item{estimate}{Observed effect size on the requested scale.}
#'   \item{p.value}{Bootstrap p-value.}
#'   \item{alternative}{Character string describing the alternative hypothesis.}
#'   \item{method}{Description of the test performed.}
#'   \item{null.value}{Null hypothesis value(s) on the requested scale.}
#'   \item{data.name}{Character string giving the name(s) of the data.}
#'   \item{call}{The matched call.}
#'   \item{model.param}{The model parameter(s) used for null generation
#'     (for diagnostics). For two-sample: the Lehmann gamma parameter(s).
#'     For paired: the sign probability parameter(s).}
#'   \item{boot.dist}{Bootstrap null distribution (if `keep_boot = TRUE` and
#'     standard alternative).}
#'   \item{boot.dist.low}{Bootstrap distribution under lower bound (if
#'     `keep_boot = TRUE` and TOST).}
#'   \item{boot.dist.high}{Bootstrap distribution under upper bound (if
#'     `keep_boot = TRUE` and TOST).}
#'
#' @examples
#' \donttest{
#' # Example 1: Two-sided test
#' set.seed(42)
#' x <- rnorm(30, mean = 0)
#' y <- rnorm(30, mean = 0.5)
#' boot_ses_test(x = x, y = y, ses = "rb",
#'               mu = 0, alternative = "two.sided", B = 599)
#'
#' # Example 2: Equivalence test with rank-biserial
#' boot_ses_test(x = x, y = y, ses = "rb",
#'               mu = c(-0.4, 0.4), alternative = "equivalence", B = 599)
#'
#' # Example 3: Paired samples
#' pre  <- c(4.5, 5.2, 3.8, 6.1, 4.9, 5.7, 3.6, 5.0, 4.3, 6.5,
#'           4.1, 5.5, 3.9, 6.0, 4.7, 5.3, 3.7, 5.1, 4.4, 6.3)
#' post <- c(5.1, 4.9, 4.5, 5.8, 5.5, 5.2, 4.3, 5.4, 4.0, 6.2,
#'           4.8, 5.3, 4.2, 5.7, 5.1, 5.0, 4.1, 5.3, 4.2, 6.1)
#' boot_ses_test(x = pre, y = post, paired = TRUE,
#'               ses = "rb", mu = 0, alternative = "two.sided", B = 599)
#'
#' # Example 4: Using formula interface
#' data(mtcars)
#' boot_ses_test(formula = mpg ~ am, data = mtcars,
#'               ses = "rb", mu = 0,
#'               alternative = "two.sided", B = 599)
#' }
#'
#' @references
#' Berger, R.L. and Boos, D.D. (1994). P values maximized over a confidence set
#' for the nuisance parameter. *Journal of the American Statistical Association*,
#' 89, 1012-1016.
#'
#' Lehmann, E.L. (1975). *Nonparametrics: Statistical Methods Based on Ranks*.
#' Holden-Day.
#'
#' @family Robust tests
#' @seealso [ses_calc()] for asymptotic inference, [boot_ses_calc()] for
#'   bootstrap confidence intervals, [brunner_munzel()] with
#'   `test_method = "perm"` for robust TOST on the probability scale.
#' @name boot_ses_test
#' @export boot_ses_test

boot_ses_test <- function(x, ...,
                          paired = FALSE,
                          ses = c("rb", "cstat", "odds", "logodds"),
                          alpha = 0.05,
                          mu = NULL,
                          alternative = c("two.sided", "less", "greater",
                                          "equivalence", "minimal.effect"),
                          B = 2000L,
                          keep_boot = TRUE) {
  UseMethod("boot_ses_test")
}

#' @rdname boot_ses_test
#' @method boot_ses_test default
#' @export
boot_ses_test.default <- function(x,
                                  y = NULL,
                                  paired = FALSE,
                                  ses = c("rb", "cstat", "odds", "logodds"),
                                  alpha = 0.05,
                                  mu = NULL,
                                  alternative = c("two.sided", "less", "greater",
                                                  "equivalence", "minimal.effect"),
                                  B = 2000L,
                                  keep_boot = TRUE,
                                  ...) {

  ses <- match.arg(ses)
  alternative <- match.arg(alternative)

  # Validate inputs --------
  if (is.null(y)) {
    stop("boot_ses_test requires two samples (x and y). One-sample designs are not supported.")
  }

  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a numeric value between 0 and 1")
  }

  B <- as.integer(B)
  if (is.na(B) || B < 100L) {
    stop("B must be an integer >= 100")
  }

  # Set default mu based on ses scale
  if (is.null(mu)) {
    mu <- switch(ses,
                 "rb" = 0,
                 "cstat" = 0.5,
                 "odds" = 1,
                 "logodds" = 0)
  }

  # Handle equivalence/minimal.effect bounds --------
  if (alternative %in% c("equivalence", "minimal.effect")) {
    if (length(mu) == 1) {
      # Create symmetric bounds
      mu <- switch(ses,
                   "rb" = c(-abs(mu), abs(mu)),
                   "cstat" = c(0.5 - abs(mu - 0.5), 0.5 + abs(mu - 0.5)),
                   "odds" = c(1 / mu, mu),
                   "logodds" = c(-abs(mu), abs(mu)))
      # Ensure correct ordering
      mu <- sort(mu)
    }
    if (length(mu) != 2) {
      stop("For equivalence or minimal.effect testing, mu must be a single value (for symmetric bounds) or a vector of two values.")
    }
    low_bound <- min(mu)
    high_bound <- max(mu)
  } else {
    if (length(mu) > 1) {
      warning("mu has length > 1; only the first element will be used")
      mu <- mu[1]
    }
  }

  # Validate mu is within valid range for ses scale --------
  validate_mu_range <- function(val, ses_type) {
    switch(ses_type,
           "rb" = {
             if (any(val <= -1 | val >= 1))
               stop("mu must be in the open interval (-1, 1) for rb scale")
           },
           "cstat" = {
             if (any(val <= 0 | val >= 1))
               stop("mu must be in the open interval (0, 1) for cstat scale")
           },
           "odds" = {
             if (any(val <= 0))
               stop("mu must be positive for odds scale")
           },
           "logodds" = {
             # logodds is unbounded, no validation needed
           })
  }
  validate_mu_range(mu, ses)

  # Data name
  dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))

  # Handle NA and paired data --------
  if (paired) {
    data <- data.frame(x = x, y = y)
    data <- na.omit(data)
    x <- data$x
    y <- data$y
  } else {
    x <- na.omit(x)
    y <- na.omit(y)
  }

  n_x <- length(x)
  n_y <- length(y)
  n <- if (paired) n_x else min(n_x, n_y)

  # Small n warning --------
  if (n < 20) {
    warning(
      "Sample size is small (n = ", n, "). ",
      "The parametric bootstrap may not be reliable for small samples. ",
      "Results should be treated as exploratory."
    )
  }

  # Compute observed effect size --------
  obs_rb <- rbs_calc(x = x, y = y, mu = 0, paired = paired)

  # Complete separation warning --------
  if (abs(obs_rb) >= 0.999) {
    warning(
      "Complete or near-complete separation detected (rb = ", round(obs_rb, 3), "). ",
      "The parametric bootstrap result may be unreliable. ",
      "Interpret with extreme caution."
    )
  }

  # Scale conversion helpers --------
  ses_to_rb <- function(val, ses_type) {
    switch(ses_type,
           "rb" = val,
           "cstat" = cstat_to_rb(val),
           "odds" = (val - 1) / (val + 1),
           "logodds" = tanh(val / 2))
  }

  rb_to_ses_val <- function(rb_val, ses_type) {
    switch(ses_type,
           "rb" = rb_val,
           "cstat" = rb_to_cstat(rb_val),
           "odds" = rb_to_odds(rb_val),
           "logodds" = log(rb_to_odds(rb_val)))
  }

  obs_es <- rb_to_ses_val(obs_rb, ses)

  # Bootstrap null distribution generators --------

  # Two-sample: Lehmann alternative (proportional hazards) model
  # Under G(t) = 1 - (1-F(t))^(1/gamma), P(X > Y) = gamma/(1+gamma)
  # So gamma = target_cstat / (1 - target_cstat)
  boot_rb_twosample <- function(x, y, target_rb, B) {
    pooled <- c(x, y)
    n_x_local <- length(x)
    n_y_local <- length(y)
    target_cstat <- (1 + target_rb) / 2
    gamma_val <- target_cstat / (1 - target_cstat)

    boot_dist <- vapply(seq_len(B), function(b) {
      u_y <- runif(n_y_local)
      u_x <- runif(n_x_local)
      # Y resampled from pooled empirical CDF
      y_boot <- quantile(pooled, u_y, type = 1)
      # X from Lehmann-transformed CDF: quantile at 1 - (1-u)^gamma
      x_boot <- quantile(pooled, 1 - (1 - u_x)^gamma_val, type = 1)
      rbs_calc(x_boot, y_boot, mu = 0, paired = FALSE)
    }, numeric(1))

    list(dist = boot_dist, model_param = gamma_val)
  }

  # Paired: sign-randomization model
  # Resample |d_i| with replacement, assign sign + with P = (1 + target_rb) / 2
  boot_rb_paired <- function(x, y, target_rb, B) {
    d <- x - y
    abs_d <- abs(d)
    # Remove zero differences (consistent with signed-rank convention)
    abs_d <- abs_d[abs_d != 0]
    n_d <- length(abs_d)
    p_sign <- (1 + target_rb) / 2

    if (n_d == 0) {
      return(list(dist = rep(0, B), model_param = p_sign))
    }

    boot_dist <- vapply(seq_len(B), function(b) {
      # Resample absolute differences
      idx <- sample(n_d, replace = TRUE)
      abs_d_boot <- abs_d[idx]
      # Assign signs with target probability
      signs <- sample(c(1, -1), n_d, replace = TRUE, prob = c(p_sign, 1 - p_sign))
      d_boot <- abs_d_boot * signs
      # Compute signed-rank rb
      d_nz <- d_boot[d_boot != 0]
      n_nz <- length(d_nz)
      if (n_nz == 0) return(0)
      rr <- rank(abs(d_nz))
      T_plus <- sum(rr[d_nz > 0])
      S <- n_nz * (n_nz + 1) / 2
      2 * T_plus / S - 1
    }, numeric(1))

    list(dist = boot_dist, model_param = p_sign)
  }

  # Dispatch to appropriate bootstrap generator
  boot_rb_under_null <- function(x, y, target_rb, B, paired) {
    if (paired) {
      boot_rb_paired(x, y, target_rb, B)
    } else {
      boot_rb_twosample(x, y, target_rb, B)
    }
  }

  # Compute p-values --------
  if (alternative %in% c("equivalence", "minimal.effect")) {
    mu_rb_low <- ses_to_rb(low_bound, ses)
    mu_rb_high <- ses_to_rb(high_bound, ses)

    # Warn if bounds are near +/-1 on rb scale
    if (any(abs(c(mu_rb_low, mu_rb_high)) > 0.95)) {
      warning(
        "One or more equivalence bounds are close to +/-1 on the rb scale. ",
        "The parametric bootstrap approximation may be unreliable in this region."
      )
    }

    boot_low_res <- boot_rb_under_null(x, y, target_rb = mu_rb_low, B = B, paired = paired)
    boot_high_res <- boot_rb_under_null(x, y, target_rb = mu_rb_high, B = B, paired = paired)

    boot_low <- boot_low_res$dist
    boot_high <- boot_high_res$dist
    model_params <- c(boot_low_res$model_param, boot_high_res$model_param)
    names(model_params) <- c("lower", "upper")

    if (alternative == "equivalence") {
      # IU test: reject non-equivalence if inside bounds
      p_low <- mean(boot_low >= obs_rb)    # H1: rb > lower bound
      p_high <- mean(boot_high <= obs_rb)   # H1: rb < upper bound
      pvalue <- max(p_low, p_high)
    } else {
      # UI test: reject equivalence if outside bounds
      p_low <- mean(boot_low <= obs_rb)
      p_high <- mean(boot_high >= obs_rb)
      pvalue <- min(p_low, p_high)
    }

    # Null values on user scale
    null_val <- c(low_bound, high_bound)
    names(null_val) <- c("lower bound", "upper bound")

  } else {
    mu_rb <- ses_to_rb(mu, ses)

    # Warn if null is near +/-1 on rb scale
    if (abs(mu_rb) > 0.95) {
      warning(
        "The null value is close to +/-1 on the rb scale. ",
        "The parametric bootstrap approximation may be unreliable in this region."
      )
    }

    boot_res <- boot_rb_under_null(x, y, target_rb = mu_rb, B = B, paired = paired)
    boot_dist <- boot_res$dist
    model_params <- boot_res$model_param
    names(model_params) <- "null"

    pvalue <- switch(alternative,
                     "two.sided" = mean(abs(boot_dist - mu_rb) >= abs(obs_rb - mu_rb)),
                     "greater" = mean(boot_dist >= obs_rb),
                     "less" = mean(boot_dist <= obs_rb))

    # Null value on user scale
    ses_name_est <- switch(ses,
                           "rb" = "rank-biserial",
                           "cstat" = "concordance",
                           "odds" = "WMW odds",
                           "logodds" = "WMW log-odds")
    null_val <- mu
    names(null_val) <- ses_name_est
  }

  # Build output --------
  if (is.null(y)) {
    sample_type <- "One Sample"
  } else if (paired) {
    sample_type <- "Paired Sample"
  } else {
    sample_type <- "Two Sample"
  }

  ses_label <- switch(ses,
                      "rb" = "Rank-Biserial Correlation",
                      "cstat" = "Concordance",
                      "odds" = "WMW Odds",
                      "logodds" = "WMW Log-Odds")

  method_desc <- paste0("Parametric Bootstrap ", sample_type, " ", ses_label, " Test")

  estimate <- obs_es
  names(estimate) <- ses_label

  rval <- list(
    estimate = estimate,
    p.value = pvalue,
    alternative = alternative,
    method = method_desc,
    null.value = null_val,
    data.name = dname,
    call = match.call(),
    model.param = model_params
  )

  # Add bootstrap distributions if requested
  if (keep_boot) {
    if (alternative %in% c("equivalence", "minimal.effect")) {
      rval$boot.dist.low <- rb_to_ses_val(boot_low, ses)
      rval$boot.dist.high <- rb_to_ses_val(boot_high, ses)
    } else {
      rval$boot.dist <- rb_to_ses_val(boot_dist, ses)
    }
  }

  # TODO 1 — Confidence intervals: Not provided. Adding CIs would require test
  # inversion across a grid of null values (Berger & Boos, 1994), which is
  # computationally expensive (B x grid size resamples). Use boot_ses_calc() or
  # ses_calc() for interval estimation.

  # TODO 2 — Rank-dependent sign model for paired: The current sign-randomization
  # assigns signs independently of rank. A model where P(sign = +) depends on the
  # rank of |d_i| could improve accuracy for some data configurations.

  # TODO 3 — Copula-based generation: A normal copula model could provide an
  # alternative mechanism, especially for paired data. Would require numerical
  # calibration rather than the Spearman-to-copula formula (which applies to
  # Spearman correlation, not the rank-biserial).

  class(rval) <- "htest"
  return(rval)
}

#' @rdname boot_ses_test
#' @method boot_ses_test formula
#' @export
boot_ses_test.formula <- function(formula,
                                  data,
                                  subset,
                                  na.action, ...) {

  if (missing(formula)
      || (length(formula) != 3L)
      || (length(attr(terms(formula[-2L]), "term.labels")) != 1L))
    stop("'formula' missing or incorrect")

  # Check for paired argument in ... and warn user
  dots <- list(...)
  if ("paired" %in% names(dots)) {
    if (isTRUE(dots$paired)) {
      message("Using 'paired = TRUE' with the formula interface is not recommended. Please ensure your data is sorted appropriately to make the correct paired comparison.")
    }
  }

  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval(m$data, parent.frame())))
    m$data <- as.data.frame(data)
  m[[1L]] <- quote(stats::model.frame)
  m$... <- NULL
  mf <- eval(m, parent.frame())
  DNAME <- paste(names(mf), collapse = " by ")
  names(mf) <- NULL
  response <- attr(attr(mf, "terms"), "response")
  g <- factor(mf[[-response]])
  if (nlevels(g) != 2L)
    stop("grouping factor must have exactly 2 levels")
  DATA <- setNames(split(mf[[response]], g), c("x", "y"))
  y <- do.call("boot_ses_test", c(DATA, list(...)))
  y$data.name <- DNAME
  y
}
