#' @title Hodges-Lehmann Test
#' @description
#' `r lifecycle::badge('maturing')`
#'
#' Performs location tests based on the Hodges-Lehmann estimator with permutation-based
#' or asymptotic inference. This function supports one-sample, two-sample (independent),
#' and paired designs, as well as equivalence and minimal effect testing.
#'
#' @section Purpose:
#' The Hodges-Lehmann estimator provides a robust alternative to the mean for testing
#' location differences. It has a breakdown point of approximately 29.3%, meaning it
#' remains stable even when nearly 30% of the data are outliers. This function offers:
#'
#'   * Exact permutation tests for small samples
#'   * Randomization tests (permutation with replacement) for larger samples
#'   * Asymptotic tests using kernel density estimation
#'   * Support for equivalence and minimal effect testing
#'   * An interface that mirrors `wilcox.test` and `perm_t_test`
#'
#' @inheritParams perm_t_test
#' @param x a (non-empty) numeric vector of data values.
#' @param y an optional numeric vector of data values.
#' @param mu a number or vector specifying the null hypothesis value(s):
#'   * For standard alternatives (two.sided, less, greater): a single value (default: 0)
#'   * For equivalence/minimal.effect: either a single value (symmetric bounds will be
#'     created) or a vector of two values representing the lower and upper bounds
#' @param paired a logical indicating whether this is a paired test. If `TRUE`, `x` and
#'   `y` must have the same length, and differences (`x - y`) are analyzed.
#' @param alpha significance level (default = 0.05). For standard alternatives, the confidence
#'   interval has level 1-alpha. For equivalence/minimal.effect, the confidence interval
#'   has level 1-2*alpha (90% when alpha = 0.05).
#' @param R the number of permutations. Default is `NULL`, which uses the asymptotic test.
#'   If `R >= max_perms` (the maximum number of possible permutations), exact permutation
#'   is computed. Otherwise, Monte Carlo (permutation with replacement; i.e., randomization testing) sampling is used.
#' @param scale the scale estimator for standardizing the test statistic in permutation
#'   tests. Options are:
#'   * `"S2"` (default): Median of absolute pairwise differences from the median-corrected
#'     combined sample. Recommended when samples come from the same distribution under H0.
#'   * `"S1"`: Pooled within-sample absolute differences. Preferred when sample sizes are
#'     small or unequal.
#' @param p_method the method for computing permutation p-values. Options are:
#'   * `NULL` (default): Automatically selects "exact" for exact permutation tests and
#'     "plusone" for randomization tests.
#'   * `"exact"`: Uses b/R where b is the count of permutation statistics at least as
#'     extreme as observed. Appropriate when all permutations are enumerated.
#'   * `"plusone"`: Uses (b+1)/(R+1), which guarantees p > 0 and provides exact Type I
#'     error control for randomization tests where permutations are sampled with replacement
#'     (Phipson & Smyth, 2010).
#' @param keep_perm logical. If `TRUE` (default), the permutation distribution of the test
#'   statistic and effects are stored in the output. Set to `FALSE` for large datasets
#'   to save memory.
#' @param ... further arguments (currently ignored).
#'
#' @details
#' ## Hodges-Lehmann Estimators
#'
#' **One-sample/paired (HL1):** The median of all pairwise averages (Walsh averages):
#' \deqn{\hat{\theta}_{HL1} = \text{med}\left\{\frac{X_i + X_j}{2} : 1 \leq i \leq j \leq n\right\}}
#'
#' This estimator is consistent with the pseudomedian returned by [stats::wilcox.test()]
#' for one-sample and paired tests.
#'
#' **Two-sample (HL2):** The median of all pairwise differences between samples:
#' \deqn{\hat{\Delta}_{HL2} = \text{med}\{Y_j - X_i : i = 1, \ldots, m; j = 1, \ldots, n\}}
#'
#' This estimator is consistent with the location shift estimate returned by
#' [stats::wilcox.test()] for two-sample tests.
#'
#' ## Test Methods
#'
#' **Asymptotic test (R = NULL):** Uses kernel density estimation to estimate the
#' variance of the Hodges-Lehmann estimator (note: this generates confidence intervals that will differ from [stats::wilcox.test()]). The test statistic follows an approximate
#' normal distribution. This method may have issues with very heavy-tailed distributions,
#' very skewed distributions, or small sample sizes (n < 30 per group). In these cases,
#' consider using the permutation test instead.
#'
#' **Exact permutation test (R >= max_perms):** Enumerates all possible permutations
#' and provides exact p-values. For one-sample/paired tests, there are 2^n possible
#' sign-flipping permutations. For two-sample tests, there are choose(m+n, m) possible
#' group reassignments.
#'
#' **Randomization test (R < max_perms):** Samples R permutations randomly
#' (with replacement) and computes approximate p-values. The (b+1)/(R+1) formula
#' guarantees exact Type I error control (Phipson & Smyth, 2010).
#'
#' ## Scale Estimators
#'
#' For permutation tests, the test statistic is standardized using a robust scale
#' estimator (Fried & Dehling, 2011):
#'
#' **S1 (pooled within-sample):**
#' \deqn{S^{(1)}_{m,n} = \text{med}\{|X_i - X_j| : 1 \leq i < j \leq m, |Y_i - Y_j| : 1 \leq i < j \leq n\}}
#'
#' **S2 (median-corrected joint):**
#' \deqn{S^{(2)}_{m,n} = \text{med}\{|Z_i - Z_j| : 1 \leq i < j \leq m+n\}}
#'
#' where Z is the median-corrected combined sample.
#'
#' ## Alternatives
#'
#' The function supports five alternative hypotheses:
#' * `"two.sided"`: Tests whether the location parameter differs from `mu`
#' * `"less"`: Tests whether the location parameter is less than `mu`
#' * `"greater"`: Tests whether the location parameter is greater than `mu`
#' * `"equivalence"`: Tests whether the location parameter lies within the bounds
#'   specified by `mu` (TOST procedure)
#' * `"minimal.effect"`: Tests whether the location parameter lies outside the bounds
#'   specified by `mu`
#'
#' @return A list with class `"htest"` containing the following components:
#'
#'   - `statistic`: the value of the test statistic.
#'   - `p.value`: the p-value for the test.
#'   - `conf.int`: a confidence interval for the location parameter.
#'   - `estimate`: the Hodges-Lehmann estimate(s).
#'   - `null.value`: the specified hypothesized value(s).
#'   - `alternative`: a character string describing the alternative hypothesis.
#'   - `method`: a character string indicating the test type.
#'   - `data.name`: a character string giving the name(s) of the data.
#'   - `call`: the matched call.
#'   - `R`: the requested number of permutations (NULL for asymptotic).
#'   - `R.used`: the actual number of permutations used.
#'   - `perm.stat`: (if `keep_perm = TRUE`) the permutation distribution of test statistics.
#'   - `perm.eff`: (if `keep_perm = TRUE`) the permutation distribution of effects.
#'
#' @examples
#' # Two-sample test (asymptotic)
#' set.seed(123)
#' x <- rnorm(30, mean = 0)
#' y <- rnorm(30, mean = 0.5)
#' hodges_lehmann(x, y)
#'
#' # Two-sample test with permutation
#' hodges_lehmann(x, y, R = 1999)
#'
#' # One-sample test
#' x <- rnorm(20, mean = 0.5)
#' hodges_lehmann(x, mu = 0)
#'
#' # Paired test
#' before <- c(5.1, 4.8, 6.2, 5.7, 6.0, 5.5, 4.9, 5.8)
#' after <- c(5.6, 5.2, 6.7, 6.1, 6.5, 5.8, 5.3, 6.2)
#' hodges_lehmann(before, after, paired = TRUE)
#'
#' # Equivalence test
#' hodges_lehmann(x, y, alternative = "equivalence", mu = c(-1, 1), R = 999)
#'
#' # Formula interface
#' hodges_lehmann(extra ~ group, data = sleep)
#'
#' @references
#' Hodges, J. L., & Lehmann, E. L. (1963). Estimates of location based on rank tests.
#' *Annals of Mathematical Statistics*, 34, 598-611.
#'
#' Fried, R., & Dehling, H. (2011). Robust nonparametric tests for the two-sample
#' location problem. *Statistical Methods & Applications*, 20, 409-422.
#'
#' Lehmann, E. L. (1963). Nonparametric confidence intervals for a shift parameter.
#' *Annals of Mathematical Statistics*, 34, 1507-1512.
#'
#' Phipson, B., & Smyth, G. K. (2010). Permutation P-values should never be zero:
#' calculating exact P-values when permutations are randomly drawn.
#' *Statistical Applications in Genetics and Molecular Biology*, 9(1), Article 39.
#'
#' @family Robust tests
#' @name hodges_lehmann
#' @export hodges_lehmann

hodges_lehmann <- function(x, ...) {
  UseMethod("hodges_lehmann")
}

# =============================================================================
# Internal helper functions
# =============================================================================

#' @keywords internal
#' @noRd
# One-sample Hodges-Lehmann estimator (median of Walsh averages)
# Uses i <= j to match wilcox.test pseudomedian
hl1_est <- function(x) {
  n <- length(x)
  # Create all pairwise averages including diagonal (i <= j)
  # This gives n(n+1)/2 Walsh averages
  pairs <- outer(x, x, "+") / 2
  median(pairs[upper.tri(pairs, diag = TRUE)])
}

#' @keywords internal
#' @noRd
# Two-sample Hodges-Lehmann estimator (median of pairwise differences)
hl2_est <- function(x, y) {
  # All m*n pairwise differences: x_i - y_j
  # This matches the sign convention of wilcox.test, t.test, and perm_t_test
  diffs <- outer(x, y, "-")
  median(as.vector(diffs))
}

#' @keywords internal
#' @noRd
# Scale estimator S1: pooled within-sample absolute differences
scale_S1 <- function(x, y = NULL) {
  if (is.null(y)) {
    # One-sample case
    n <- length(x)
    if (n < 2) return(NA_real_)
    abs_diff <- abs(outer(x, x, "-"))
    return(median(abs_diff[lower.tri(abs_diff)]))
  }

  # Two-sample case: pool absolute differences from both samples
  m <- length(x)
  n <- length(y)

  abs_diff_x <- if (m >= 2) {
    mat <- abs(outer(x, x, "-"))
    mat[lower.tri(mat)]
  } else {
    numeric(0)
  }


  abs_diff_y <- if (n >= 2) {
    mat <- abs(outer(y, y, "-"))
    mat[lower.tri(mat)]
  } else {
    numeric(0)
  }

  pooled <- c(abs_diff_x, abs_diff_y)
  if (length(pooled) == 0) return(NA_real_)
  median(pooled)
}

#' @keywords internal
#' @noRd
# Scale estimator S2: median-corrected joint sample
scale_S2 <- function(x, y = NULL) {
  if (is.null(y)) {
    # One-sample: center by median
    x_centered <- x - median(x)
    z <- x_centered
  } else {
    # Two-sample: center each group by its median, then combine
    x_centered <- x - median(x)
    y_centered <- y - median(y)
    z <- c(x_centered, y_centered)
  }

  n <- length(z)
  if (n < 2) return(NA_real_)

  abs_diff_z <- abs(outer(z, z, "-"))
  median(abs_diff_z[lower.tri(abs_diff_z)])
}

#' @keywords internal
#' @noRd
# Generate sign-flipping permutations for one-sample/paired test
hl_perm_signs <- function(n, R) {
  max_perms <- 2^n

  if (is.null(R) || R >= max_perms) {
    # Exact: all sign combinations
    signs <- as.matrix(expand.grid(rep(list(c(-1, 1)), n)))
    list(signs = signs, exact = TRUE, R.used = max_perms, max_perms = max_perms)
  } else {
    # Randomization: sample with replacement
    signs <- matrix(sample(c(-1, 1), size = n * R, replace = TRUE),
                    nrow = R, ncol = n)
    list(signs = signs, exact = FALSE, R.used = R, max_perms = max_perms)
  }
}

#' @keywords internal
#' @noRd
# Generate group reassignment permutations for two-sample test
hl_perm_groups <- function(m, n, R) {
  n_total <- m + n
  max_perms <- choose(n_total, m)

  if (is.null(R) || R >= max_perms) {
    # Exact: all combinations
    idx <- utils::combn(n_total, m)
    list(idx_x = t(idx), exact = TRUE, R.used = max_perms, max_perms = max_perms)
  } else {
    # Randomization: sample with replacement
    idx_x <- t(replicate(R, sample(n_total, m)))
    list(idx_x = idx_x, exact = FALSE, R.used = R, max_perms = max_perms)
  }
}

#' @keywords internal
#' @noRd
# Compute permutation p-value
hl_perm_pval <- function(b, R, p_method) {
  if (p_method == "exact") {
    b / R
  } else {
    # plusone
    (b + 1) / (R + 1)
  }
}

#' @keywords internal
#' @noRd
# Kernel density estimate at zero for asymptotic variance
kde_at_zero <- function(diffs) {
  if (length(diffs) < 2) return(NA_real_)
  dens <- stats::density(diffs)
  # Find density value at x closest to 0
  idx <- which.min(abs(dens$x))
  dens$y[idx]
}

# =============================================================================
# Main implementation
# =============================================================================

#' @rdname hodges_lehmann
#' @method hodges_lehmann default
#' @export
hodges_lehmann.default <- function(x,
                                   y = NULL,
                                   alternative = c("two.sided",
                                                   "less",
                                                   "greater",
                                                   "equivalence",
                                                   "minimal.effect"),
                                   mu = 0,
                                   paired = FALSE,
                                   alpha = 0.05,
                                   R = NULL,
                                   scale = c("S2", "S1"),
                                   p_method = NULL,
                                   keep_perm = TRUE,
                                   ...) {

  alternative <- match.arg(alternative)
  scale <- match.arg(scale)

  # Input validation
  if (!is.numeric(x)) {
    stop("'x' must be numeric")
  }

  if (!is.null(y) && !is.numeric(y)) {
    stop("'y' must be numeric")
  }

  if (!missing(alpha) && (length(alpha) != 1 || !is.finite(alpha) ||
                          alpha < 0 || alpha > 1)) {
    stop("'alpha' must be a single number between 0 and 1")
  }

  if (!is.null(R) && (length(R) != 1 || !is.finite(R) || R < 1)) {
    stop("'R' must be NULL (for asymptotic) or a positive integer")
  }
  if (!is.null(R)) R <- as.integer(R)

  # Handle mu for equivalence/minimal.effect
  if (alternative %in% c("equivalence", "minimal.effect")) {
    if (length(mu) == 1) {
      if (mu == 0) {
        stop("For equivalence or minimal.effect testing, 'mu' must specify bounds (e.g., mu = c(-0.5, 0.5))")
      }
      mu <- c(-abs(mu), abs(mu))
    }
    if (length(mu) != 2) {
      stop("For equivalence or minimal.effect testing, 'mu' must be a vector of length 2")
    }
    mu <- sort(mu)
    low_bound <- mu[1]
    high_bound <- mu[2]
  } else {
    if (length(mu) != 1 || !is.finite(mu)) {
      stop("'mu' must be a single finite number for this alternative")
    }
  }

  # Data name extraction
  if (!is.null(y)) {
    dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
  } else {
    dname <- deparse(substitute(x))
  }

  # Check lengths for paired test before handling missing values
  if (paired && !is.null(y)) {
    if (length(x) != length(y)) {
      stop("'x' and 'y' must have the same length for paired test")
    }
  }

  # Handle missing values
  if (!is.null(y)) {
    if (paired) {
      ok <- complete.cases(x, y)
      x <- x[ok]
      y <- y[ok]
    } else {
      x <- x[!is.na(x)]
      y <- y[!is.na(y)]
    }
  } else {
    x <- x[!is.na(x)]
  }

  # Convert paired to one-sample problem
  if (paired && !is.null(y)) {
    if (length(x) != length(y)) {
      # This shouldn't happen after NA handling, but keep for safety
      stop("'x' and 'y' must have the same length for paired test after removing NAs")
    }
    x <- x - y
    y <- NULL
  }

  # Check minimum sample sizes
  if (is.null(y)) {
    # One-sample/paired
    n <- length(x)
    if (n < 2) {
      stop("need at least 2 observations")
    }
    max_perms <- 2^n
  } else {
    # Two-sample
    m <- length(x)
    n <- length(y)
    if (m < 2 || n < 2) {
      stop("need at least 2 observations in each group")
    }
    max_perms <- choose(m + n, m)
  }

  # Determine test method
  if (is.null(R)) {
    test_method <- "asymptotic"
  } else if (R >= max_perms) {
    test_method <- "exact"
  } else {
    test_method <- "randomization"
  }

  # Set p_method default
  if (is.null(p_method)) {
    p_method <- if (test_method == "exact") "exact" else "plusone"
  } else {
    p_method <- match.arg(p_method, c("exact", "plusone"))
  }

  # Compute confidence level
  if (alternative %in% c("equivalence", "minimal.effect")) {
    ci_level <- 1 - 2 * alpha
  } else {
    ci_level <- 1 - alpha
  }

  # ==========================================================================
  # One-sample / Paired test
  # ==========================================================================
  if (is.null(y)) {

    # Compute observed estimate
    estimate <- hl1_est(x)
    names(estimate) <- if (paired) "pseudomedian of differences" else "pseudomedian"

    if (test_method == "asymptotic") {
      # ----- Asymptotic test -----
      # Estimate variance using kernel density on pairwise differences
      # For one-sample: differences within sample
      diff_mat <- outer(x, x, "-")
      diffs <- diff_mat[lower.tri(diff_mat)]

      h0 <- kde_at_zero(diffs)
      if (is.na(h0) || h0 <= 0) {
        warning("Kernel density estimation failed; using permutation test instead")
        R <- 1999
        test_method <- "randomization"
      } else {
        # Asymptotic SE: from Fried & Dehling (2011)
        # SE = 1 / (sqrt(12 * n) * h(0))
        se <- 1 / (sqrt(12 * n) * h0)

        if (alternative %in% c("equivalence", "minimal.effect")) {
          z_low <- (estimate - low_bound) / se
          z_high <- (estimate - high_bound) / se
          tstat <- c(z_low, z_high)
        } else {
          tstat <- (estimate - mu) / se
        }

        # P-values
        if (alternative == "two.sided") {
          pval <- 2 * stats::pnorm(-abs(tstat))
        } else if (alternative == "less") {
          pval <- stats::pnorm(tstat)
        } else if (alternative == "greater") {
          pval <- stats::pnorm(tstat, lower.tail = FALSE)
        } else if (alternative == "equivalence") {
          p_low <- stats::pnorm(z_low, lower.tail = FALSE)
          p_high <- stats::pnorm(z_high)
          pval <- max(p_low, p_high)
        } else if (alternative == "minimal.effect") {
          p_low <- stats::pnorm(z_low)
          p_high <- stats::pnorm(z_high, lower.tail = FALSE)
          pval <- min(p_low, p_high)
        }

        # Confidence interval
        z_crit <- stats::qnorm(1 - (1 - ci_level) / 2)
        if (alternative == "two.sided") {
          cint <- estimate + c(-1, 1) * z_crit * se
        } else if (alternative == "less") {
          cint <- c(-Inf, estimate + z_crit * se)
        } else if (alternative == "greater") {
          cint <- c(estimate - z_crit * se, Inf)
        } else {
          # equivalence or minimal.effect
          z_crit_tost <- stats::qnorm(1 - alpha)
          cint <- estimate + c(-1, 1) * z_crit_tost * se
        }

        # For reporting: use the relevant test statistic
        if (alternative == "equivalence") {
          tstat_report <- if (p_low >= p_high) z_low else z_high
        } else if (alternative == "minimal.effect") {
          tstat_report <- if (p_low <= p_high) z_low else z_high
        } else {
          tstat_report <- tstat
        }

        R.used <- NA
        TSTAT <- NULL
        EFF <- NULL
      }
    }

    if (test_method %in% c("exact", "randomization")) {
      # ----- Permutation test -----
      # Center data at the observed estimate for permutation
      x_centered <- x - estimate

      # Generate permutations
      perm_result <- hl_perm_signs(n, R)
      signs_matrix <- perm_result$signs
      R.used <- perm_result$R.used
      exact_perm <- perm_result$exact

      if (exact_perm) {
        message("Computing all ", R.used, " exact permutations.")
      }

      # Compute scale estimate
      scale_est <- if (scale == "S1") scale_S1(x) else scale_S2(x)
      if (is.na(scale_est) || scale_est <= 0) {
        warning("Scale estimate is zero or NA; test may be unreliable")
        scale_est <- 1
      }

      # Compute permutation statistics
      TSTAT <- numeric(R.used)
      EFF <- numeric(R.used)

      for (i in 1:R.used) {
        x_perm <- x_centered * signs_matrix[i, ]
        hl_perm <- hl1_est(x_perm)
        TSTAT[i] <- hl_perm / scale_est
        EFF[i] <- hl_perm + estimate
      }

      # Observed test statistic
      obs_stat <- (estimate - mu) / scale_est

      # P-values
      if (alternative == "two.sided") {
        b <- sum(abs(TSTAT) >= abs(obs_stat))
        pval <- hl_perm_pval(b, R.used, p_method)
        cint <- stats::quantile(EFF, c((1 - ci_level) / 2, 1 - (1 - ci_level) / 2),
                                names = FALSE)
      } else if (alternative == "less") {
        b <- sum(TSTAT <= obs_stat)
        pval <- hl_perm_pval(b, R.used, p_method)
        cint <- c(-Inf, stats::quantile(EFF, ci_level, names = FALSE))
      } else if (alternative == "greater") {
        b <- sum(TSTAT >= obs_stat)
        pval <- hl_perm_pval(b, R.used, p_method)
        cint <- c(stats::quantile(EFF, 1 - ci_level, names = FALSE), Inf)
      } else if (alternative == "equivalence") {
        obs_stat_low <- (estimate - low_bound) / scale_est
        obs_stat_high <- (estimate - high_bound) / scale_est
        b_low <- sum(TSTAT >= obs_stat_low)
        b_high <- sum(TSTAT <= obs_stat_high)
        p_low <- hl_perm_pval(b_low, R.used, p_method)
        p_high <- hl_perm_pval(b_high, R.used, p_method)
        pval <- max(p_low, p_high)
        cint <- stats::quantile(EFF, c(alpha, 1 - alpha), names = FALSE)
        obs_stat <- if (p_low >= p_high) obs_stat_low else obs_stat_high
      } else if (alternative == "minimal.effect") {
        obs_stat_low <- (estimate - low_bound) / scale_est
        obs_stat_high <- (estimate - high_bound) / scale_est
        b_low <- sum(TSTAT <= obs_stat_low)
        b_high <- sum(TSTAT >= obs_stat_high)
        p_low <- hl_perm_pval(b_low, R.used, p_method)
        p_high <- hl_perm_pval(b_high, R.used, p_method)
        pval <- min(p_low, p_high)
        cint <- stats::quantile(EFF, c(alpha, 1 - alpha), names = FALSE)
        obs_stat <- if (p_low <= p_high) obs_stat_low else obs_stat_high
      }

      tstat_report <- obs_stat
    }

  } else {
    # ==========================================================================
    # Two-sample test
    # ==========================================================================
    m <- length(x)
    n <- length(y)

    # Compute observed estimate
    estimate <- hl2_est(x, y)
    names(estimate) <- "difference in location"

    if (test_method == "asymptotic") {
      # ----- Asymptotic test -----
      # Estimate h(0) using within-sample pairwise differences
      diff_x <- if (m >= 2) {
        mat <- outer(x, x, "-")
        mat[lower.tri(mat)]
      } else {
        numeric(0)
      }

      diff_y <- if (n >= 2) {
        mat <- outer(y, y, "-")
        mat[lower.tri(mat)]
      } else {
        numeric(0)
      }

      diffs <- c(diff_x, diff_y)

      h0 <- kde_at_zero(diffs)
      if (is.na(h0) || h0 <= 0) {
        warning("Kernel density estimation failed; using permutation test instead")
        R <- 1999
        test_method <- "randomization"
      } else {
        # Asymptotic SE from Fried & Dehling (2011)
        N <- m + n
        lambda <- m / N
        se <- 1 / (sqrt(12 * lambda * (1 - lambda) * N) * h0)

        if (alternative %in% c("equivalence", "minimal.effect")) {
          z_low <- (estimate - low_bound) / se
          z_high <- (estimate - high_bound) / se
          tstat <- c(z_low, z_high)
        } else {
          tstat <- (estimate - mu) / se
        }

        # P-values
        if (alternative == "two.sided") {
          pval <- 2 * stats::pnorm(-abs(tstat))
        } else if (alternative == "less") {
          pval <- stats::pnorm(tstat)
        } else if (alternative == "greater") {
          pval <- stats::pnorm(tstat, lower.tail = FALSE)
        } else if (alternative == "equivalence") {
          p_low <- stats::pnorm(z_low, lower.tail = FALSE)
          p_high <- stats::pnorm(z_high)
          pval <- max(p_low, p_high)
        } else if (alternative == "minimal.effect") {
          p_low <- stats::pnorm(z_low)
          p_high <- stats::pnorm(z_high, lower.tail = FALSE)
          pval <- min(p_low, p_high)
        }

        # Confidence interval
        z_crit <- stats::qnorm(1 - (1 - ci_level) / 2)
        if (alternative == "two.sided") {
          cint <- estimate + c(-1, 1) * z_crit * se
        } else if (alternative == "less") {
          cint <- c(-Inf, estimate + z_crit * se)
        } else if (alternative == "greater") {
          cint <- c(estimate - z_crit * se, Inf)
        } else {
          z_crit_tost <- stats::qnorm(1 - alpha)
          cint <- estimate + c(-1, 1) * z_crit_tost * se
        }

        # For reporting
        if (alternative == "equivalence") {
          tstat_report <- if (p_low >= p_high) z_low else z_high
        } else if (alternative == "minimal.effect") {
          tstat_report <- if (p_low <= p_high) z_low else z_high
        } else {
          tstat_report <- tstat
        }

        R.used <- NA
        TSTAT <- NULL
        EFF <- NULL
      }
    }

    if (test_method %in% c("exact", "randomization")) {
      # ----- Permutation test -----
      z <- c(x, y)

      # Generate permutations
      perm_result <- hl_perm_groups(m, n, R)
      idx_matrix <- perm_result$idx_x
      R.used <- perm_result$R.used
      exact_perm <- perm_result$exact

      if (exact_perm) {
        message("Computing all ", R.used, " exact permutations.")
      }

      # Compute scale estimate
      scale_est <- if (scale == "S1") scale_S1(x, y) else scale_S2(x, y)
      if (is.na(scale_est) || scale_est <= 0) {
        warning("Scale estimate is zero or NA; test may be unreliable")
        scale_est <- 1
      }

      # Compute permutation statistics
      TSTAT <- numeric(R.used)
      EFF <- numeric(R.used)

      for (i in 1:R.used) {
        idx_x <- idx_matrix[i, ]
        idx_y <- setdiff(1:(m + n), idx_x)

        x_perm <- z[idx_x]
        y_perm <- z[idx_y]

        hl_perm <- hl2_est(x_perm, y_perm)
        # Scale under permutation
        scale_perm <- if (scale == "S1") scale_S1(x_perm, y_perm) else scale_S2(x_perm, y_perm)
        if (is.na(scale_perm) || scale_perm <= 0) scale_perm <- scale_est

        TSTAT[i] <- hl_perm / scale_perm
        EFF[i] <- hl_perm
      }

      # Observed test statistic
      obs_stat <- (estimate - mu) / scale_est

      # P-values
      if (alternative == "two.sided") {
        b <- sum(abs(TSTAT) >= abs(obs_stat))
        pval <- hl_perm_pval(b, R.used, p_method)
        cint <- stats::quantile(EFF, c((1 - ci_level) / 2, 1 - (1 - ci_level) / 2),
                                names = FALSE)
      } else if (alternative == "less") {
        b <- sum(TSTAT <= obs_stat)
        pval <- hl_perm_pval(b, R.used, p_method)
        cint <- c(-Inf, stats::quantile(EFF, ci_level, names = FALSE))
      } else if (alternative == "greater") {
        b <- sum(TSTAT >= obs_stat)
        pval <- hl_perm_pval(b, R.used, p_method)
        cint <- c(stats::quantile(EFF, 1 - ci_level, names = FALSE), Inf)
      } else if (alternative == "equivalence") {
        obs_stat_low <- (estimate - low_bound) / scale_est
        obs_stat_high <- (estimate - high_bound) / scale_est
        b_low <- sum(TSTAT >= obs_stat_low)
        b_high <- sum(TSTAT <= obs_stat_high)
        p_low <- hl_perm_pval(b_low, R.used, p_method)
        p_high <- hl_perm_pval(b_high, R.used, p_method)
        pval <- max(p_low, p_high)
        cint <- stats::quantile(EFF, c(alpha, 1 - alpha), names = FALSE)
        obs_stat <- if (p_low >= p_high) obs_stat_low else obs_stat_high
      } else if (alternative == "minimal.effect") {
        obs_stat_low <- (estimate - low_bound) / scale_est
        obs_stat_high <- (estimate - high_bound) / scale_est
        b_low <- sum(TSTAT <= obs_stat_low)
        b_high <- sum(TSTAT >= obs_stat_high)
        p_low <- hl_perm_pval(b_low, R.used, p_method)
        p_high <- hl_perm_pval(b_high, R.used, p_method)
        pval <- min(p_low, p_high)
        cint <- stats::quantile(EFF, c(alpha, 1 - alpha), names = FALSE)
        obs_stat <- if (p_low <= p_high) obs_stat_low else obs_stat_high
      }

      tstat_report <- obs_stat
    }
  }

  # ==========================================================================
  # Build output
  # ==========================================================================

  # Null value
  if (alternative %in% c("equivalence", "minimal.effect")) {
    null.value <- mu
    names(null.value) <- rep("location", 2)
  } else {
    null.value <- mu
    names(null.value) <- "location"
  }

  # Confidence interval attributes
  attr(cint, "conf.level") <- ci_level

  # Method string
  if (test_method == "asymptotic") {
    method_prefix <- "Asymptotic"
  } else if (test_method == "exact") {
    method_prefix <- "Exact Permutation"
  } else {
    method_prefix <- "Randomization"
  }

  if (is.null(y)) {
    if (paired) {
      method <- paste(method_prefix, "Hodges-Lehmann Paired Test")
    } else {
      method <- paste(method_prefix, "Hodges-Lehmann One Sample Test")
    }
  } else {
    method <- paste(method_prefix, "Hodges-Lehmann Two Sample Test")
  }

  names(tstat_report) <- "Z"

  rval <- list(
    statistic = tstat_report,
    p.value = pval,
    conf.int = cint,
    estimate = estimate,
    null.value = null.value,
    alternative = alternative,
    method = method,
    data.name = dname,
    call = match.call(),
    R = R,
    R.used = R.used
  )

  if (keep_perm && !is.null(TSTAT)) {
    rval$perm.stat <- TSTAT
    rval$perm.eff <- EFF
  }

  class(rval) <- "htest"

  return(rval)
}

#' @rdname hodges_lehmann
#' @method hodges_lehmann formula
#' @export
hodges_lehmann.formula <- function(formula, data, subset, na.action, ...) {

  if (missing(formula) ||
      (length(formula) != 3L) ||
      (length(attr(stats::terms(formula[-2L]), "term.labels")) != 1L)) {
    stop("'formula' missing or incorrect")
  }

  # Check for paired argument in ... and warn user
  dots <- list(...)
  if ("paired" %in% names(dots)) {
    if (isTRUE(dots$paired)) {
      message("Using 'paired = TRUE' with the formula interface is not recommended. ",
              "Please ensure your data is sorted appropriately to make the correct paired comparison.")
    }
  }

  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval(m$data, parent.frame()))) {
    m$data <- as.data.frame(data)
  }

  m[[1L]] <- quote(stats::model.frame)
  m$... <- NULL
  mf <- eval(m, parent.frame())
  DNAME <- paste(names(mf), collapse = " by ")
  names(mf) <- NULL
  response <- attr(attr(mf, "terms"), "response")
  g <- factor(mf[[-response]])

  if (nlevels(g) != 2L) {
    stop("grouping factor must have exactly 2 levels")
  }

  DATA <- stats::setNames(split(mf[[response]], g), c("x", "y"))
  y <- do.call("hodges_lehmann", c(DATA, list(...)))
  y$data.name <- DNAME
  y
}
