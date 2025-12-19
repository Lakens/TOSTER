#' @title Permutation t-test
#' @description
#' `r lifecycle::badge('maturing')`
#'
#' Performs t-tests with permutation-based p-values and confidence intervals. This function supports
#' standard hypothesis testing alternatives as well as equivalence and minimal effect testing,
#' with optional trimmed means (Yuen's approach) for robust inference.
#'
#' @section Purpose:
#' Use this function when:
#'   * You need more robust inference than provided by standard t-tests

#'   * You want to perform equivalence or minimal effect testing with permutation methods
#'   * Sample sizes are small or standard parametric approaches may be kess reliable
#'   * You prefer the standard `htest` output format for compatibility with other R functions
#'   * You want an alternative that tests means, mean difference, or difference in means that can also use trimming
#'
#' @inheritParams simple_htest
#' @param mu a number or vector specifying the null hypothesis value(s):
#'     * For standard alternatives: for standard alternatives (two.sided, less, greater): a single value (default: 0
#'     * For equivalence/minimal.effect: either a single value (symmetric bounds will be created) or a vector of two values representing the lower and upper bounds
#'
#' @param paired a logical indicating whether you want a paired t-test.
#' @param var.equal a logical variable indicating whether to treat the two variances as being equal.
#'     If `TRUE` then the pooled variance is used to estimate the variance otherwise the Welch
#'     (or Satterthwaite) approximation to the degrees of freedom is used. Default is `FALSE`.
#' @param alpha significance level (default = 0.05).
#' @param trim the fraction (0 to 0.5) of observations to be trimmed from each end before
#'     computing the mean and winsorized variance. Default is 0 (no trimming).
#'     When trim > 0, the function performs Yuen's trimmed t-test.
#' @param R the number of permutations. Default is 1999.
#' @param symmetric a logical variable indicating whether to assume symmetry in the two-sided test.
#'     If `TRUE` (default) then the symmetric permutation p-value is computed, otherwise the
#'     equal-tail permutation p-value is computed. Only relevant for `alternative = "two.sided"`.
#' @param keep_perm logical. If `TRUE` (default), the permutation distribution of the test
#'     statistic and effect sizes are stored in the output. Set to `FALSE` for large datasets
#'     to save memory.
#' @param ... further arguments (currently ignored).
#'
#' @details
#' This function performs permutation-based t-tests, providing more robust inference than standard
#' parametric t-tests. It supports one-sample, two-sample (independent), and paired designs,
#' as well as five different alternative hypotheses.
#'
#' The permutation procedure follows these steps:
#'   * Calculate the test statistic from the original data
#'   * Generate R permutation samples by randomly reassigning observations to groups
#'   * Calculate the test statistic for each permutation sample
#'   * Compute the p-value by comparing the original test statistic to the permutation distribution
#'   * Calculate confidence intervals using the percentile method
#'
#' For one-sample and paired tests, permutation is performed by randomly flipping the signs of
#' the (centered) observations or difference scores.
#'
#' For two-sample tests, permutation is performed by randomly reassigning observations to the
#' two groups. The studentized permutation approach of Janssen (1997) and Chung & Romano (2013)
#' is used, which is valid even under heteroscedasticity.
#'
#' When `trim > 0`, the function uses Yuen's trimmed t-test approach:
#'   * Trimmed means are computed by removing the fraction `trim` of observations from each tail
#'   * Winsorized variances are used in place of standard variances
#'   * This provides robustness against outliers and heavy-tailed distributions
#'
#' For different alternatives, the p-values are calculated as follows:
#'   * "two.sided": Either symmetric (proportion of |T\*| >= |T|) or equal-tail
#'       (2 * min of one-sided p-values) depending on the `symmetric` argument
#'   * "less": Proportion of T\* <= T
#'   * "greater": Proportion of T\* >= T
#'   * "equivalence": Maximum of two one-sided p-values (for lower and upper bounds)
#'   * "minimal.effect": Minimum of two one-sided p-values (for lower and upper bounds)
#'
#' For two-sample tests, the test is of \eqn{\bar x - \bar y} (mean of x minus mean of y).
#' For paired samples, the test is of the difference scores (z),
#' wherein \eqn{z = x - y}, and the test is of \eqn{\bar z} (mean of the difference scores).
#' For one-sample tests, the test is of \eqn{\bar x} (mean of x).
#'
#' If the number of possible permutations is less than R, all exact permutations are computed
#' and a message is printed to the console.
#'
#' @return A list with class `"htest"` containing the following components:
#'
#'   - "p.value": the permutation p-value for the test.
#'   - "stderr": the standard error of the mean (difference).
#'   - "conf.int": a permutation percentile confidence interval appropriate to the
#'       specified alternative hypothesis.
#'   - "estimate": the estimated mean or difference in means (or trimmed means if trim > 0).
#'   - "null.value": the specified hypothesized value(s) of the mean or mean difference.
#'   - "alternative": a character string describing the alternative hypothesis.
#'   - "method": a character string indicating what type of permutation t-test was performed.
#'   - "data.name": a character string giving the name(s) of the data.
#'   - "call": the matched call.
#'   - "R": the requested number of permutations.
#'   - "R.used": the actual number of permutations used (may differ if exact permutations computed).
#'   - "perm.stat": (if keep_perm = TRUE) the permutation distribution of the test statistic.
#'   - "perm.eff": (if keep_perm = TRUE) the permutation distribution of the effect (mean differences).
#'
#' @examples
#'
#' # Example 1: Basic two-sample test with formula notation
#' data(sleep)
#' result <- perm_t_test(extra ~ group, data = sleep)
#' result
#'
#' # Example 2: One-sample permutation t-test
#' set.seed(123)
#' x <- rnorm(20, mean = 0.5, sd = 1)
#' perm_t_test(x, mu = 0, R = 199)
#'
#' # Example 3: Paired samples test
#' before <- c(5.1, 4.8, 6.2, 5.7, 6.0, 5.5, 4.9, 5.8)
#' after <- c(5.6, 5.2, 6.7, 6.1, 6.5, 5.8, 5.3, 6.2)
#' perm_t_test(x = before, y = after,
#'             paired = TRUE,
#'             alternative = "less",
#'             R = 199)
#'
#' # Example 4: Equivalence testing with permutation t-test
#' data(mtcars)
#' perm_t_test(mpg ~ am, data = mtcars,
#'             alternative = "equivalence",
#'             mu = c(-3, 3),
#'             R = 199)
#'
#' # Example 5: Trimmed permutation t-test (Yuen's approach)
#' # Useful when data may contain outliers
#' set.seed(42)
#' x <- c(rnorm(18), 10, 15)  # Two outliers
#' y <- rnorm(20)
#' perm_t_test(x, y, trim = 0.1, R = 999)
#'
#' # Example 6: Minimal effect testing
#' perm_t_test(mpg ~ am, data = mtcars,
#'             alternative = "minimal.effect",
#'             mu = c(-3, 3),
#'             R = 199)
#'
#' @references
#' Efron, B., & Tibshirani, R. J. (1993). An Introduction to the Bootstrap. Chapman and Hall/CRC.
#'
#' Janssen, A. (1997). Studentized permutation tests for non-i.i.d. hypotheses and the
#' generalized Behrens-Fisher problem. Statistics & Probability Letters, 36, 9-21.
#'
#' Chung, E., & Romano, J. P. (2013). Exact and asymptotically robust permutation tests.
#' The Annals of Statistics, 41(2), 484-507.
#'
#' Yuen, K. K. (1974). The two-sample trimmed t for unequal population variances.
#' Biometrika, 61(1), 165-170.
#'
#' @family Robust tests
#' @name perm_t_test
#' @export perm_t_test

#' @keywords internal
#' @noRd
# Helper function to compute trimmed mean
trimmed_mean <- function(x, trim = 0) {
  mean(x, trim = trim, na.rm = TRUE)
}

#' @keywords internal
#' @noRd
# Helper function to compute winsorized variance
winsorized_var <- function(x, trim = 0) {

  n <- length(x)
  if (trim == 0) {
    return(var(x))
  }

  g <- floor(trim * n)
  if (g == 0) {
    return(var(x))
  }

  x_sorted <- sort(x)
  # Winsorize: replace trimmed values with the adjacent non-trimmed values
  x_wins <- x_sorted
  x_wins[1:g] <- x_sorted[g + 1]
  x_wins[(n - g + 1):n] <- x_sorted[n - g]

  # Return winsorized variance
  var(x_wins)
}

#' @keywords internal
#' @noRd
# Helper function to compute effective sample size after trimming
effective_n <- function(n, trim) {
  n - 2 * floor(trim * n)
}

#' @keywords internal
#' @noRd
# Helper function to compute degrees of freedom for trimmed t-test (Yuen-Welch)
yuen_welch_df <- function(nx, ny, vx_wins, vy_wins, trim) {
  hx <- effective_n(nx, trim)
  hy <- effective_n(ny, trim)

  # Winsorized standard errors
  dx <- (nx - 1) * vx_wins / (hx * (hx - 1))
  dy <- (ny - 1) * vy_wins / (hy * (hy - 1))

  # Welch-Satterthwaite degrees of freedom
  df <- (dx + dy)^2 / (dx^2 / (hx - 1) + dy^2 / (hy - 1))

  return(df)
}

#' @keywords internal
#' @noRd
# Helper function to compute pooled winsorized df
pooled_wins_df <- function(nx, ny, trim) {
  hx <- effective_n(nx, trim)
  hy <- effective_n(ny, trim)
  hx + hy - 2
}

#' @keywords internal
#' @noRd
# Helper function to generate permutation samples for one-sample/paired case
perm_signs <- function(n, R) {
  # Total possible permutations for sign flipping
  max_perms <- 2^n

  if (max_perms <= R) {
    # Generate all possible sign combinations
    message("Computing all ", max_perms, " exact permutations.")
    signs <- as.matrix(expand.grid(rep(list(c(-1, 1)), n)))
    return(list(signs = signs, exact = TRUE, R.used = max_perms))
  } else {
    # Monte Carlo sampling
    signs <- matrix(sample(c(-1, 1), size = n * R, replace = TRUE), nrow = R, ncol = n)
    return(list(signs = signs, exact = FALSE, R.used = R))
  }
}

#' @keywords internal
#' @noRd
# Helper function to generate permutation samples for two-sample case
perm_groups <- function(nx, ny, R) {
  n_total <- nx + ny
  max_perms <- choose(n_total, nx)

  if (max_perms <= R) {
    # Generate all possible combinations
    message("Computing all ", max_perms, " exact permutations.")
    idx <- utils::combn(n_total, nx)
    # Convert to matrix where each row is a permutation
    perms <- t(idx)
    return(list(idx_x = perms, exact = TRUE, R.used = max_perms))
  } else {
    # Monte Carlo sampling - ensure unique permutations
    perms <- matrix(NA, nrow = R, ncol = nx)
    for (i in 1:R) {
      perms[i, ] <- sort(sample(n_total, nx))
    }
    # Remove duplicates
    perms <- unique(perms)
    iter <- 0
    max_iter <- 10
    while (nrow(perms) < R && iter < max_iter) {
      needed <- R - nrow(perms)
      new_perms <- matrix(NA, nrow = needed * 2, ncol = nx)
      for (i in 1:(needed * 2)) {
        new_perms[i, ] <- sort(sample(n_total, nx))
      }
      perms <- unique(rbind(perms, new_perms))
      iter <- iter + 1
    }
    if (nrow(perms) > R) {
      perms <- perms[1:R, , drop = FALSE]
    }
    return(list(idx_x = perms, exact = FALSE, R.used = nrow(perms)))
  }
}

perm_t_test <- function(x, ...) {
  UseMethod("perm_t_test")
}

#' @rdname perm_t_test
#' @method perm_t_test default
#' @export

perm_t_test.default <- function(x,
                                y = NULL,
                                alternative = c("two.sided",
                                                "less",
                                                "greater",
                                                "equivalence",
                                                "minimal.effect"),
                                mu = 0,
                                paired = FALSE,
                                var.equal = FALSE,
                                alpha = 0.05,
                                trim = 0,
                                R = 1999,
                                symmetric = TRUE,
                                keep_perm = TRUE,
                                ...) {

  alternative <- match.arg(alternative)

  # Input validation
  if (!missing(alpha) && (length(alpha) != 1 || !is.finite(alpha) ||
                          alpha < 0 || alpha > 1)) {
    stop("'alpha' must be a single number between 0 and 1")
  }

  if (!missing(trim) && (length(trim) != 1 || !is.finite(trim) ||
                         trim < 0 || trim >= 0.5)) {
    stop("'trim' must be a single number between 0 and 0.5 (exclusive)")
  }

  if (!missing(R) && (length(R) != 1 || !is.finite(R) || R < 1)) {
    stop("'R' must be a positive integer")
  }
  R <- as.integer(R)

  # Handle mu for equivalence/minimal.effect
  if (alternative %in% c("equivalence", "minimal.effect")) {
    if (length(mu) == 1) {
      if (mu == 0) {
        stop("For equivalence or minimal.effect testing, 'mu' must specify two bounds (e.g., mu = c(-0.5, 0.5))")
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
    if (paired) {
      stop("'y' is missing for paired test")
    }
  }

  # Convert paired to one-sample problem
  if (paired && !is.null(y)) {
    x <- x - y
    y <- NULL
  }

  nx <- length(x)

  # Check minimum sample size based on trimming
  # Need at least 2 observations remaining after trimming for variance calculation
  if (trim > 0) {
    g <- floor(trim * nx)
    effective <- nx - 2 * g
    if (effective < 2) {
      min_n <- ceiling(2 / (1 - 2 * trim)) + 1
      stop("Sample size too small for specified trimming proportion. ",
           "With trim = ", trim, ", need at least ", min_n, " observations, but only have ", nx, ".")
    }
  } else {
    if (nx < 2) {
      stop("not enough 'x' observations")
    }
  }

  # Compute statistics for one-sample/paired case
  if (is.null(y)) {
    # One sample or paired
    mx <- trimmed_mean(x, trim)
    vx <- winsorized_var(x, trim)
    hx <- effective_n(nx, trim)

    if (trim > 0) {
      # Yuen's approach: use winsorized SE
      stderr <- sqrt((nx - 1) * vx / (hx * (hx - 1)))
    } else {
      stderr <- sqrt(vx / nx)
    }

    df <- hx - 1

    if (stderr < 10 * .Machine$double.eps * abs(mx)) {
      stop("data are essentially constant")
    }

    if (alternative %in% c("equivalence", "minimal.effect")) {
      tstat_low <- (mx - low_bound) / stderr
      tstat_high <- (mx - high_bound) / stderr
      tstat <- c(tstat_low, tstat_high)
    } else {
      tstat <- (mx - mu) / stderr
    }

    method <- if (paired) {
      if (trim > 0) paste0(method_prefix," Paired Yuen t-test") else paste0(method_prefix," Paired t-test")
    } else {
      if (trim > 0) paste0(method_prefix," One Sample Yuen t-test") else paste0(method_prefix," One Sample t-test")
    }

    estimate <- setNames(mx, if (paired) {
      if (trim > 0) "trimmed mean of the differences" else "mean of the differences"
    } else {
      if (trim > 0) "trimmed mean of x" else "mean of x"
    })

    # Generate permutation distribution
    # For one-sample: flip signs of centered data
    x_centered <- x - mx

    perm_result <- perm_signs(nx, R)
    signs_matrix <- perm_result$signs
    R_used <- perm_result$R.used
    exact_perm <- perm_result$exact
    if (exact_perm) {
      method_prefix <- paste("Exact Permutation", ...)
    } else {
      method_prefix <- paste("Monte Carlo Permutation", ...)
    }
    # Compute permutation statistics
    TSTAT <- numeric(R_used)
    EFF <- numeric(R_used)

    for (i in 1:R_used) {
      x_perm <- abs(x_centered) * signs_matrix[i, ]
      mx_perm <- trimmed_mean(x_perm, trim)
      vx_perm <- winsorized_var(x_perm, trim)

      if (trim > 0) {
        stderr_perm <- sqrt((nx - 1) * vx_perm / (hx * (hx - 1)))
      } else {
        stderr_perm <- sqrt(vx_perm / nx)
      }

      if (stderr_perm > 10 * .Machine$double.eps) {
        TSTAT[i] <- mx_perm / stderr_perm
      } else {
        TSTAT[i] <- 0
      }
      EFF[i] <- mx_perm + mx
    }

  } else {
    # Two-sample case
    ny <- length(y)

    # Check minimum sample size for y
    if (trim > 0) {
      g <- floor(trim * ny)
      effective <- ny - 2 * g
      if (effective < 2) {
        min_n <- ceiling(2 / (1 - 2 * trim)) + 1
        stop("Sample size too small for specified trimming proportion. ",
             "With trim = ", trim, ", need at least ", min_n, " observations in each group, but y has only ", ny, ".")
      }
    } else {
      if (!var.equal && ny < 2) {
        stop("not enough 'y' observations")
      }
      if (var.equal && (nx + ny) < 3) {
        stop("not enough observations")
      }
    }

    mx <- trimmed_mean(x, trim)
    my <- trimmed_mean(y, trim)
    vx <- winsorized_var(x, trim)
    vy <- winsorized_var(y, trim)
    hx <- effective_n(nx, trim)
    hy <- effective_n(ny, trim)

    if (var.equal) {
      # Pooled variance approach
      if (trim > 0) {
        # Pooled winsorized variance
        df <- pooled_wins_df(nx, ny, trim)
        v_pooled <- ((hx - 1) * vx + (hy - 1) * vy) / df
        stderr <- sqrt(v_pooled * ((nx - 1) / (hx * (hx - 1)) + (ny - 1) / (hy * (hy - 1))))
      } else {
        df <- nx + ny - 2
        v_pooled <- ((nx - 1) * vx + (ny - 1) * vy) / df
        stderr <- sqrt(v_pooled * (1/nx + 1/ny))
      }
    } else {
      # Welch/Yuen-Welch approach
      if (trim > 0) {
        dx <- (nx - 1) * vx / (hx * (hx - 1))
        dy <- (ny - 1) * vy / (hy * (hy - 1))
        stderr <- sqrt(dx + dy)
        df <- yuen_welch_df(nx, ny, vx, vy, trim)
      } else {
        stderrx <- sqrt(vx / nx)
        stderry <- sqrt(vy / ny)
        stderr <- sqrt(stderrx^2 + stderry^2)
        df <- stderr^4 / (stderrx^4 / (nx - 1) + stderry^4 / (ny - 1))
      }
    }

    if (stderr < 10 * .Machine$double.eps * max(abs(mx), abs(my))) {
      stop("data are essentially constant")
    }

    diff_means <- mx - my

    if (alternative %in% c("equivalence", "minimal.effect")) {
      tstat_low <- (diff_means - low_bound) / stderr
      tstat_high <- (diff_means - high_bound) / stderr
      tstat <- c(tstat_low, tstat_high)
    } else {
      tstat <- (diff_means - mu) / stderr
    }

    method <- paste(method_prefix,
                    if (!var.equal) "Welch",
                    if (trim > 0) "Yuen",
                    "Two Sample t-test")
    method <- gsub("\\s+", " ", method)  # Clean up extra spaces

    if (trim > 0) {
      estimate <- c(mx, my)
      names(estimate) <- c("trimmed mean of x", "trimmed mean of y")
    } else {
      estimate <- c(mx, my)
      names(estimate) <- c("mean of x", "mean of y")
    }

    # Generate permutation distribution
    z <- c(x, y)

    perm_result <- perm_groups(nx, ny, R)
    idx_x_matrix <- perm_result$idx_x
    R_used <- perm_result$R.used
    exact_perm <- perm_result$exact

    # Compute permutation statistics
    # Following the studentized permutation approach:
    # Permute labels, compute test statistic, compare to observed
    TSTAT <- numeric(R_used)
    EFF <- numeric(R_used)

    for (i in 1:R_used) {
      idx_x <- idx_x_matrix[i, ]
      idx_y <- setdiff(1:(nx + ny), idx_x)

      x_perm <- z[idx_x]
      y_perm <- z[idx_y]

      mx_perm <- trimmed_mean(x_perm, trim)
      my_perm <- trimmed_mean(y_perm, trim)
      vx_perm <- winsorized_var(x_perm, trim)
      vy_perm <- winsorized_var(y_perm, trim)

      if (var.equal) {
        if (trim > 0) {
          df_perm <- pooled_wins_df(nx, ny, trim)
          v_pooled_perm <- ((hx - 1) * vx_perm + (hy - 1) * vy_perm) / df_perm
          stderr_perm <- sqrt(v_pooled_perm * ((nx - 1) / (hx * (hx - 1)) + (ny - 1) / (hy * (hy - 1))))
        } else {
          df_perm <- nx + ny - 2
          v_pooled_perm <- ((nx - 1) * vx_perm + (ny - 1) * vy_perm) / df_perm
          stderr_perm <- sqrt(v_pooled_perm * (1/nx + 1/ny))
        }
      } else {
        if (trim > 0) {
          dx_perm <- (nx - 1) * vx_perm / (hx * (hx - 1))
          dy_perm <- (ny - 1) * vy_perm / (hy * (hy - 1))
          stderr_perm <- sqrt(dx_perm + dy_perm)
        } else {
          stderrx_perm <- sqrt(vx_perm / nx)
          stderry_perm <- sqrt(vy_perm / ny)
          stderr_perm <- sqrt(stderrx_perm^2 + stderry_perm^2)
        }
      }

      # Compute test statistic under permutation
      # Under the null hypothesis, the group labels are exchangeable
      # so we compute the raw difference standardized by permutation SE
      if (stderr_perm > 10 * .Machine$double.eps) {
        TSTAT[i] <- (mx_perm - my_perm) / stderr_perm
      } else {
        TSTAT[i] <- 0
      }

      # Effect for CI: permutation difference centered around observed
      # This gives us the distribution of the effect under repeated sampling
      EFF[i] <- (mx_perm - my_perm) + diff_means
    }
  }

  # Compute confidence level based on alternative
  if (alternative %in% c("equivalence", "minimal.effect")) {
    conf.level <- 1 - alpha * 2
  } else {
    conf.level <- 1 - alpha
  }

  # Compute p-values and confidence intervals
  if (alternative == "less") {
    if (is.null(y)) {
      perm.pval <- mean(TSTAT <= tstat)
    } else {
      perm.pval <- mean(TSTAT <= tstat)
    }
    perm.pval <- max(perm.pval, 1 / R_used)
    perm.cint <- c(-Inf, quantile(EFF, conf.level, names = FALSE))

  } else if (alternative == "greater") {
    perm.pval <- mean(TSTAT >= tstat)
    perm.pval <- max(perm.pval, 1 / R_used)
    perm.cint <- c(quantile(EFF, 1 - conf.level, names = FALSE), Inf)

  } else if (alternative == "two.sided") {
    if (symmetric) {
      perm.pval <- mean(abs(TSTAT) >= abs(tstat))
    } else {
      perm.pval <- 2 * min(mean(TSTAT <= tstat), mean(TSTAT > tstat))
    }
    perm.pval <- max(perm.pval, 1 / R_used)
    perm.cint <- quantile(EFF, c(alpha / 2, 1 - alpha / 2), names = FALSE)

  } else if (alternative == "equivalence") {
    # Two one-sided tests
    p_low <- mean(TSTAT >= tstat[1])   # Test against lower bound
    p_high <- mean(TSTAT <= tstat[2])  # Test against upper bound
    perm.pval <- max(p_low, p_high)
    perm.pval <- max(perm.pval, 1 / R_used)
    perm.cint <- quantile(EFF, c(alpha, 1 - alpha), names = FALSE)

  } else if (alternative == "minimal.effect") {
    # Two one-sided tests (opposite direction)
    p_low <- mean(TSTAT <= tstat[1])   # Test against lower bound
    p_high <- mean(TSTAT >= tstat[2])  # Test against upper bound
    perm.pval <- min(p_low, p_high)
    perm.pval <- max(perm.pval, 1 / R_used)
    perm.cint <- quantile(EFF, c(alpha, 1 - alpha), names = FALSE)
  }

  # Set up null value
  if (alternative %in% c("equivalence", "minimal.effect")) {
    null.value <- mu
    names(null.value) <- rep("difference in means", 2)
  } else {
    null.value <- mu
    names(null.value) <- if (paired || !is.null(y)) "difference in means" else "mean"
  }

  attr(perm.cint, "conf.level") <- conf.level

  # Build output
  rval <- list(
    p.value = perm.pval,
    stderr = stderr,
    conf.int = perm.cint,
    estimate = estimate,
    null.value = null.value,
    alternative = alternative,
    method = method,
    data.name = dname,
    call = match.call(),
    R = R,
    R.used = R_used
  )

  if (keep_perm) {
    rval$perm.stat <- TSTAT
    rval$perm.eff <- EFF
  }

  class(rval) <- "htest"

  return(rval)
}

#' @rdname perm_t_test
#' @method perm_t_test formula
#' @export

perm_t_test.formula <- function(formula, data, subset, na.action, ...) {
  if (missing(formula) ||
      (length(formula) != 3L) ||
      (length(attr(terms(formula[-2L]), "term.labels")) != 1L)) {
    stop("'formula' missing or incorrect")
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

  DATA <- setNames(split(mf[[response]], g), c("x", "y"))
  y <- do.call("perm_t_test", c(DATA, list(...)))
  y$data.name <- DNAME
  y
}
