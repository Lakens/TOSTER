#' @title Bootstrapped t-test
#' @description
#' `r lifecycle::badge('stable')`
#'
#' Performs t-tests with bootstrapped p-values and confidence intervals, with optional
#' trimmed means (Yuen's approach) for robust inference. This function supports
#' standard hypothesis testing alternatives as well as equivalence and minimal effect testing,
#' all with the familiar `htest` output structure.
#'
#' @section Purpose:
#' Use this function when:
#'   * You need more robust inference than provided by standard t-tests
#'   * Your data don't meet the assumptions of normality or homogeneity
#'   * You want to perform equivalence or minimal effect testing with bootstrap methods
#'   * Sample sizes are small or standard parametric approaches may be unreliable
#'   * You prefer the standard `htest` output format for compatibility with other R functions
#'
#' @inheritParams simple_htest
#' @inheritParams boot_t_TOST
#' @param alternative the alternative hypothesis:
#'     * "two.sided": different from mu (default)
#'     * "less": less than mu
#'     * "greater": greater than mu
#'     * "equivalence": between specified bounds
#'     * "minimal.effect": outside specified bounds
#'
#' @param mu a number or vector specifying the null hypothesis value(s):
#'     * For standard alternatives: a single value (default = 0)
#'     * For equivalence/minimal.effect: two values representing the lower and upper bounds
#'
#' @param tr the fraction (0 to 0.5) of observations to be trimmed from
#'     each end before computing the mean and winsorized variance.
#'     Default is 0 (no trimming). When tr > 0, the function performs
#'     a bootstrapped Yuen's trimmed t-test.
#'
#' @details
#' This function performs bootstrapped t-tests, providing more robust inference than standard
#' parametric t-tests. It supports one-sample, two-sample (independent), and paired designs,
#' as well as five different alternative hypotheses.
#'
#' The bootstrap procedure follows these steps:
#'   * Calculate the test statistic from the original data
#'   * Generate R bootstrap samples by resampling with replacement
#'   * Calculate the test statistic for each bootstrap sample
#'   * Compute the p-value by comparing the original test statistic to the bootstrap distribution
#'   * Calculate confidence intervals using the specified bootstrap method
#'
#' ## Bootstrap Confidence Interval Methods
#'
#' Four bootstrap confidence interval methods are available via the `boot_ci` argument:
#'   - **Studentized bootstrap ("stud")**: Uses the bootstrap distribution of pivotal
#'     t-statistics to account for variability in standard error estimates. This is the
#'     default and usually provides the most accurate coverage.
#'   - **Basic bootstrap ("basic")**: Reflects the bootstrap distribution of estimates
#'     around the observed value.
#'   - **Percentile bootstrap ("perc")**: Uses percentiles of the bootstrap distribution directly.
#'   - **Bias-corrected and accelerated ("bca")**: Corrects for both bias and skewness in the
#'     bootstrap distribution using jackknife-based acceleration. Most accurate when the
#'     bootstrap distribution is skewed, but computationally more expensive.
#'
#' ## Bootstrap P-values
#'
#' The p-value is computed using the method that matches the selected `boot_ci`,
#' ensuring that p < alpha if and only if the corresponding confidence interval
#' excludes the null value (CI inversion principle). Previously, all bootstrap
#' CI methods used the studentized (pivot) p-value, which could produce p-values
#' inconsistent with non-studentized CIs.
#'
#' For different alternatives, the p-values are calculated as follows:
#'   * "two.sided": Two-tailed p-value from the bootstrap distribution
#'   * "less": One-sided p-value for the hypothesis that the true value is less than the null
#'   * "greater": One-sided p-value for the hypothesis that the true value is greater than the null
#'   * "equivalence": Maximum of two one-sided p-values (for lower and upper bounds)
#'   * "minimal.effect": Minimum of two one-sided p-values (for lower and upper bounds)
#'
#'
#' For two-sample tests, the test is of \eqn{\bar x - \bar y} (mean of x minus mean of y).
#' For paired samples, the test is of the difference scores (z),
#' wherein \eqn{z = x - y}, and the test is of \eqn{\bar z} (mean of the difference scores).
#' For one-sample tests, the test is of \eqn{\bar x} (mean of x).
#'
#' When `tr > 0`, the function uses Yuen's trimmed t-test approach: trimmed means
#' are computed by removing the fraction `tr` of observations from each tail,
#' and winsorized variances are used in place of standard variances. This provides
#' robustness against outliers and heavy-tailed distributions. The bootstrap
#' procedure recomputes trimmed means and winsorized standard errors for each
#' bootstrap replicate.
#'
#' Unlike the `t_TOST` function, this function returns a standard `htest` object for
#' compatibility with other R functions, while still providing the benefits of bootstrapping.
#'
#' For detailed information on calculation methods, see `vignette("robustTOST")`.
#'
#' @return A list with class `"htest"` containing the following components:
#'
#'   - "statistic": the observed t-statistic (note: the p-value is derived from the bootstrap
#'       distribution, not from this statistic and the degrees of freedom).
#'   - "parameter": the degrees of freedom for the t-statistic.
#'   - "p.value": the bootstrapped p-value for the test.
#'   - "stderr": the bootstrapped standard error.
#'   - "conf.int": a bootstrapped confidence interval for the mean appropriate to the specified alternative hypothesis.
#'   - "estimate": the estimated mean or difference in means.
#'   - "null.value": the specified hypothesized value(s) of the mean or mean difference.
#'   - "alternative": a character string describing the alternative hypothesis.
#'   - "method": a character string indicating what type of bootstrapped t-test was performed.
#'   - "boot": the bootstrap samples of the mean or mean difference.
#'   - "data.name": a character string giving the name(s) of the data.
#'   - "call": the matched call.
#'
#' @examples
#'
#' # Example 1: Basic two-sample test with formula notation
#' data(sleep)
#' result <- boot_t_test(extra ~ group, data = sleep)
#' result  # Standard htest output format
#'
#' # Example 2: One-sample bootstrapped t-test
#' set.seed(123)
#' x <- rnorm(20, mean = 0.5, sd = 1)
#' boot_t_test(x, mu = 0, R = 999) # Using fewer replicates for demonstration
#'
#' # Example 3: Paired samples test with percentile bootstrap CI
#' before <- c(5.1, 4.8, 6.2, 5.7, 6.0, 5.5, 4.9, 5.8)
#' after <- c(5.6, 5.2, 6.7, 6.1, 6.5, 5.8, 5.3, 6.2)
#' boot_t_test(x = before, y = after,
#'             paired = TRUE,
#'             alternative = "less",  # Testing if before < after
#'             boot_ci = "perc",
#'             R = 999)
#'
#' # Example 4: Equivalence testing with bootstrapped t-test
#' # Testing if the effect is within ±0.5 units
#' data(mtcars)
#' boot_t_test(mpg ~ am, data = mtcars,
#'             alternative = "equivalence",
#'             mu = c(-0.5, 0.5),
#'             boot_ci = "stud",
#'             R = 999)
#'
#' # Example 5: Minimal effect testing with bootstrapped t-test
#' # Testing if the effect is outside ±3 units
#' boot_t_test(mpg ~ am, data = mtcars,
#'             alternative = "minimal.effect",
#'             mu = c(-3, 3),
#'             R = 999)
#'
#' # Example 6: Bootstrapped Yuen's trimmed t-test (10% trimming)
#' boot_t_test(extra ~ group, data = sleep, tr = 0.1, R = 999)
#'
#' @references
#' Efron, B., & Tibshirani, R. J. (1994). An introduction to the bootstrap. CRC press.
#'
#' Yuen, K. K. (1974). The two-sample trimmed t for unequal population variances.
#' Biometrika, 61(1), 165-170.
#'
#' @family Robust tests
#' @name boot_t_test
#' @export boot_t_test

# TODO: add xname and yname arguments to allow user-specified group labels
boot_t_test <- function(x, ...){
  UseMethod("boot_t_test")
}

#' @rdname boot_t_test
#' @method boot_t_test default
#' @export

boot_t_test.default <- function(x,
                                y = NULL,
                                var.equal = FALSE,
                                paired = FALSE,
                                alternative = c("two.sided",
                                                "less",
                                                "greater",
                                                "equivalence",
                                                "minimal.effect"),
                                mu = 0,
                                alpha = 0.05,
                                tr = 0,
                                boot_ci = c("stud","basic","perc","bca"),
                                R = 1999, ...){
  alternative = match.arg(alternative)
  boot_ci = match.arg(boot_ci)

  if(!missing(alpha) && (length(alpha) != 1 || !is.finite(alpha) ||
                         alpha < 0 || alpha > 1)) {
    stop("'alpha' must be a single number between 0 and 1")
  }

  if (!missing(tr) && (length(tr) != 1 || !is.finite(tr) ||
                       tr < 0 || tr >= 0.5)) {
    stop("'tr' must be a single number between 0 and 0.5 (exclusive)")
  }


  if (!is.null(y)) {
    dname <- paste(deparse(substitute(x)), "and",
                   deparse(substitute(y)))
  }
  else {
    dname <- deparse(substitute(x))
  }

  # When tr == 0, use simple_htest for initial statistics (backward compat)
  # When tr > 0, compute initial stats directly since t.test() doesn't support trimming
  if (tr == 0) {
    null_test = simple_htest(x = x,
                             y = y,
                             test = "t.test",
                             var.equal = var.equal,
                             paired = paired,
                             alternative = alternative,
                             mu = mu,
                             alpha = 0.05)
    mu = null_test$null.value
  } else {
    # Handle mu for equivalence/minimal.effect the same way simple_htest does
    if (alternative %in% c("equivalence", "minimal.effect")) {
      if (length(mu) == 1) {
        if (mu == 0) {
          stop("mu cannot be zero if alternative is equivalence or minimal.effect")
        }
        mu = c(mu, -1 * mu)
      }
      mu = sort(mu)
      names(mu) = rep(paste0("trimmed mean difference (tr = ", tr, ")"), 2)
    } else {
      names(mu) = paste0("trimmed mean difference (tr = ", tr, ")")
    }
  }

  m_vec <- rep(NA, times=length(R)) # mean difference vector
  m_se_vec <- rep(NA, times=length(R)) # mean difference vector
  if(alternative %in% c("equivalence","minimal.effect")){
    conf.level = 1-alpha*2
  } else {
    conf.level = 1-alpha
  }

  # CI method label for method string
  ci_label <- switch(boot_ci,
                     "basic" = "(basic)",
                     "perc" = " (percentile)",
                     "bca" = " (BCa)",
                     "stud" = " (studentized)")

  if(!is.null(y)){
    dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
    if (paired) {
      i1 <- y
      i2 <- x
      data <- data.frame(i1 = i1, i2 = i2)
      data <- na.omit(data)
      y <- data$i1
      x <- data$i2
    }
    yok <- !is.na(y)
    xok <- !is.na(x)
    y <- y[yok]

  }else{
    dname <- deparse(substitute(x))
    xok <- !is.na(x)
    yok <- NULL
  }
  x <- x[xok]
  if(paired && !is.null(y)){
    x <- x - y
    y <- NULL
  }
  nx <- length(x)

  # Sample size check for trimming
  if (tr > 0) {
    g <- floor(tr * nx)
    effective_x <- nx - 2 * g
    if (effective_x < 2) {
      min_n <- ceiling(2 / (1 - 2 * tr)) + 1
      stop("Sample size too small for specified trimming proportion. ",
           "With tr = ", tr, ", need at least ", min_n,
           " observations, but only have ", nx, ".")
    }
  }

  # Compute location and scale for x
  if (tr > 0) {
    mx <- trimmed_mean(x, tr)
    vx <- winsorized_var(x, tr)
  } else {
    mx <- mean(x)
    vx <- var(x)
  }

  if (is.null(y)) {
    if (nx < 2)
      stop("not enough 'x' observations")

    if (tr > 0) {
      hx <- effective_n(nx, tr)
      df <- hx - 1
      stderr <- sqrt((nx - 1) * vx / (hx * (hx - 1)))
    } else {
      df <- nx - 1
      stderr <- sqrt(vx/nx)
    }

    if (stderr < 10 * .Machine$double.eps * abs(mx)){
      stop("data are essentially constant")
    }
    tstat <- (mx - mu)/stderr

    # Method string
    if (paired) {
      method <- if (tr > 0) "Bootstrapped Paired Yuen t-test" else "Bootstrapped Paired t-test"
    } else {
      method <- if (tr > 0) "Bootstrapped One Sample Yuen t-test" else "Bootstrapped One Sample t-test"
    }

    method = paste0(method, " ", ci_label)

    # Estimate labels
    XNAME <- "x"
    YNAME <- "y"
    if (tr > 0) {
      estimate <- setNames(mx, if (paired)
        paste0("trimmed mean of the differences (z = ", XNAME, " - ", YNAME, ", tr = ", tr, ")")
        else paste0("trimmed mean of ", XNAME))
    } else {
      estimate <- setNames(mx, if (paired)
        ttest_estimate_label(type = "t", xname = XNAME, yname = YNAME, paired = TRUE)
        else ttest_estimate_label(type = "t", xname = XNAME, yname = NULL, paired = FALSE))
    }

    if (tr == 0) {
      x.cent <- x - mx
      X <- matrix(sample(x.cent, size = nx*R, replace = TRUE), nrow = R)

      MX <- rowMeans(X)
      VX <- rowSums((X - MX) ^ 2) / (nx - 1)

      STDERR <- sqrt(VX/nx)
      TSTAT <- MX/STDERR
      EFF <- MX + mx

      for(i in 1:nrow(X)){
        dat = X[i,] + mx
        m_vec[i] <- mean(dat, na.rm=TRUE)
        m_se_vec[i] <- sd(dat, na.rm = TRUE)/sqrt(length(na.omit(dat)))
      }
    } else {
      # Trimmed path
      X <- matrix(sample(x, size = nx*R, replace = TRUE), nrow = R)
      TSTAT <- numeric(R)
      hx <- effective_n(nx, tr)

      for (i in 1:R) {
        dat <- X[i, ]
        mx_boot <- trimmed_mean(dat, tr)
        vx_boot <- winsorized_var(dat, tr)
        se_boot <- sqrt((nx - 1) * vx_boot / (hx * (hx - 1)))

        # Centered trimmed mean (under the null)
        mx_centered <- trimmed_mean(dat - mx, tr)

        if (se_boot > 10 * .Machine$double.eps) {
          TSTAT[i] <- mx_centered / se_boot
        } else {
          TSTAT[i] <- 0
        }

        m_vec[i] <- mx_boot
        m_se_vec[i] <- se_boot
      }
    }
  }

  if(!is.null(y) && !paired){
    ny <- length(y)

    # Sample size check for y with trimming
    if (tr > 0) {
      g_y <- floor(tr * ny)
      effective_y <- ny - 2 * g_y
      if (effective_y < 2) {
        min_n <- ceiling(2 / (1 - 2 * tr)) + 1
        stop("Sample size too small for specified trimming proportion. ",
             "With tr = ", tr, ", need at least ", min_n,
             " observations, but only have ", ny, ".")
      }
    }

    if(nx < 1 || (!var.equal && nx < 2))
      stop("not enough 'x' observations")
    if(ny < 1 || (!var.equal && ny < 2))
      stop("not enough 'y' observations")
    if(var.equal && nx + ny < 3)
      stop("not enough observations")

    if (tr > 0) {
      my <- trimmed_mean(y, tr)
      vy <- winsorized_var(y, tr)
    } else {
      my <- mean(y)
      vy <- var(y)
    }

    # Method string
    if (tr > 0) {
      method <- paste("Bootstrapped",
                      if (!var.equal) "Welch",
                      "Yuen Two Sample t-test")
      method <- gsub("\\s+", " ", method)
    } else {
      method <- paste("Bootstrapped", paste(if (!var.equal) "Welch", "Two Sample t-test "),
                      ci_label)
    }

    # Estimate labels
    XNAME <- "x"
    YNAME <- "y"
    if (tr > 0) {
      estimate <- c(mx, my, mx - my)
      names(estimate) <- c(paste0("trimmed mean of ", XNAME),
                            paste0("trimmed mean of ", YNAME),
                            paste0("trimmed mean difference (", XNAME, " - ", YNAME, ", tr = ", tr, ")"))
    } else {
      est_labels <- ttest_estimate_label(type = "t", xname = XNAME, yname = YNAME, paired = FALSE)
      estimate <- c(mx, my)
      names(estimate) <- est_labels
    }

    if(var.equal){
      if (tr > 0) {
        hx <- effective_n(nx, tr)
        hy <- effective_n(ny, tr)
        df <- pooled_wins_df(nx, ny, tr)
        v_pooled <- ((hx - 1) * vx + (hy - 1) * vy) / df
        stderr <- sqrt(v_pooled * ((nx - 1) / (hx * (hx - 1)) +
                                    (ny - 1) / (hy * (hy - 1))))
      } else {
        df <- nx + ny - 2
        v <- 0
        if (nx > 1){
          v <- v + (nx - 1) * vx
        }
        if (ny > 1){
          v <- v + (ny - 1) * vy
        }
        v <- v/df
        stderr <- sqrt(v * (1/nx + 1/ny))
      }

      z <- c(x, y)
      if (tr > 0) {
        mz <- trimmed_mean(z, tr)
      } else {
        mz <- mean(z)
      }

      if (tr == 0) {
        x.cent <- x - mx + mz
        y.cent <- y - my + mz
        X <- matrix(sample(x.cent, size = nx*R, replace = TRUE), nrow = R)
        Y <- matrix(sample(y.cent, size = ny*R, replace = TRUE), nrow = R)

        MX <- rowMeans(X)
        MY <- rowMeans(Y)
        V <- (rowSums((X-MX)^2) + rowSums((Y-MY)^2))/df
        STDERR <- sqrt(V*(1/nx + 1/ny))
        EFF <- (MX + mx) - (MY + my)

        for(i in 1:nrow(X)){
          dat_x = X[i,] + mx - mz
          dat_y = Y[i,] + my - mz
          m_vec[i] <- mean(dat_x, na.rm=TRUE) - mean(dat_y, na.rm=TRUE)
          m_se_vec[i] <- sqrt(sd(dat_x, na.rm=TRUE)^2/length(na.omit(dat_x)) + sd(dat_y, na.rm=TRUE)^2/length(na.omit(dat_y)))
        }
      } else {
        # Trimmed path - equal variance
        X <- matrix(sample(x, size = nx*R, replace = TRUE), nrow = R)
        Y <- matrix(sample(y, size = ny*R, replace = TRUE), nrow = R)
        TSTAT <- numeric(R)
        hx <- effective_n(nx, tr)
        hy <- effective_n(ny, tr)
        df_pool <- pooled_wins_df(nx, ny, tr)

        for (i in 1:R) {
          dat_x <- X[i, ]
          dat_y <- Y[i, ]

          mx_boot <- trimmed_mean(dat_x, tr)
          my_boot <- trimmed_mean(dat_y, tr)
          vx_boot <- winsorized_var(dat_x, tr)
          vy_boot <- winsorized_var(dat_y, tr)

          v_pooled_boot <- ((hx - 1) * vx_boot + (hy - 1) * vy_boot) / df_pool
          se_boot <- sqrt(v_pooled_boot * ((nx - 1) / (hx * (hx - 1)) +
                                            (ny - 1) / (hy * (hy - 1))))

          # Centered under the null
          mx_cent <- trimmed_mean(dat_x - mx + mz, tr)
          my_cent <- trimmed_mean(dat_y - my + mz, tr)

          if (se_boot > 10 * .Machine$double.eps) {
            TSTAT[i] <- (mx_cent - my_cent) / se_boot
          } else {
            TSTAT[i] <- 0
          }

          m_vec[i] <- mx_boot - my_boot
          m_se_vec[i] <- se_boot
        }
      }
    }else{
      if (tr > 0) {
        hx <- effective_n(nx, tr)
        hy <- effective_n(ny, tr)
        dx <- (nx - 1) * vx / (hx * (hx - 1))
        dy <- (ny - 1) * vy / (hy * (hy - 1))
        stderr <- sqrt(dx + dy)
        df <- yuen_welch_df(nx, ny, vx, vy, tr)
      } else {
        stderrx <- sqrt(vx/nx)
        stderry <- sqrt(vy/ny)
        stderr <- sqrt(stderrx^2 + stderry^2)
        df <- stderr^4/(stderrx^4/(nx - 1) + stderry^4/(ny - 1))
      }

      z <- c(x, y)
      if (tr > 0) {
        mz <- trimmed_mean(z, tr)
      } else {
        mz <- mean(z)
      }

      if (tr == 0) {
        x.cent <- x - mx + mz
        y.cent <- y - my + mz
        X <- matrix(sample(x.cent, size = nx*R, replace = TRUE), nrow = R)
        Y <- matrix(sample(y.cent, size = ny*R, replace = TRUE), nrow = R)

        MX <- rowMeans(X)
        MY <- rowMeans(Y)
        VX <- rowSums((X-MX)^2)/(nx-1)
        VY <- rowSums((Y-MY)^2)/(ny-1)
        STDERR <- sqrt(VX/nx + VY/ny)
        EFF <- (MX + mx) - (MY + my)

        for(i in 1:nrow(X)){
          dat_x = X[i,] + mx - mz
          dat_y = Y[i,] + my - mz
          m_vec[i] <- mean(dat_x, na.rm=TRUE) - mean(dat_y, na.rm=TRUE)
          m_se_vec[i] <- sqrt(sd(dat_x, na.rm=TRUE)^2/length(na.omit(dat_x)) + sd(dat_y, na.rm=TRUE)^2/length(na.omit(dat_y)))
        }
      } else {
        # Trimmed path - Welch/Yuen
        X <- matrix(sample(x, size = nx*R, replace = TRUE), nrow = R)
        Y <- matrix(sample(y, size = ny*R, replace = TRUE), nrow = R)
        TSTAT <- numeric(R)
        hx <- effective_n(nx, tr)
        hy <- effective_n(ny, tr)

        for (i in 1:R) {
          dat_x <- X[i, ]
          dat_y <- Y[i, ]

          mx_boot <- trimmed_mean(dat_x, tr)
          my_boot <- trimmed_mean(dat_y, tr)
          vx_boot <- winsorized_var(dat_x, tr)
          vy_boot <- winsorized_var(dat_y, tr)

          dx_boot <- (nx - 1) * vx_boot / (hx * (hx - 1))
          dy_boot <- (ny - 1) * vy_boot / (hy * (hy - 1))
          se_boot <- sqrt(dx_boot + dy_boot)

          mx_cent <- trimmed_mean(dat_x - mx + mz, tr)
          my_cent <- trimmed_mean(dat_y - my + mz, tr)

          if (se_boot > 10 * .Machine$double.eps) {
            TSTAT[i] <- (mx_cent - my_cent) / se_boot
          } else {
            TSTAT[i] <- 0
          }

          m_vec[i] <- mx_boot - my_boot
          m_se_vec[i] <- se_boot
        }
      }
    }
    if (stderr < 10 * .Machine$double.eps * max(abs(mx), abs(my))){
      stop("data are essentially constant")
    }

    tstat <- (mx - my - mu)/stderr

    if (tr == 0) {
      TSTAT <- (MX-MY)/STDERR
    }
  }

  if(is.null(y)){
    diff = mx
  }else{
    diff = mx-my
  }

  # Jackknife for BCa (if needed)
  if (boot_ci == "bca") {
    if (is.null(y)) {
      # One-sample or paired (paired already converted x <- x - y above)
      n_jack <- nx
      jack_est <- numeric(n_jack)
      for (j in seq_len(n_jack)) {
        jack_est[j] <- trimmed_mean(x[-j], tr)
      }
    } else {
      # Two-sample: pooled jackknife (delete one from combined)
      n_jack <- nx + ny
      jack_est <- numeric(n_jack)
      for (j in seq_len(nx)) {
        jack_est[j] <- trimmed_mean(x[-j], tr) - trimmed_mean(y, tr)
      }
      for (j in seq_len(ny)) {
        jack_est[nx + j] <- trimmed_mean(x, tr) - trimmed_mean(y[-j], tr)
      }
    }
  }

  if(alternative %in% c("equivalence", "minimal.effect")){
    tstat_l = (diff-min(mu))/stderr
    tstat_u = (diff-max(mu))/stderr
  }

  # se0 for studentized CI: use the observed stderr
  se0_val <- stderr

  # Pre-compute BCa parameters for p-value use
  z0 <- NULL; acc <- NULL
  if (boot_ci == "bca") {
    bca_par <- bca_params(m_vec, diff, jack_est)
    z0 <- bca_par$z0; acc <- bca_par$acc
  }

  # CI computation (alpha depends on alternative)
  ci_alpha <- if (alternative == "two.sided") alpha else alpha * 2
  boot.cint <- switch(boot_ci,
                      "stud" = stud(m_vec,
                                    boots_se = m_se_vec,
                                    t0 = diff,
                                    se0 = se0_val,
                                    alpha = ci_alpha),
                      "basic" = basic(m_vec, t0 = diff, ci_alpha),
                      "perc" = perc(m_vec, ci_alpha),
                      "bca" = bca_ci(boots_est = m_vec, t0 = diff,
                                     jack_est = jack_est, alpha = ci_alpha))

  # P-value computation (method-consistent with CI)
  if (alternative %in% c("two.sided", "greater", "less")) {
    boot.pval <- boot_pvalue(bvec = m_vec, est = diff, null = mu,
                             alternative = alternative, boot_ci = boot_ci,
                             tvec = TSTAT, se_obs = stderr,
                             z0 = z0, acc = acc, nboot = R)
  } else if (alternative == "equivalence") {
    p_l <- boot_pvalue(bvec = m_vec, est = diff, null = min(mu),
                       alternative = "greater", boot_ci = boot_ci,
                       tvec = TSTAT, se_obs = stderr,
                       z0 = z0, acc = acc, nboot = R)
    p_u <- boot_pvalue(bvec = m_vec, est = diff, null = max(mu),
                       alternative = "less", boot_ci = boot_ci,
                       tvec = TSTAT, se_obs = stderr,
                       z0 = z0, acc = acc, nboot = R)
    boot.pval <- max(p_l, p_u)
  } else if (alternative == "minimal.effect") {
    p_l <- boot_pvalue(bvec = m_vec, est = diff, null = max(mu),
                       alternative = "greater", boot_ci = boot_ci,
                       tvec = TSTAT, se_obs = stderr,
                       z0 = z0, acc = acc, nboot = R)
    p_u <- boot_pvalue(bvec = m_vec, est = diff, null = min(mu),
                       alternative = "less", boot_ci = boot_ci,
                       tvec = TSTAT, se_obs = stderr,
                       z0 = z0, acc = acc, nboot = R)
    boot.pval <- min(p_l, p_u)
  }

  boot.se = sd(m_vec, na.rm = TRUE)
  attr(boot.cint, "conf.level") <- conf.level

  # Name the statistic and parameter
  # For equivalence/minimal.effect, report the t-statistic corresponding to the p-value
  if (alternative == "equivalence") {
    # p-value is max of the two one-sided p-values
    # Report the t-statistic for the "binding" bound (the one with higher p-value)
    if (p_l >= p_u) {
      tstat_report <- tstat_l
    } else {
      tstat_report <- tstat_u
    }
  } else if (alternative == "minimal.effect") {
    # p-value is min of the two one-sided p-values
    # Report the t-statistic for the "binding" bound (the one with lower p-value)
    if (p_l <= p_u) {
      tstat_report <- tstat_l
    } else {
      tstat_report <- tstat_u
    }
  } else {
    tstat_report <- tstat
  }

  names(tstat_report) <- "t-observed"
  names(df) <- "df"

  # Build estimate and data.name from appropriate sources
  if (tr > 0) {
    rval_estimate <- estimate
    rval_null <- mu
    rval_dname <- dname
  } else {
    rval_estimate <- null_test$estimate
    rval_null <- null_test$null.value
    rval_dname <- null_test$data.name
  }

  # Compute sample_size
  if (tr > 0) {
    # When tr > 0, compute directly from data
    # Note: for paired, x <- x - y has already been done so y is NULL
    if (is.null(y)) {
      sample_size <- c(n = nx)
    } else {
      sample_size <- c(nx = nx, ny = ny)
    }
  } else {
    # When tr == 0, inherit from simple_htest which now includes sample_size
    sample_size <- null_test$sample_size
  }

  rval = list(
    statistic = tstat_report,
    parameter = df,
    p.value = boot.pval,
    stderr = boot.se,
    conf.int = boot.cint,
    estimate = rval_estimate,
    null.value = rval_null,
    alternative = alternative,
    method = method,
    boot = m_vec,
    data.name = rval_dname,
    call = match.call(),
    sample_size = sample_size
  )

  class(rval) = "htest"

  return(rval)
}

#' @rdname boot_t_test
#' @method boot_t_test formula
#' @export
#'
boot_t_test.formula <- function (formula, data, subset, na.action, ...){
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
  y <- do.call("boot_t_test", c(DATA, list(...)))
  y$data.name <- DNAME

  # Resolve actual group labels from factor levels
  XNAME <- levels(g)[1]
  YNAME <- levels(g)[2]
  xq <- quote_if_numeric(XNAME)
  yq <- quote_if_numeric(YNAME)

  is_paired <- isTRUE(dots$paired)
  tr_val <- if (!is.null(dots$tr)) dots$tr else 0

  if (tr_val > 0) {
    # Trimmed labels: substitute group names
    nms <- names(y$estimate)
    nms <- gsub("\\bof x\\b", paste0("of ", xq), nms)
    nms <- gsub("\\bof y\\b", paste0("of ", yq), nms)
    nms <- gsub("\\(x - y", paste0("(", xq, " - ", yq), nms)
    nms <- gsub("\\(z = x - y", paste0("(z = ", xq, " - ", yq), nms)
    names(y$estimate) <- nms
  } else {
    if (!is_paired && length(y$estimate) >= 2) {
      est_labels <- ttest_estimate_label(type = "t", xname = XNAME, yname = YNAME, paired = FALSE)
      diff_label <- paste0("mean difference (", xq, " - ", yq, ")")
      names(y$estimate) <- c(est_labels, diff_label)
    } else if (is_paired) {
      names(y$estimate) <- ttest_estimate_label(type = "t", xname = XNAME, yname = YNAME, paired = TRUE)
    }
  }

  # Relabel sample_size names
  if (!is.null(y$sample_size) && length(y$sample_size) == 2) {
    names(y$sample_size) <- levels(g)
  }

  y
}
