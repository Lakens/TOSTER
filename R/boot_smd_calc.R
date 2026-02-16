#' @title Bootstrapped Standardized Mean Difference (SMD) Calculation
#' @description
#' `r lifecycle::badge('maturing')`
#'
#' Calculates standardized mean differences (SMDs) with bootstrap confidence intervals,
#' with optional hypothesis testing.
#'
#' @section Purpose:
#' Use this function when:
#'
#'   * You need more robust confidence intervals for standardized mean differences
#'   * You want to account for non-normality or heterogeneity in your effect size estimates
#'   * Sample sizes are small or standard error approximations may be unreliable
#'   * You prefer resampling-based confidence intervals over parametric approximations
#'   * You need to quantify uncertainty in SMD estimates more accurately
#'   * You want to test hypotheses about effect size magnitudes using bootstrap methods
#'
#' @inheritParams boot_t_TOST
#' @param mu null value to adjust the calculation. If non-zero, the function calculates x-y-mu (default = 0).
#' @param ... further arguments to be passed to or from methods.
#' @param output a character string specifying the output format:
#'     - "htest": (default) Returns an object of class "htest" compatible with standard R output.
#'     - "data.frame": Returns a data frame for backward compatibility.
#' @param null.value a number or vector specifying the null hypothesis value(s) on the SMD scale:
#'     - For standard alternatives: a single value (default = 0)
#'     - For equivalence/minimal.effect: two values representing the lower and upper bounds
#' @param alternative a character string specifying the alternative hypothesis:
#'     - "none": (default) No hypothesis test is performed; only effect size and CI are returned.
#'     - "two.sided": Test whether SMD differs from null.value
#'     - "less": Test whether SMD is less than null.value
#'     - "greater": Test whether SMD is greater than null.value
#'     - "equivalence": Test whether SMD is between specified bounds
#'     - "minimal.effect": Test whether SMD is outside specified bounds
#' @param denom a character string specifying the denominator for standardization:
#'     - "auto": (default) Uses the standard denominator based on design and other arguments
#'       (glass, rm_correction, var.equal).
#'     - "z": SD of differences (Cohen's d_z). Valid for paired and one-sample designs.
#'     - "rm": Repeated-measures corrected (Cohen's d_rm). Valid for paired designs only.
#'     - "pooled": Pooled SD (Cohen's d_s). Valid for independent samples only.
#'     - "avg": Root-mean-square SD (Cohen's d_av). Valid for independent samples only.
#'     - "glass1": First group's (x) SD (Glass's delta). Valid for paired and independent designs.
#'     - "glass2": Second group's (y) SD (Glass's delta). Valid for paired and independent designs.
#'
#'     When set to any value other than "auto", this overrides the glass, rm_correction,
#'     and var.equal arguments. The bias_correction argument is not affected.
#'
#' @details
#' This function calculates bootstrapped confidence intervals for standardized mean differences.
#' It is an extension of the `smd_calc()` function that uses resampling to provide more robust
#' confidence intervals, especially for small sample sizes or when data violate assumptions
#' of parametric methods.
#'
#' The function implements the following bootstrap approach:
#'   * Calculate the raw SMD and its standard error using the original data
#'   * Create R bootstrap samples by resampling with replacement from the original data
#'   * Calculate the SMD and its standard error for each bootstrap sample
#'   * Calculate confidence intervals using the specified method
#'
#' Three bootstrap confidence interval methods are available:
#'   - **Studentized bootstrap ("stud")**: Accounts for the variability in standard error estimates. Usually provides the most accurate coverage probability and is set as the default.
#'   - **Basic bootstrap ("basic")**: Uses the empirical distribution of bootstrap estimates. Simple approach that works well for symmetric distributions.
#'   - **Percentile bootstrap ("perc")**: Uses percentiles of the bootstrap distribution directly. More robust to skewness in the bootstrap distribution.
#'
#' The function supports various SMD variants:
#'   * Classic standardized mean difference (bias_correction = FALSE)
#'   * Bias-corrected version (bias_correction = TRUE)
#'   * Glass's delta: Uses only one group's standard deviation as the denominator (glass = "glass1" or "glass2")
#'   * Repeated measures d: Accounts for correlation in paired designs (rm_correction = TRUE)
#'
#' The function supports three study designs:
#'   * One-sample design: Standardizes the difference between the sample mean and zero (or other specified value)
#'   * Two-sample independent design: Standardizes the difference between two group means
#'   * Paired samples design: Standardizes the mean difference between paired observations
#'
#' The `denom` parameter provides a direct way to select the standardization denominator.
#' When `denom` is not "auto", it takes precedence over the `glass`, `rm_correction`, and
#' `var.equal` arguments, which are overridden as needed. A message is emitted if any
#' explicitly provided arguments are overridden. The `bias_correction` argument is always
#' respected regardless of `denom`.
#'
#' For detailed information on calculation methods, see `vignette("SMD_calcs")`.
#'
#' @return
#' If `output = "htest"` (default), returns a list with class `"htest"` containing:
#'   - estimate: The SMD estimate (Cohen's d, Hedges' g, or Glass's delta)
#'   - stderr: Standard error estimated from the bootstrap distribution
#'   - conf.int: Bootstrap confidence interval with conf.level attribute
#'   - alternative: A character string describing the alternative hypothesis
#'   - method: A character string indicating what type of test was performed
#'   - note: A character string describing the bootstrap CI method used
#'   - boot: The bootstrap distribution of SMD estimates
#'   - data.name: A character string giving the name(s) of the data
#'   - call: The matched call
#'   - statistic: z-statistic (only if alternative != "none")
#'   - p.value: Bootstrap p-value (only if alternative != "none")
#'   - null.value: The specified hypothesized value(s) (only if alternative != "none")
#'
#' If `output = "data.frame"`, returns a data frame containing:
#'   - estimate: The SMD calculated from the original data
#'   - bias: Estimated bias (difference between original estimate and median of bootstrap estimates)
#'   - SE: Standard error estimated from the bootstrap distribution
#'   - lower.ci: Lower bound of the bootstrap confidence interval
#'   - upper.ci: Upper bound of the bootstrap confidence interval
#'   - conf.level: Confidence level (1-alpha)
#'   - boot_ci: The bootstrap confidence interval method used
#'
#' @examples
#' # Example 1: Independent groups comparison with studentized bootstrap CI
#' set.seed(123)
#' group1 <- rnorm(30, mean = 100, sd = 15)
#' group2 <- rnorm(30, mean = 110, sd = 18)
#'
#' # Use fewer bootstrap replicates for a quick example
#' result <- boot_smd_calc(x = group1, y = group2,
#'                       boot_ci = "stud",
#'                       R = 999)
#'
#' # Example 2: Using formula notation with basic bootstrap and Hedges' g
#' df <- data.frame(
#'   value = c(group1, group2),
#'   group = factor(rep(c("A", "B"), each = 30))
#' )
#' result <- boot_smd_calc(formula = value ~ group,
#'                       data = df,
#'                       boot_ci = "basic",
#'                       bias_correction = TRUE,
#'                       R = 999)
#'
#' # Example 3: Paired samples with percentile bootstrap
#' set.seed(456)
#' before <- rnorm(30)
#' after <- rnorm(30)
#' result <- boot_smd_calc(x = before,
#'                       y = after,
#'                       paired = TRUE,
#'                       boot_ci = "perc",
#'                       R = 999)
#'
#' # Example 4: Glass's delta with homogeneous variances
#' set.seed(456)
#' control <- rnorm(25, mean = 50, sd = 10)
#' treatment <- rnorm(25, mean = 60, sd = 10)
#' result <- boot_smd_calc(x = control,
#'                       y = treatment,
#'                       glass = "glass1",
#'                       boot_ci = "stud",
#'                       R = 999)
#'
#' # Example 5: Two-sided hypothesis test
#' result <- boot_smd_calc(x = group1, y = group2,
#'                       alternative = "two.sided",
#'                       null.value = 0,
#'                       R = 999)
#'
#' # Example 6: Equivalence test with bootstrap
#' result <- boot_smd_calc(x = group1, y = group2,
#'                       alternative = "equivalence",
#'                       null.value = c(-0.5, 0.5),
#'                       R = 999)
#'
#' # Example 7: Legacy data.frame output
#' result <- boot_smd_calc(x = group1, y = group2,
#'                       output = "data.frame",
#'                       R = 999)
#'
#' @family effect sizes
#' @name boot_smd_calc
#' @export boot_smd_calc

# Bootstrap -------

#smd_calc <- setClass("smd_calc")
boot_smd_calc <- function(x, ...,
                          paired = FALSE,
                          var.equal = FALSE,
                          alpha = 0.05,
                          bias_correction = TRUE,
                          rm_correction = FALSE,
                          glass = NULL,
                          denom = c("auto", "z", "rm", "pooled", "avg",
                                    "glass1", "glass2"),
                          boot_ci = c("stud","basic","perc","bca"),
                          R = 1999,
                          output = c("htest", "data.frame"),
                          null.value = 0,
                          alternative = c("none", "two.sided", "less", "greater",
                                          "equivalence", "minimal.effect")){
  UseMethod("boot_smd_calc")
}



#' @rdname boot_smd_calc
#' @method boot_smd_calc default
#' @export
boot_smd_calc.default = function(x,
                                 y = NULL,
                                 paired = FALSE,
                                 var.equal = FALSE,
                                 alpha = 0.05,
                                 mu = 0,
                                 bias_correction = TRUE,
                                 rm_correction = FALSE,
                                 glass = NULL,
                                 denom = c("auto", "z", "rm", "pooled", "avg",
                                           "glass1", "glass2"),
                                 boot_ci = c("stud","basic","perc","bca"),
                                 R = 1999,
                                 output = c("htest", "data.frame"),
                                 null.value = 0,
                                 alternative = c("none", "two.sided", "less", "greater",
                                                 "equivalence", "minimal.effect"),
                                 ...) {
  denom = match.arg(denom)
  boot_ci = match.arg(boot_ci)
  output = match.arg(output)
  alternative = match.arg(alternative)

  # Capture explicit-pass status before any modifications
  var.equal_explicit <- !missing(var.equal)
  rm_correction_explicit <- !missing(rm_correction)
  glass_explicit <- !missing(glass)

  # Determine sample_type early for denom resolution
  if(is.null(y)){
    sample_type = "One Sample"
  } else if(paired == TRUE) {
    sample_type = "Paired Sample"
  } else {
    sample_type = "Two Sample"
  }

  # Resolve denom ONCE - messages emitted here, not in loop
  resolved <- resolve_denom(
    denom = denom,
    sample_type = sample_type,
    var.equal = var.equal,
    rm_correction = rm_correction,
    glass = if (glass_explicit) glass else NULL,
    var.equal_explicit = var.equal_explicit,
    rm_correction_explicit = rm_correction_explicit,
    glass_explicit = glass_explicit
  )

  for (msg in resolved$messages) message(msg)

  # Apply resolved values for all subsequent smd_calc calls
  var.equal <- resolved$var.equal
  rm_correction <- resolved$rm_correction
  glass <- resolved$glass

  # Handle equivalence/minimal.effect bounds
  if (alternative %in% c("equivalence", "minimal.effect")) {
    if (length(null.value) != 2) {
      stop("For equivalence or minimal.effect testing, null.value must be a vector of two values (lower and upper bounds)")
    }
    low_bound <- min(null.value)
    high_bound <- max(null.value)
    conf.level <- 1 - alpha * 2
  } else {
    if (length(null.value) > 1) {
      warning("null.value has length > 1; only the first element will be used")
      null.value <- null.value[1]
    }
    low_bound <- null.value
    high_bound <- null.value
    conf.level <- 1 - alpha
  }

  # Get data name
  if (!is.null(y)) {
    dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
  } else {
    dname <- deparse(substitute(x))
  }

  if(paired == TRUE && !missing(y)){
    i1 <- x
    i2 <- y
    data <- data.frame(x = i1, y = i2)
    data <- na.omit(data)
    raw_smd = smd_calc(x = data$x,
                       y = data$y,
                       paired = paired,
                       var.equal = var.equal,
                       alpha = alpha,
                       mu = mu,
                       bias_correction = bias_correction,
                       rm_correction = rm_correction,
                       glass = glass,
                       smd_ci = "z",
                       output = "data.frame")

    boots = c()
    boots_se = c()

    for(i in 1:R){
      sampler = sample(1:nrow(data), replace = TRUE)
      res_boot = smd_calc(x = data$x[sampler],
                          y = data$y[sampler],
                          paired = paired,
                          var.equal = var.equal,
                          alpha = alpha,
                          mu = mu,
                          bias_correction = bias_correction,
                          rm_correction = rm_correction,
                          glass = glass,
                          smd_ci = "z",
                          output = "data.frame")
      boots = c(boots, res_boot$estimate)
      boots_se = c(boots_se, res_boot$SE)
    }


  } else if(!missing(y)){

    i1 <- na.omit(x)
    i2 <- na.omit(y)
    data <- data.frame(values = c(i1,i2),
                       group = c(rep("x",length(i1)),
                                 rep("y",length(i2))))

    raw_smd = smd_calc(x = i1,
                       y = i2,
                       paired = paired,
                       var.equal = var.equal,
                       alpha = alpha,
                       mu = mu,
                       bias_correction = bias_correction,
                       rm_correction = rm_correction,
                       glass = glass,
                       smd_ci = "z",
                       output = "data.frame")

    boots = c()
    boots_se = c()
    for(i in 1:R){
      sampler = sample(1:nrow(data), replace = TRUE)
      boot_dat = data[sampler,]
      x_boot = subset(boot_dat,
                      group == "x")
      y_boot = subset(boot_dat,
                      group == "y")
      res_boot = smd_calc(x = x_boot$values,
                          y = y_boot$values,
                          paired = paired,
                          var.equal = var.equal,
                          alpha = alpha,
                          mu = mu,
                          bias_correction = bias_correction,
                          rm_correction = rm_correction,
                          glass = glass,
                          smd_ci = "z",
                          output = "data.frame")
      boots = c(boots, res_boot$estimate)
      boots_se = c(boots_se, res_boot$SE)
    }

  } else {

    x1 = na.omit(x)
    n1 = length(x1)
    raw_smd = smd_calc(x = x1,
                       paired = paired,
                       var.equal = var.equal,
                       alpha = alpha,
                       mu = mu,
                       bias_correction = bias_correction,
                       rm_correction = rm_correction,
                       glass = glass,
                       smd_ci = "z",
                       output = "data.frame")

    boots = c()
    boots_se = c()
    for(i in 1:R){
      sampler = sample(1:length(x1), replace = TRUE)
      x_boot = x1[sampler]

      res_boot = smd_calc(x = x_boot,
                          paired = paired,
                          var.equal = var.equal,
                          alpha = alpha,
                          mu = mu,
                          bias_correction = bias_correction,
                          rm_correction = rm_correction,
                          glass = glass,
                          smd_ci = "z",
                          output = "data.frame")
      boots = c(boots, res_boot$estimate)
      boots_se = c(boots_se, res_boot$SE)
    }

  }

  # Jackknife for BCa (if needed)
  if (boot_ci == "bca") {
    if (paired == TRUE && !missing(y)) {
      # Paired: delete one pair at a time
      n_jack <- nrow(data)
      jack_est <- numeric(n_jack)
      for (j in seq_len(n_jack)) {
        res_jack <- smd_calc(x = data$x[-j],
                             y = data$y[-j],
                             paired = paired,
                             var.equal = var.equal,
                             alpha = alpha,
                             mu = mu,
                             bias_correction = bias_correction,
                             rm_correction = rm_correction,
                             glass = glass,
                             smd_ci = "z",
                             output = "data.frame")
        jack_est[j] <- res_jack$estimate
      }
    } else if (!missing(y)) {
      # Two-sample: pooled jackknife (delete one from combined)
      n1 <- length(i1)
      n2 <- length(i2)
      n_total <- n1 + n2
      jack_est <- numeric(n_total)
      for (j in seq_len(n1)) {
        res_jack <- smd_calc(x = i1[-j],
                             y = i2,
                             paired = paired,
                             var.equal = var.equal,
                             alpha = alpha,
                             mu = mu,
                             bias_correction = bias_correction,
                             rm_correction = rm_correction,
                             glass = glass,
                             smd_ci = "z",
                             output = "data.frame")
        jack_est[j] <- res_jack$estimate
      }
      for (j in seq_len(n2)) {
        res_jack <- smd_calc(x = i1,
                             y = i2[-j],
                             paired = paired,
                             var.equal = var.equal,
                             alpha = alpha,
                             mu = mu,
                             bias_correction = bias_correction,
                             rm_correction = rm_correction,
                             glass = glass,
                             smd_ci = "z",
                             output = "data.frame")
        jack_est[n1 + j] <- res_jack$estimate
      }
    } else {
      # One-sample: delete one observation at a time
      n_jack <- length(x1)
      jack_est <- numeric(n_jack)
      for (j in seq_len(n_jack)) {
        res_jack <- smd_calc(x = x1[-j],
                             paired = paired,
                             var.equal = var.equal,
                             alpha = alpha,
                             mu = mu,
                             bias_correction = bias_correction,
                             rm_correction = rm_correction,
                             glass = glass,
                             smd_ci = "z",
                             output = "data.frame")
        jack_est[j] <- res_jack$estimate
      }
    }
  }

  ci_alpha <- if(alternative %in% c("equivalence", "minimal.effect")) alpha * 2 else alpha
  ci = switch(boot_ci,
              "stud" = stud(boots_est = boots, boots_se = boots_se,
                            se0=raw_smd$SE[1L], t0 = raw_smd$estimate[1L],
                            alpha = ci_alpha),
              "perc" = perc(boots, ci_alpha),
              "basic" = basic(boots, t0 = raw_smd$estimate[1L],
                              alpha = ci_alpha),
              "bca" = bca_ci(boots_est = boots, t0 = raw_smd$estimate[1L],
                             jack_est = jack_est, alpha = ci_alpha))

  # Determine SMD label
  smd_label <- if (bias_correction) {
    if (!is.null(glass) && glass %in% c("glass1", "glass2")) {
      paste0("Glass's ", ifelse(glass == "glass1", "delta1", "delta2"))
    } else {
      "Hedges' g"
    }
  } else {
    "Cohen's d"
  }

  # Bootstrap SE (for reporting in stderr field)
  boot_se <- sd(boots, na.rm = TRUE)

  # Studentized bootstrap pivot: centered at observed estimate, scaled by bootstrap SE
  # Analogous to TSTAT in boot_t_test
  TSTAT <- (boots - raw_smd$estimate[1L]) / boots_se

  # Compute p-value using studentized bootstrap distribution
  if (alternative != "none") {
    est_val <- raw_smd$estimate[1L]
    se_obs <- raw_smd$SE[1L]

    if (alternative %in% c("equivalence", "minimal.effect")) {
      z_obs_l <- (est_val - low_bound) / se_obs
      z_obs_u <- (est_val - high_bound) / se_obs
    } else {
      z_obs <- (est_val - null.value) / se_obs
    }

    if (alternative == "two.sided") {
      boot.pval <- 2 * min(mean(TSTAT <= z_obs, na.rm = TRUE),
                           mean(TSTAT > z_obs, na.rm = TRUE))

    } else if (alternative == "less") {
      boot.pval <- mean(TSTAT < z_obs, na.rm = TRUE)

    } else if (alternative == "greater") {
      boot.pval <- mean(TSTAT > z_obs, na.rm = TRUE)

    } else if (alternative == "equivalence") {
      p_l <- mean(TSTAT > z_obs_l, na.rm = TRUE)
      p_u <- mean(TSTAT < z_obs_u, na.rm = TRUE)
      boot.pval <- max(p_l, p_u)

    } else if (alternative == "minimal.effect") {
      p_l <- mean(TSTAT < z_obs_l, na.rm = TRUE)
      p_u <- mean(TSTAT > z_obs_u, na.rm = TRUE)
      boot.pval <- min(p_l, p_u)
    }

    # Report z-observed (analogous to t-observed in boot_t_test)
    if (alternative %in% c("equivalence", "minimal.effect")) {
      if (alternative == "equivalence") {
        z_stat <- if (p_l >= p_u) z_obs_l else z_obs_u
      } else {
        z_stat <- if (p_l <= p_u) z_obs_l else z_obs_u
      }
    } else {
      z_stat <- z_obs
    }
  }

  # Build output
  if (output == "data.frame") {
    effsize = data.frame(
      estimate = raw_smd$estimate,
      bias = raw_smd$estimate - median(boots, na.rm=TRUE),
      SE = sd(boots),
      lower.ci = ci[1],
      upper.ci = ci[2],
      conf.level = conf.level,
      boot_ci = boot_ci,
      row.names = row.names(raw_smd)
    )
    return(effsize)

  } else {
    # htest output
    estimate <- raw_smd$estimate
    names(estimate) <- row.names(raw_smd)

    conf.int <- c(ci[1], ci[2])
    attr(conf.int, "conf.level") <- conf.level

    method_suffix <- if (alternative != "none") "test" else "estimate with CI"
    method_desc <- paste0("Bootstrapped ", sample_type, " ", smd_label, " ", method_suffix)

    note_text <- paste0("Bootstrap CI: ", boot_ci)

    rval <- list(
      estimate = estimate,
      stderr = boot_se,
      conf.int = conf.int,
      alternative = alternative,
      method = method_desc,
      note = note_text,
      boot = boots,
      data.name = dname,
      call = match.call()
    )

    if (alternative != "none") {
      names(z_stat) <- "z-observed"

      if (alternative %in% c("equivalence", "minimal.effect")) {
        null_val <- c(low_bound, high_bound)
        names(null_val) <- c("lower bound", "upper bound")
      } else {
        null_val <- null.value
        names(null_val) <- smd_label
      }

      rval$statistic <- z_stat
      rval$p.value <- boot.pval
      rval$null.value <- null_val
    }

    class(rval) <- "htest"
    return(rval)
  }

}

#' @rdname boot_smd_calc
#' @method boot_smd_calc formula
#' @export

boot_smd_calc.formula = function(formula,
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
  y <- do.call("boot_smd_calc", c(DATA, list(...)))

  # Update data.name for htest output
  if (inherits(y, "htest")) {
    y$data.name <- DNAME
  }

  y

}
