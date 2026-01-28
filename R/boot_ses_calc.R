#' @title Bootstrapped Standardized Effect Size (SES) Calculation
#' @description
#' `r lifecycle::badge('maturing')`
#'
#' Calculates non-SMD standardized effect sizes with bootstrap confidence intervals and
#' optional hypothesis testing. This function provides robust confidence intervals for
#' rank-based and probability-based effect size measures through resampling methods.
#'
#' @section Purpose:
#' Use this function when:
#'   - You need more robust confidence intervals for non-parametric effect sizes
#'   - You prefer resampling-based confidence intervals over asymptotic approximations
#'   - You need to quantify uncertainty in rank-based effect sizes more accurately
#'   - You want to perform hypothesis testing with bootstrap methods
#'
#' @inheritParams wilcox_TOST
#' @inheritParams boot_t_TOST
#' @inheritParams boot_smd_calc
#' @param ses a character string specifying the effect size measure to calculate:
#'     - "rb": rank-biserial correlation (default)
#'     - "odds": Wilcoxon-Mann-Whitney odds
#'     - "logodds": Wilcoxon-Mann-Whitney log-odds
#'     - "cstat": concordance statistic (C-statistic/AUC)
#' @param mu number indicating the value around which asymmetry (for one-sample or paired samples)
#'   or shift (for independent samples) is to be estimated (default = 0).
#' @param se_method a character string specifying the method for computing standard errors
#'   within each bootstrap sample:
#'     - "agresti": (default) Uses the Agresti/Lehmann placement-based variance estimation.
#'       This method has better asymptotic properties.
#'     - "fisher": Uses the legacy Fisher z-transformation method. Retained for backward
#'       compatibility.
#' @param output a character string specifying the output format:
#'     - "htest": (default) Returns an object of class "htest" compatible with standard R output.
#'     - "data.frame": Returns a data frame with effect size estimates and confidence intervals.
#' @param alternative a character string specifying the alternative hypothesis:
#'     - "two.sided": effect differs from null.value (default)
#'     - "less": effect is less than null.value
#'     - "greater": effect is greater than null.value
#'     - "equivalence": effect is between specified bounds
#'     - "minimal.effect": effect is outside specified bounds
#' @param null.value a number or vector specifying the null hypothesis value(s):
#'     - For standard alternatives: a single value (default = 0 for rb/logodds, 0.5 for cstat, 1 for odds)
#'     - For equivalence/minimal.effect: two values representing the lower and upper bounds
#' @param ... further arguments to be passed to or from methods.
#'
#' @details
#' This function calculates bootstrapped confidence intervals for rank-based and probability-based
#' effect size measures. It extends the `ses_calc()` function by using resampling
#' to provide more robust confidence intervals, especially for small sample sizes.
#'
#' The function implements the following bootstrap approach:
#'   - Calculate the raw effect size using the original data
#'   - Create R bootstrap samples by resampling with replacement from the original data
#'   - Calculate the effect size for each bootstrap sample
#'   - Apply the Fisher z-transformation to stabilize variance for rank-biserial correlation values
#'   - Calculate confidence intervals using the specified method
#'   - Back-transform the confidence intervals to the original scale
#'   - Convert to the requested effect size measure (if not rank-biserial)
#'
#' ## Standard Error Methods
#'
#' Two methods are available for computing standard errors within each bootstrap sample:
#'
#'   - **Agresti method** (`se_method = "agresti"`): Uses the Agresti/Lehmann placement-based
#'     variance estimation. This method computes the variance of the concordance probability
#'     and propagates it to other effect size scales using the delta method.
#'
#'   - **Fisher method** (`se_method = "fisher"`): Uses the legacy formula based on the
#'     Wilcoxon statistic variance. This is retained for backward compatibility.
#'
#' ## Bootstrap Confidence Interval Methods
#'
#' Three bootstrap confidence interval methods are available:
#'   - **Basic bootstrap ("basic")**: Uses the empirical distribution of bootstrap estimates
#'   - **Studentized bootstrap ("stud")**: Accounts for the variability in standard error estimates
#'   - **Percentile bootstrap ("perc")**: Uses percentiles of the bootstrap distribution directly
#'
#' ## Hypothesis Testing
#'
#' When an alternative other than "two.sided" is specified, or when null.value is not the
#' default, the function performs bootstrap hypothesis testing. For equivalence and minimal
#' effect testing, specify null.value as a vector of two values (lower and upper bounds).
#'
#' For different alternatives, the p-values are calculated as follows:
#'   * "two.sided": Proportion of bootstrap statistics at least as extreme as the observed statistic
#'   * "less": Proportion of bootstrap statistics less than or equal to the observed statistic
#'   * "greater": Proportion of bootstrap statistics greater than or equal to the observed statistic
#'   * "equivalence": Maximum of two one-sided p-values (for lower and upper bounds)
#'   * "minimal.effect": Minimum of two one-sided p-values (for lower and upper bounds)
#'
#' The function supports three study designs:
#'   - One-sample design: Compares a single sample to a specified value
#'   - Two-sample independent design: Compares two independent groups
#'   - Paired samples design: Compares paired observations
#'
#' Note that extreme values (perfect separation between groups) can produce infinite values during
#' the bootstrapping process. The function will issue a warning if this occurs.
#'
#' For detailed information on calculation methods, see `vignette("robustTOST")`.
#'
#' @return
#' If `output = "htest"` (default), returns a list with class `"htest"` containing:
#'   - estimate: The effect size estimate calculated from the original data
#'   - stderr: Standard error estimated from the bootstrap distribution
#'   - conf.int: Bootstrap confidence interval with conf.level attribute
#'   - statistic: The observed z-statistic (for hypothesis testing)
#'   - p.value: The bootstrapped p-value for the test
#'   - null.value: The specified hypothesized value(s)
#'   - alternative: A character string describing the alternative hypothesis
#'   - method: A character string indicating what type of test was performed
#'   - boot: The bootstrap samples of the effect size (on the requested scale)
#'   - data.name: A character string giving the name(s) of the data
#'   - call: The matched call
#'
#' If `output = "data.frame"`, returns a data frame containing:
#'   - estimate: The effect size estimate
#'   - SE: Standard error from the bootstrap distribution
#'   - lower.ci: Lower bound of the bootstrap confidence interval
#'   - upper.ci: Upper bound of the bootstrap confidence interval
#'   - conf.level: Confidence level (1-alpha or 1-2*alpha for equivalence)
#'   - boot_ci: The bootstrap CI method used
#'
#' @examples
#' # Example 1: Independent groups comparison with basic bootstrap CI
#' set.seed(123)
#' group1 <- c(1.2, 2.3, 3.1, 4.6, 5.2, 6.7)
#' group2 <- c(3.5, 4.8, 5.6, 6.9, 7.2, 8.5)
#'
#' # Use fewer bootstrap replicates for a quick example
#' result <- boot_ses_calc(x = group1, y = group2,
#'                         ses = "rb",
#'                         boot_ci = "basic",
#'                         R = 99)
#'
#' # Example 2: Hypothesis testing (two-sided)
#' result <- boot_ses_calc(x = group1, y = group2,
#'                         ses = "rb",
#'                         alternative = "two.sided",
#'                         null.value = 0,
#'                         R = 99)
#'
#' # Example 3: Equivalence testing
#' result <- boot_ses_calc(x = group1, y = group2,
#'                         ses = "rb",
#'                         alternative = "equivalence",
#'                         null.value = c(-0.3, 0.3),
#'                         R = 99)
#'
#' # Example 4: Paired samples
#' data(sleep)
#' with(sleep, boot_ses_calc(x = extra[group == 1],
#'                           y = extra[group == 2],
#'                           paired = TRUE,
#'                           ses = "rb",
#'                           alternative = "greater",
#'                           R = 99))
#'
#' # Example 5: Using formula notation
#' data(mtcars)
#' result <- boot_ses_calc(formula = mpg ~ am,
#'                         data = mtcars,
#'                         ses = "cstat",
#'                         boot_ci = "perc",
#'                         R = 99)
#'
#' @references
#' Agresti, A. (1980). Generalized odds ratios for ordinal data. *Biometrics*, 36, 59-67.
#'
#' Bamber, D. (1975). The area above the ordinal dominance graph and the area below the receiver
#' operating characteristic graph. *Journal of Mathematical Psychology*, 12, 387-415.
#'
#' Efron, B., & Tibshirani, R. J. (1994). An introduction to the bootstrap. CRC press.
#'
#' Lehmann, E.L. (1975). *Nonparametrics: Statistical Methods Based on Ranks*. Holden-Day.
#'
#' @family effect sizes
#' @name boot_ses_calc
#' @export boot_ses_calc

boot_ses_calc <- function(x, ...,
                          paired = FALSE,
                          ses = "rb",
                          alpha = 0.05,
                          mu = 0,
                          boot_ci = c("basic","stud","perc"),
                          R = 1999,
                          se_method = c("agresti", "fisher"),
                          output = c("htest", "data.frame"),
                          alternative = c("two.sided", "less", "greater",
                                          "equivalence", "minimal.effect"),
                          null.value = NULL){
  UseMethod("boot_ses_calc")
}


#' @rdname boot_ses_calc
#' @method boot_ses_calc default
#' @export
boot_ses_calc.default = function(x,
                                 y = NULL,
                                 paired = FALSE,
                                 ses = c("rb","odds","logodds","cstat"),
                                 alpha = 0.05,
                                 mu = 0,
                                 boot_ci = c("basic","stud", "perc"),
                                 R = 1999,
                                 se_method = c("agresti", "fisher"),
                                 output = c("htest", "data.frame"),
                                 alternative = c("two.sided", "less", "greater",
                                                 "equivalence", "minimal.effect"),
                                 null.value = NULL,
                                 ...) {
  boot_ci = match.arg(boot_ci)
  ses = match.arg(ses)
  se_method = match.arg(se_method)
  output = match.arg(output)
  alternative = match.arg(alternative)

  # Set default null.value based on effect size type
  if (is.null(null.value)) {
    null.value <- switch(ses,
                         "rb" = 0,
                         "cstat" = 0.5,
                         "odds" = 1,
                         "logodds" = 0)
  }

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

  # Helper function to compute SE for a bootstrap sample
  compute_boot_se <- function(x_boot, y_boot, paired, se_method, mu_val) {
    if (se_method == "agresti") {
      est_results <- ses_compute_agresti(x = x_boot, y = y_boot, paired = paired, mu = mu_val)
      if (is.null(est_results)) {
        return(NA)
      }
      rb <- est_results$rb
      se_rb <- est_results$se_rb
      if (abs(rb) >= 1) {
        return(NA)
      }
      se_z <- se_rb / (1 - rb^2)
      return(se_z)
    } else {
      if (paired || is.null(y_boot)) {
        n1 <- if(is.null(y_boot)) length(x_boot) else length(x_boot)
        n2 <- NULL
      } else {
        n1 <- length(x_boot)
        n2 <- length(y_boot)
      }
      return(se_fisher_z(n1, n2, paired = paired || is.null(y_boot)))
    }
  }

  # Helper function to check for complete separation
  check_complete_separation <- function(p_init) {
    if (p_init <= 0 || p_init >= 1) {
      stop(
        "Complete separation detected (concordance probability = ", p_init, "). ",
        "Bootstrap confidence intervals are not appropriate because resampling ",
        "will produce a degenerate distribution with near-zero variance. ",
        "Use ses_calc() with se_method = 'agresti' for asymptotic inference instead."
      )
    }
  }

  # paired ----
  if(paired == TRUE && !is.null(y)){
    i1 <- x
    i2 <- y

    data <- data.frame(x = i1, y = i2)
    data <- na.omit(data)
    nd = nrow(data)

    if (se_method == "agresti") {
      est_results <- ses_compute_agresti(x = data$x, y = data$y, paired = TRUE, mu = mu)
      raw_rb <- est_results$rb
      raw_SE <- est_results$se_rb / (1 - raw_rb^2)
    } else {
      maxw <- (nd^2 + nd) / 2
      raw_SE = sqrt((2 * nd^3 + 3 * nd^2 + nd) / 6) / maxw
    }

    # Check for complete separation (paired case)
    p_init <- rb_to_cstat(rbs_calc(x = data$y, y = data$x, mu = mu, paired = TRUE))
    check_complete_separation(p_init)

    raw_ses = ses_calc(x = data$x,
                       y = data$y,
                       paired = paired,
                       ses = "rb",
                       mu = mu,
                       alpha = alpha,
                       se_method = se_method,
                       output = "data.frame")

    boots = c()
    boots_se = c()
    for(i in 1:R){
      sampler = sample(1:nrow(data), replace = TRUE)
      res_boot = ses_calc(x = data$x[sampler],
                          y = data$y[sampler],
                          paired = paired,
                          ses = "rb",
                          mu = mu,
                          alpha = alpha,
                          se_method = se_method,
                          output = "data.frame")
      boots = c(boots, atanh(res_boot$estimate))

      rfSE <- compute_boot_se(data$x[sampler], data$y[sampler], paired = TRUE, se_method, mu)
      boots_se = c(boots_se, rfSE)
    }


  } else if(!is.null(y)){
    # two sample -----
    i1 <- na.omit(x)
    i2 <- na.omit(y)
    data <- data.frame(values = c(i1,i2),
                       group = c(rep("x",length(i1)),
                                 rep("y",length(i2))))
    n1 = length(i1)
    n2 = length(i2)

    # Check for complete separation (two-sample case)
    p_init <- mean(sapply(i1, function(xi) mean(i2 < xi) + 0.5 * mean(i2 == xi)))
    check_complete_separation(p_init)

    if (se_method == "agresti") {
      est_results <- ses_compute_agresti(x = i1, y = i2, paired = FALSE, mu = mu)
      raw_rb <- est_results$rb
      raw_SE <- est_results$se_rb / (1 - raw_rb^2)
    } else {
      raw_SE = sqrt((n1 + n2 + 1) / (3 * n1 * n2))
    }

    raw_ses = ses_calc(x = i1,
                       y = i2,
                       paired = paired,
                       ses = "rb",
                       mu = mu,
                       alpha = alpha,
                       se_method = se_method,
                       output = "data.frame")

    boots = c()
    boots_se = c()
    for(i in 1:R){
      sampler = sample(1:nrow(data), replace = TRUE)
      boot_dat = data[sampler,]
      x_boot = subset(boot_dat, group == "x")
      y_boot = subset(boot_dat, group == "y")
      res_boot = ses_calc(x = x_boot$values,
                          y = y_boot$values,
                          paired = paired,
                          ses = "rb",
                          mu = mu,
                          alpha = alpha,
                          se_method = se_method,
                          output = "data.frame")
      boots = c(boots, atanh(res_boot$estimate))

      rfSE <- compute_boot_se(x_boot$values, y_boot$values, paired = FALSE, se_method, mu)
      boots_se = c(boots_se, rfSE)
    }

  } else {
    # one-sample -----
    x1 = na.omit(x)
    n1 = nd =  length(x1)

    # Check for complete separation (one-sample case)
    # One-sample: compare x to mu using signed rank approach
    d <- x1 - mu
    d_nonzero <- d[d != 0]
    if (length(d_nonzero) > 0) {
      # Compute concordance probability from signed ranks
      p_init <- rb_to_cstat(rbs_calc(x = rep(mu, length(x1)), y = x1, mu = 0, paired = TRUE))
      check_complete_separation(p_init)
    }

    if (se_method == "agresti") {
      est_results <- ses_compute_agresti(x = x1, y = NULL, paired = FALSE, mu = mu)
      raw_rb <- est_results$rb
      raw_SE <- est_results$se_rb / (1 - raw_rb^2)
    } else {
      maxw <- (nd^2 + nd) / 2
      raw_SE = sqrt((2 * nd^3 + 3 * nd^2 + nd) / 6) / maxw
    }

    raw_ses = ses_calc(x = x1,
                       paired = paired,
                       ses = "rb",
                       mu = mu,
                       alpha = alpha,
                       se_method = se_method,
                       output = "data.frame")

    boots = c()
    boots_se = c()
    for(i in 1:R){
      sampler = sample(1:length(x1), replace = TRUE)
      x_boot = x1[sampler]

      res_boot = ses_calc(x = x_boot,
                          paired = paired,
                          ses = "rb",
                          mu = mu,
                          alpha = alpha,
                          se_method = se_method,
                          output = "data.frame")
      boots = c(boots, atanh(res_boot$estimate))

      rfSE <- compute_boot_se(x_boot, NULL, paired = FALSE, se_method, mu)
      boots_se = c(boots_se, rfSE)
    }
  }

  if(any(is.infinite(boots))){
    message("Bootstrapped results contain extreme results (i.e., no overlap), caution advised interpreting confidence intervals.")
  }

  # Get CI on Fisher Z scale
  zci = switch(boot_ci,
               "stud" = stud(boots_est = boots, boots_se = boots_se,
                             se0=raw_SE, t0 = atanh(raw_ses$estimate[1L]),
                             alpha = if(alternative %in% c("equivalence", "minimal.effect")) alpha * 2 else alpha),
               "perc" = perc(boots, if(alternative %in% c("equivalence", "minimal.effect")) alpha * 2 else alpha),
               "basic" = basic(boots, t0 = atanh(raw_ses$estimate),
                               if(alternative %in% c("equivalence", "minimal.effect")) alpha * 2 else alpha))

  # Transform back to rb scale
  rci = tanh(zci)
  rboots = tanh(boots)

  # Transform to requested effect size scale
  boots_transformed = switch(ses,
                  "rb" = rboots,
                  "cstat" = rb_to_cstat(rboots),
                  "odds" = rb_to_odds(rboots),
                  "logodds" = log(rb_to_odds(rboots)))

  # Handle infinite values in bootstrap distribution
  if(any(is.infinite(boots_transformed))){
    sum_inf = sum(is.infinite(boots_transformed))
    if(sum_inf/R > .1){
      message("More than 10% of bootstrap estimates contain infinite values, bias and SE calculations may be affected.")
    } else{
      message(paste0("A total of ", sum_inf, " bootstrapped estimates of ", R,
                     " samples are infinite values. Bias and SE estimates are affected; proceed with caution."))
      upper_inf = max(boots_transformed[is.finite(boots_transformed)])
      lower_inf = min(boots_transformed[is.finite(boots_transformed)])
      boots_transformed[boots_transformed == Inf] <- upper_inf
      boots_transformed[boots_transformed == -Inf] <- lower_inf
    }
  }

  # Transform CI to desired estimate scale
  ci = switch(ses,
              "rb" = rci,
              "cstat" = rb_to_cstat(rci),
              "odds" = rb_to_odds(rci),
              "logodds" = log(rb_to_odds(rci)))

  # Transform point estimate to desired scale
  est_val = switch(ses,
                "rb" = raw_ses$estimate,
                "cstat" = rb_to_cstat(raw_ses$estimate),
                "odds" = rb_to_odds(raw_ses$estimate),
                "logodds" = log(rb_to_odds(raw_ses$estimate)))

  ses_name <- switch(ses,
                     "rb" = "Rank-Biserial Correlation",
                     "odds" = "WMW Odds",
                     "logodds" = "WMW Log-Odds",
                     "cstat" = "Concordance")

  # Bootstrap SE on the transformed scale
  boot_se = sd(boots_transformed, na.rm = TRUE)

  # Compute p-value using bootstrap distribution
  # Center the bootstrap distribution at the null value
  # For rb/logodds: null is 0; for cstat: null is 0.5; for odds: null is 1

  if (alternative == "two.sided") {
    # Two-sided test: proportion of bootstrap values at least as extreme as observed
    # Centered at null
    boot_centered <- boots_transformed - null.value
    obs_centered <- est_val - null.value
    boot.pval <- 2 * min(mean(boot_centered <= obs_centered),
                         mean(boot_centered > obs_centered))
  } else if (alternative == "less") {
    # One-sided: proportion of bootstrap values less than or equal to observed
    boot_centered <- boots_transformed - null.value
    obs_centered <- est_val - null.value
    boot.pval <- mean(boot_centered <= obs_centered)
  } else if (alternative == "greater") {
    # One-sided: proportion of bootstrap values greater than or equal to observed
    boot_centered <- boots_transformed - null.value
    obs_centered <- est_val - null.value
    boot.pval <- mean(boot_centered >= obs_centered)
  } else if (alternative == "equivalence") {
    # Equivalence: max of two one-sided p-values
    # Test 1: H0: effect <= low_bound vs H1: effect > low_bound
    # Test 2: H0: effect >= high_bound vs H1: effect < high_bound
    boot_centered_low <- boots_transformed - low_bound
    boot_centered_high <- boots_transformed - high_bound
    obs_centered_low <- est_val - low_bound
    obs_centered_high <- est_val - high_bound

    p_low <- mean(boot_centered_low >= obs_centered_low)  # greater than low bound
    p_high <- mean(boot_centered_high <= obs_centered_high)  # less than high bound
    boot.pval <- max(p_low, p_high)
  } else if (alternative == "minimal.effect") {
    # Minimal effect: min of two one-sided p-values
    # Test 1: H0: effect >= low_bound vs H1: effect < low_bound
    # Test 2: H0: effect <= high_bound vs H1: effect > high_bound
    boot_centered_low <- boots_transformed - low_bound
    boot_centered_high <- boots_transformed - high_bound
    obs_centered_low <- est_val - low_bound
    obs_centered_high <- est_val - high_bound

    p_low <- mean(boot_centered_low <= obs_centered_low)  # less than low bound
    p_high <- mean(boot_centered_high >= obs_centered_high)  # greater than high bound
    boot.pval <- min(p_low, p_high)
  }

  # Compute z-statistic (for reference)
  if (alternative %in% c("equivalence", "minimal.effect")) {
    z_low <- (est_val - low_bound) / boot_se
    z_high <- (est_val - high_bound) / boot_se
    if (alternative == "equivalence") {
      z_stat <- if (abs(z_low) < abs(z_high)) z_low else z_high
    } else {
      z_stat <- if (abs(z_low) < abs(z_high)) z_low else z_high
    }
  } else {
    z_stat <- (est_val - null.value) / boot_se
  }

  # Determine sample type for method string
  if(is.null(y)){
    sample_type = "One Sample"
  } else if(paired == TRUE) {
    sample_type = "Paired Sample"
  } else {
    sample_type = "Two Sample"
  }

  method_desc <- paste0("Bootstrapped ", sample_type, " ", ses_name, " test")

  # Note: bootstrap methodology details
  note_text <- paste0("Bootstrap CI: ", boot_ci,
                      "; SE method: ", if (se_method == "agresti") "Agresti/Lehmann placement" else "Fisher z-transform")

  # Build output
  if (output == "data.frame") {
    effsize = data.frame(
      estimate = est_val,
      SE = boot_se,
      lower.ci = ci[1],
      upper.ci = ci[2],
      conf.level = conf.level,
      boot_ci = boot_ci,
      row.names = ses_name
    )
    return(effsize)

  } else {
    # htest output
    # Set up return structure
    estimate <- est_val
    names(estimate) <- ses_name

    attr(ci, "conf.level") <- conf.level

    names(z_stat) <- "z"

    if (alternative %in% c("equivalence", "minimal.effect")) {
      null_val <- c(low_bound, high_bound)
      names(null_val) <- c("lower bound", "upper bound")
    } else {
      null_val <- null.value
      names(null_val) <- ses_name
    }

    rval = list(
      statistic = z_stat,
      p.value = boot.pval,
      stderr = boot_se,
      conf.int = ci,
      estimate = estimate,
      null.value = null_val,
      alternative = alternative,
      method = method_desc,
      note = note_text,
      boot = boots_transformed,
      data.name = dname,
      call = match.call()
    )

    class(rval) = "htest"

    return(rval)
  }
}

#' @rdname boot_ses_calc
#' @method boot_ses_calc formula
#' @export

boot_ses_calc.formula = function(formula,
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
  y <- do.call("boot_ses_calc", c(DATA, list(...)))
  y$data.name <- DNAME
  y
}
