#' @title Standardized Mean Difference (SMD) Calculation
#' @description
#' `r lifecycle::badge('stable')`
#'
#' Calculates standardized mean difference (SMD) effect sizes and their confidence intervals
#' from raw data, with optional hypothesis testing.
#'
#' @section Purpose:
#' Use this function when:
#'   * You need to calculate standardized effect sizes (Cohen's d, Hedges' g, Glass's delta)
#'   * You want confidence intervals for your effect size estimates
#'   * You need effect sizes for meta-analysis or reporting
#'   * You want to compare effect sizes across different studies or measures
#'   * You want to test hypotheses about effect size magnitudes (e.g., equivalence testing)
#'
#' @inheritParams t_TOST
#' @inheritParams boot_t_TOST
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
#' @param test_method a character string specifying the method for hypothesis testing:
#'     - "z": Use z-statistic (normal distribution)
#'     - "t": Use t-statistic with degrees of freedom from the SMD calculation
#'
#' @details
#' This function calculates standardized mean differences (SMD) for various study designs:
#'
#'   * One-sample design: Standardizes the difference between the sample mean and zero (or other specified value)
#'   * Two-sample independent design: Standardizes the difference between two group means
#'   * Paired samples design: Standardizes the mean difference between paired observations
#'
#' The function supports multiple SMD variants:
#'   * Cohen's d: Classic standardized mean difference (bias_correction = FALSE)
#'   * Hedges' g: Bias-corrected version of Cohen's d (bias_correction = TRUE)
#'   * Glass's delta: Uses only one group's standard deviation as the denominator (glass = "glass1" or "glass2")
#'   * Repeated measures d: Accounts for correlation in paired designs (rm_correction = TRUE)
#'
#' Different confidence interval calculation methods are available:
#'   * "nct": Uses the noncentral t-distribution (most accurate in most cases)
#'   * "goulet": Uses the Goulet-Pelletier method
#'   * "t": Uses the central t-distribution
#'   * "z": Uses the normal distribution
#'
#' For detailed information on calculation methods, see `vignette("SMD_calcs")`.
#'
#' @return
#' If `output = "htest"` (default), returns a list with class `"htest"` containing:
#'   - estimate: The SMD estimate (Cohen's d, Hedges' g, or Glass's delta)
#'   - stderr: Standard error of the estimate
#'   - conf.int: Confidence interval with conf.level attribute
#'   - alternative: A character string describing the alternative hypothesis
#'   - method: A character string indicating what type of test was performed
#'   - data.name: A character string giving the name(s) of the data
#'   - statistic: Test statistic (only if alternative != "none")
#'   - parameter: Degrees of freedom (only if test_method = "t" and alternative != "none")
#'   - p.value: P-value for the test (only if alternative != "none")
#'   - null.value: The specified hypothesized value(s) (only if alternative != "none")
#'
#' If `output = "data.frame"`, returns a data frame containing:
#'   - estimate: The SMD estimate
#'   - SE: Standard error of the estimate
#'   - lower.ci: Lower bound of the confidence interval
#'   - upper.ci: Upper bound of the confidence interval
#'   - conf.level: Confidence level (1-alpha)
#'
#' @examples
#' # Example 1: Independent groups comparison (Cohen's d)
#' set.seed(123)
#' group1 <- rnorm(30, mean = 100, sd = 15)
#' group2 <- rnorm(30, mean = 110, sd = 18)
#' smd_calc(x = group1, y = group2, bias_correction = FALSE)
#'
#' # Example 2: Independent groups with formula notation (Hedges' g)
#' df <- data.frame(
#'   value = c(group1, group2),
#'   group = factor(rep(c("A", "B"), each = 30))
#' )
#' smd_calc(formula = value ~ group, data = df)
#'
#' # Example 3: Paired samples with repeated measures correction
#' before <- c(5.1, 4.8, 6.2, 5.7, 6.0, 5.5, 4.9, 5.8)
#' after <- c(5.6, 5.2, 6.7, 6.1, 6.5, 5.8, 5.3, 6.2)
#' smd_calc(x = before, y = after, paired = TRUE, rm_correction = TRUE)
#'
#' # Example 4: Glass's delta (using only first group's SD)
#' smd_calc(x = group1, y = group2, glass = "glass1")
#'
#' # Example 5: Two-sided test against null of 0
#' smd_calc(x = group1, y = group2,
#'          alternative = "two.sided", null.value = 0)
#'
#' # Example 6: Equivalence test (TOST)
#' smd_calc(x = group1, y = group2,
#'          alternative = "equivalence", null.value = c(-0.5, 0.5))
#'
#' # Example 7: Using t-distribution for test
#' smd_calc(x = group1, y = group2,
#'          alternative = "two.sided", null.value = 0,
#'          test_method = "t", smd_ci = "t")
#'
#' # Example 8: Legacy data.frame output
#' smd_calc(x = group1, y = group2, output = "data.frame")
#'
#' @family effect sizes
#' @name smd_calc
#' @export smd_calc

#smd_calc <- setClass("smd_calc")
smd_calc <- function(x, ...,
                     paired = FALSE,
                     var.equal = FALSE,
                     alpha = 0.05,
                     bias_correction = TRUE,
                     rm_correction = FALSE,
                     glass = NULL,
                     smd_ci = c("nct", "goulet", "t", "z"),
                     output = c("htest", "data.frame"),
                     null.value = 0,
                     alternative = c("none", "two.sided", "less", "greater",
                                     "equivalence", "minimal.effect"),
                     test_method = c("z", "t")){
  UseMethod("smd_calc")
}

#' @rdname smd_calc
#' @importFrom stats sd cor na.omit setNames t.test terms nlm optim optimize pt pnorm
#' @method smd_calc default
#' @export

# @method smd_calc default
smd_calc.default = function(x,
                            y = NULL,
                            paired = FALSE,
                            var.equal = FALSE,
                            alpha = 0.05,
                            mu = 0,
                            bias_correction = TRUE,
                            rm_correction = FALSE,
                            glass = NULL,
                            smd_ci = c("nct", "goulet", "t", "z"),
                            output = c("htest", "data.frame"),
                            null.value = 0,
                            alternative = c("none", "two.sided", "less", "greater",
                                            "equivalence", "minimal.effect"),
                            test_method = c("z", "t"),
                            ...) {

  if(is.null(glass)){
    glass = "no"
  }
  smd_ci = match.arg(smd_ci)
  output = match.arg(output)
  alternative = match.arg(alternative)
  test_method = match.arg(test_method)

  if(bias_correction){
    smd_type = 'g'
  } else {
    smd_type = 'd'
  }

  if(is.null(y)){
    sample_type = "One Sample"
  } else if(paired == TRUE) {
    sample_type = "Paired Sample"
  } else {
    sample_type = "Two Sample"
  }

  if(glass == "glass1" || glass == "glass2"){
    if(glass == "glass1"){
      denom = "glass1"
    }

    if(glass == "glass2"){
      denom = "glass2"
    }
  } else{
    if(sample_type != "Two Sample" ){
      if(rm_correction){
        denom = "rm"
      } else {
        denom = "z"
      }
    } else{
      denom = "d"
    }
  }

  if(!is.numeric(alpha) || alpha <=0 || alpha >=1){
    stop("The alpha must be a numeric value between 0 and 1")
  }

  if (!is.null(y)) {
    dname <- paste(deparse(substitute(x)), "and",
                   deparse(substitute(y)))
  }
  else {
    dname <- deparse(substitute(x))
  }

  # Handle equivalence/minimal.effect bounds
  if (alternative %in% c("equivalence", "minimal.effect")) {
    if (length(null.value) != 2) {
      stop("For equivalence or minimal.effect testing, null.value must be a vector of two values")
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

  # Warn if test_method and smd_ci mismatch
  if (alternative != "none" && test_method != smd_ci && smd_ci %in% c("nct", "goulet", "t", "z")) {
    warning("test_method ('", test_method, "') differs from smd_ci ('", smd_ci,
            "'). Consider aligning these for consistency.")
  }

  if(paired == TRUE && !missing(y)){
    i1 <- x
    i2 <- y
    data <- data.frame(i1 = i1, i2 = i2)
    data <- na.omit(data)
    colnames(data) = c("i1", "i2")
    data2 =  data
    data2$diff = data2$i2 - data2$i1

    n <- nrow(data)
    i1 <- data$i1
    i2 <- data$i2
    m1 <- mean(i1)
    m2 <- mean(i2) + mu
    sd1  <- sd(i1)
    sd2  <- sd(i2)
    r12 <- cor(i1, i2)

    # Calculate Cohens d
    cohen_res = d_est_pair(
      n = n,
      m1 = m1,
      m2 = m2,
      sd1 = sd1,
      sd2 = sd2,
      r12 = r12,
      type = smd_type,
      denom = denom,
      alpha = alpha/2,
      smd_ci = smd_ci
    )

  } else if(!missing(y)){

    x1 = na.omit(x)
    y1 = na.omit(y)
    n1 = length(x1)
    n2 = length(y1)
    m1 = mean(x1)
    m2 = mean(y1) + mu
    sd1 = sd(x1)
    sd2 = sd(y1)

    cohen_res = d_est_ind(
      n1 = n1,
      n2 = n2,
      m1 = m1,
      m2 = m2,
      sd1 = sd1,
      sd2 = sd2,
      type = smd_type,
      var.equal = var.equal,
      alpha = alpha/2,
      denom = denom,
      smd_ci = smd_ci
    )

  } else {

    x1 = na.omit(x)
    n1 = length(x1)
    m1 = mean(x1) + mu
    sd1 = sd(x1)

    cohen_res = d_est_one(
      n = n1,
      mu = m1,
      sd = sd1,
      type = smd_type,
      testValue = 0,
      alpha = alpha/2,
      smd_ci = smd_ci
    )

  }

  # Legacy data.frame output
  if (output == "data.frame") {
    effsize = data.frame(
      estimate = c(cohen_res$d),
      SE = c(cohen_res$d_sigma),
      lower.ci = c(cohen_res$dlow),
      upper.ci = c( cohen_res$dhigh),
      conf.level = c((1-alpha)),
      row.names = c(cohen_res$smd_label)
    )
    return(effsize)
  }

  # htest output
  est_val <- cohen_res$d
  se_val <- cohen_res$d_sigma
  df <- cohen_res$d_df

  # Determine SMD label
  smd_label <- if (bias_correction) {
    if (glass != "no" && glass %in% c("glass1", "glass2")) {
      paste0("Glass's ", ifelse(glass == "glass1", "delta1", "delta2"))
    } else {
      "Hedges' g"
    }
  } else {
    "Cohen's d"
  }

  # Construct method description
  method_suffix <- if (alternative != "none") "test" else "estimate with CI"
  method_desc <- paste0(sample_type, " ", smd_label, " ", method_suffix)

  # Set up estimate with name
  estimate <- est_val
  names(estimate) <- smd_label

  # Recalculate CI for equivalence (90% CI when alpha=0.05)
  if (alternative %in% c("equivalence", "minimal.effect")) {
    # Need to recalculate CI with conf.level = 1 - 2*alpha
    if(paired == TRUE && !missing(y)){
      cohen_res_ci = d_est_pair(
        n = n, m1 = m1, m2 = m2, sd1 = sd1, sd2 = sd2, r12 = r12,
        type = smd_type, denom = denom, alpha = alpha, smd_ci = smd_ci
      )
    } else if(!missing(y)){
      cohen_res_ci = d_est_ind(
        n1 = n1, n2 = n2, m1 = m1, m2 = m2, sd1 = sd1, sd2 = sd2,
        type = smd_type, var.equal = var.equal, alpha = alpha, denom = denom, smd_ci = smd_ci
      )
    } else {
      cohen_res_ci = d_est_one(
        n = n1, mu = m1, sd = sd1, type = smd_type, testValue = 0,
        alpha = alpha, smd_ci = smd_ci
      )
    }
    conf.int <- c(cohen_res_ci$dlow, cohen_res_ci$dhigh)
  } else {
    conf.int <- c(cohen_res$dlow, cohen_res$dhigh)
  }
  attr(conf.int, "conf.level") <- conf.level

  # Build basic htest structure
  rval <- list(
    estimate = estimate,
    stderr = se_val,
    conf.int = conf.int,
    alternative = alternative,
    method = method_desc,
    data.name = dname
  )

  # Add hypothesis test components only if alternative != "none"
  if (alternative != "none") {
    if (alternative %in% c("equivalence", "minimal.effect")) {
      # Two one-sided tests
      stat_low <- (est_val - low_bound) / se_val
      stat_high <- (est_val - high_bound) / se_val

      if (test_method == "t") {
        p_low <- pt(stat_low, df = df, lower.tail = FALSE)
        p_high <- pt(stat_high, df = df, lower.tail = TRUE)
      } else {
        p_low <- pnorm(stat_low, lower.tail = FALSE)
        p_high <- pnorm(stat_high, lower.tail = TRUE)
      }

      if (alternative == "equivalence") {
        p_val <- max(p_low, p_high)
        test_stat <- if (abs(stat_low) < abs(stat_high)) stat_low else stat_high
      } else {  # minimal.effect
        p_val <- min(p_low, p_high)
        test_stat <- if (abs(stat_low) < abs(stat_high)) stat_low else stat_high
      }

      null_val <- c(low_bound, high_bound)
      names(null_val) <- c("lower bound", "upper bound")

    } else {
      # Standard alternatives
      test_stat <- (est_val - null.value) / se_val

      if (test_method == "t") {
        p_val <- switch(alternative,
                        "two.sided" = 2 * pt(-abs(test_stat), df = df),
                        "less" = pt(test_stat, df = df),
                        "greater" = pt(test_stat, df = df, lower.tail = FALSE))
      } else {
        p_val <- switch(alternative,
                        "two.sided" = 2 * pnorm(-abs(test_stat)),
                        "less" = pnorm(test_stat),
                        "greater" = pnorm(test_stat, lower.tail = FALSE))
      }

      null_val <- null.value
      names(null_val) <- smd_label
    }

    names(test_stat) <- test_method
    rval$statistic <- test_stat
    if (test_method == "t") {
      df_named <- df
      names(df_named) <- "df"
      rval$parameter <- df_named
    }
    rval$p.value <- p_val
    rval$null.value <- null_val
  }

  class(rval) <- "htest"
  return(rval)

}

#' @rdname smd_calc
#' @method smd_calc formula
#' @export

smd_calc.formula = function(formula,
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
  y <- do.call("smd_calc", c(DATA, list(...)))

  # Update data.name for htest output
  if (inherits(y, "htest")) {
    y$data.name <- DNAME
  }

  y

}
