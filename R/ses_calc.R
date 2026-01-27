#' @title Standardized Effect Size (SES) Calculation
#' @description
#' `r lifecycle::badge('stable')`
#'
#' Calculates non-SMD standardized effect sizes for group comparisons. This function focuses on
#' rank-based and probability-based effect size measures, which are especially useful for
#' non-parametric analyses and when data do not meet normality assumptions.
#'
#' @section Purpose:
#' Use this function when:
#'   - You want to report non-parametric effect size measures
#'   - You need to quantify the magnitude of differences using ranks or probabilities
#'   - Your outcome variable is ordinal
#'   - You want to complement results from Wilcoxon-Mann-Whitney type tests
#'
#' @inheritParams t_TOST
#' @param ses a character string specifying the effect size measure to calculate:
#'     - "rb": rank-biserial correlation (default)
#'     - "odds": Wilcoxon-Mann-Whitney odds
#'     - "logodds": Wilcoxon-Mann-Whitney log-odds
#'     - "cstat": concordance statistic (C-statistic, equivalent to the area under the ROC curve)
#'
#' @param alpha alpha level for confidence interval calculation (default = 0.05).
#' @param mu number indicating the value around which asymmetry (for one-sample or paired samples)
#'   or shift (for independent samples) is to be estimated (default = 0).
#' @param se_method a character string specifying the method for computing standard errors and
#'   confidence intervals:
#'     - "agresti": (default) Uses the Agresti/Lehmann placement-based variance estimation with
#'       confidence intervals computed on the log-odds scale and back-transformed. This method
#'       has better asymptotic properties and faster convergence to normality.
#'     - "fisher": Uses the legacy Fisher z-transformation method for confidence intervals.
#'       This method is retained for backward compatibility.
#' @param output a character string specifying the output format:
#'     - "htest": (default) Returns an object of class "htest" compatible with standard R output.
#'     - "data.frame": Returns a data frame with effect size estimates and confidence intervals.
#' @param null.value a number or vector specifying the null hypothesis value(s) on the scale of the
#'   effect size estimate. Only used when `alternative != "none"`.
#'     - For standard alternatives: a single value (default = 0 for rb/logodds, 0.5 for cstat, 1 for odds)
#'     - For equivalence/minimal.effect: two values representing the lower and upper bounds
#' @param alternative a character string specifying the alternative hypothesis for optional
#'   hypothesis testing:
#'     - "none": (default) No hypothesis test is performed; only effect size and CI are returned.
#'     - "two.sided": Test whether effect differs from null.value
#'     - "less": Test whether effect is less than null.value
#'     - "greater": Test whether effect is greater than null.value
#'     - "equivalence": Test whether effect is between specified bounds
#'     - "minimal.effect": Test whether effect is outside specified bounds
#' @param ... further arguments to be passed to or from methods.
#'
#' @details
#' This function calculates standardized effect sizes that are not standardized mean differences (SMDs).
#' These effect sizes are particularly useful for non-parametric analyses or when data violate
#' assumptions of normality.
#'
#' The available effect size measures are:
#'
#'   - **Rank-biserial correlation ("rb")**: A correlation coefficient based on ranks,
#'     ranging from -1 to 1. It can be interpreted as the difference between the proportion of
#'     favorable pairs and the proportion of unfavorable pairs. For independent samples, this is
#'     equivalent to Cliff's delta.
#'
#'   - **Wilcoxon-Mann-Whitney odds ("odds")**: The ratio of the probability that a
#'     randomly selected observation from group 1 exceeds a randomly selected observation from
#'     group 2, to the probability of the reverse. Values range from 0 to infinity, with 1
#'     indicating no effect.
#'
#'   - **Wilcoxon-Mann-Whitney log-odds ("logodds")**: The natural logarithm of the
#'     WMW odds. This transforms the odds scale to range from negative infinity to positive
#'     infinity, with 0 indicating no effect.
#'
#'   - **Concordance statistic ("cstat")**: The probability that a randomly selected
#'     observation from group 1 exceeds a randomly selected observation from group 2. Also known
#'     as the common language effect size or the area under the ROC curve. Values range from 0 to
#'     1, with 0.5 indicating no effect.
#'
#' ## Standard Error Methods
#'
#' Two methods are available for computing standard errors and confidence intervals:
#'
#'   - **Agresti method** (`se_method = "agresti"`): This method computes the variance of the
#'     concordance probability using the Lehmann/Agresti formula based on placement values.
#'     Confidence intervals are constructed on the log-odds scale (which has the best asymptotic

#'     properties) and then back-transformed to the desired effect size scale. This ensures that
#'     confidence intervals respect the natural bounds of each effect size measure.
#'
#'   - **Fisher method** (`se_method = "fisher"`): This legacy method uses Fisher's z-transformation
#'     (arctanh) for the rank-biserial correlation. Confidence intervals for other effect sizes
#'     are obtained by simple transformation of the rank-biserial CI bounds.
#'
#' The function supports three study designs:
#'   - One-sample design: Compares a single sample to a specified value
#'   - Two-sample independent design: Compares two independent groups
#'   - Paired samples design: Compares paired observations
#'
#' For detailed information on calculation methods, see `vignette("robustTOST")`.
#'
#' @return
#' If `output = "htest"` (default), returns a list with class `"htest"` containing:
#'   - estimate: The effect size estimate
#'   - stderr: Standard error of the estimate
#'   - conf.int: Confidence interval with conf.level attribute
#'   - null.value: The null hypothesis value (if alternative != "none")
#'   - alternative: The alternative hypothesis
#'   - method: Description of the method used
#'   - data.name: Names of the input data
#'   - statistic: Test statistic (only if alternative != "none")
#'   - parameter: Degrees of freedom or other parameters (only if alternative != "none")
#'   - p.value: P-value for the test (only if alternative != "none")
#'
#' If `output = "data.frame"`, returns a data frame containing:
#'   - estimate: The effect size estimate
#'   - SE: Standard error of the estimate
#'   - lower.ci: Lower bound of the confidence interval
#'   - upper.ci: Upper bound of the confidence interval
#'   - conf.level: Confidence level (1-alpha)
#'   - se_method: The SE method used
#'
#' @examples
#' # Example 1: Independent groups comparison (rank-biserial correlation)
#' set.seed(123)
#' group1 <- c(1.2, 2.3, 3.1, 4.6, 5.2, 6.7)
#' group2 <- c(3.5, 4.8, 5.6, 6.9, 7.2, 8.5)
#' ses_calc(x = group1, y = group2, ses = "rb")
#'
#' # Example 2: Using formula notation to calculate WMW odds
#' data(mtcars)
#' ses_calc(formula = mpg ~ am, data = mtcars, ses = "odds")
#'
#' # Example 3: Paired samples with concordance statistic
#' data(sleep)
#' with(sleep, ses_calc(x = extra[group == 1],
#'                      y = extra[group == 2],
#'                      paired = TRUE,
#'                      ses = "cstat"))
#'
#' # Example 4: With hypothesis testing
#' ses_calc(x = group1, y = group2, ses = "rb",
#'          alternative = "two.sided", null.value = 0)
#'
#' # Example 5: Return as data frame (legacy format)
#' ses_calc(x = group1, y = group2, ses = "rb", output = "data.frame")
#'
#' # Example 6: Using Fisher method for backward compatibility
#' ses_calc(x = group1, y = group2, ses = "rb", se_method = "fisher")
#'
#' @references
#' Agresti, A. (1980). Generalized odds ratios for ordinal data. *Biometrics*, 36, 59-67.
#'
#' Bamber, D. (1975). The area above the ordinal dominance graph and the area below the receiver
#' operating characteristic graph. *Journal of Mathematical Psychology*, 12, 387-415.
#'
#' Kerby, D. S. (2014). The simple difference formula: An approach to teaching nonparametric
#' correlation. *Comprehensive Psychology*, 3, 11-IT.
#'
#' Lehmann, E.L. (1975). *Nonparametrics: Statistical Methods Based on Ranks*. Holden-Day.
#'
#' @family effect sizes
#' @name ses_calc
#' @export ses_calc


#ses_calc <- setClass("ses_calc")
ses_calc <- function(x, ...,
                     paired = FALSE,
                     ses = "rb",
                     alpha = 0.05,
                     se_method = c("agresti", "fisher"),
                     output = c("htest", "data.frame"),
                     null.value = NULL,
                     alternative = c("none", "two.sided", "less", "greater",
                                     "equivalence", "minimal.effect")){
  UseMethod("ses_calc")
}

#' @rdname ses_calc
#' @importFrom stats sd cor na.omit setNames wilcox.test terms pnorm qnorm plogis qlogis
#' @method ses_calc default
#' @export

# @method ses_calc default
ses_calc.default = function(x,
                          y = NULL,
                          paired = FALSE,
                          ses = c("rb","odds","logodds","cstat"),
                          alpha = 0.05,
                          mu = 0,
                          se_method = c("agresti", "fisher"),
                          output = c("htest", "data.frame"),
                          null.value = NULL,
                          alternative = c("none", "two.sided", "less", "greater",
                                          "equivalence", "minimal.effect"),
                          ...) {

  ses = match.arg(ses)
  se_method = match.arg(se_method)
  output = match.arg(output)
  alternative = match.arg(alternative)

  # Set default null.value based on effect size type (consistent with boot_ses_calc)
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

  if(is.null(y)){
    sample_type = "One Sample"
  } else if(paired == TRUE) {
    sample_type = "Paired Sample"
  } else {
    sample_type = "Two Sample"
  }

  if (!is.null(y)) {
    dname <- paste(deparse(substitute(x)), "and",
                   deparse(substitute(y)))
  } else {
    dname <- deparse(substitute(x))
  }

  if(!is.numeric(alpha) || alpha <=0 || alpha >=1){
    stop("The alpha must be a numeric value between 0 and 1")
  }

  # Handle NA removal
  if (is.null(y)) {
    x <- na.omit(x)
    n1 <- length(x)
    n2 <- NULL
  } else if (paired) {
    data <- data.frame(x = x, y = y)
    data <- na.omit(data)
    x <- data$x
    y <- data$y
    n1 <- nrow(data)
    n2 <- NULL
  } else {
    x <- na.omit(x)
    y <- na.omit(y)
    n1 <- length(x)
    n2 <- length(y)
  }

  # Compute point estimate using existing function
  # Note: rbs_calc for paired/one-sample expects x and y swapped
  # For one-sample: rbs() sets y=rep(0,n), paired=TRUE, then swaps to x=zeros, y=original
  # For paired: we pass y as first arg, x as second (see rbs() lines 92-100)
  if (is.null(y)) {
    # One-sample: match rbs() convention - zeros first, then data
    r_rbs <- rbs_calc(x = rep(0, length(x)), y = x, mu = mu, paired = TRUE)
  } else if (paired) {
    # Paired: match rbs() convention - swap x and y
    r_rbs <- rbs_calc(x = y, y = x, mu = mu, paired = TRUE)
  } else {
    # Two-sample: no swap
    r_rbs <- rbs_calc(x = x, y = y, mu = mu, paired = FALSE)
  }

  # Get concordance and other transformations
  p_hat <- rb_to_cstat(r_rbs)
  alpha_hat <- rb_to_odds(r_rbs)
  eta_hat <- log(alpha_hat)

  # Compute SE and CI based on method
  if (se_method == "agresti") {
    # New Agresti/Lehmann method
    est_results <- ses_compute_agresti(x = x, y = y, paired = paired, mu = mu)

    if (is.null(est_results)) {
      stop("Unable to compute effect size - check that data has sufficient non-zero differences")
    }

    ci_results <- ses_ci_logodds(est_results, conf.level = conf.level)

    # Extract results for requested effect size
    se_val <- switch(ses,
                     "rb" = est_results$se_rb,
                     "cstat" = est_results$se_cstat,
                     "odds" = est_results$se_odds,
                     "logodds" = est_results$se_logodds)

    ci_val <- switch(ses,
                     "rb" = ci_results$ci_rb,
                     "cstat" = ci_results$ci_cstat,
                     "odds" = ci_results$ci_odds,
                     "logodds" = ci_results$ci_logodds)

  } else {
    # Legacy Fisher method
    ci_results_fisher <- ses_ci_fisher(r_rbs, n1, n2, paired = paired || is.null(y),
                                        conf.level = conf.level)

    se_val <- switch(ses,
                     "rb" = ci_results_fisher$se_rb,
                     "cstat" = ci_results_fisher$se_cstat,
                     "odds" = ci_results_fisher$se_odds,
                     "logodds" = ci_results_fisher$se_logodds)

    ci_val <- switch(ses,
                     "rb" = ci_results_fisher$ci_rb,
                     "cstat" = ci_results_fisher$ci_cstat,
                     "odds" = ci_results_fisher$ci_odds,
                     "logodds" = ci_results_fisher$ci_logodds)
  }

  # Get point estimate for requested effect size
  est_val <- switch(ses,
                    "rb" = r_rbs,
                    "cstat" = p_hat,
                    "odds" = alpha_hat,
                    "logodds" = eta_hat)

  ses_name <- switch(ses,
                     "rb" = "Rank-Biserial Correlation",
                     "odds" = "WMW Odds",
                     "logodds" = "WMW Log-Odds",
                     "cstat" = "Concordance")

  # Build output
  if (output == "data.frame") {
    effsize = data.frame(
      estimate = est_val,
      SE = se_val,
      lower.ci = ci_val[1],
      upper.ci = ci_val[2],
      conf.level = conf.level,
      se_method = se_method,
      row.names = ses_name
    )
    return(effsize)

  } else {
    # htest output
    method_desc <- paste0(ses_name, " (", sample_type, ")")
    if (se_method == "agresti") {
      method_desc <- paste0(method_desc, "\nSE: Agresti/Lehmann; CI: log-odds transform")
    } else {
      method_desc <- paste0(method_desc, "\nSE: Fisher z-transform")
    }

    # Set up estimate with name
    estimate <- est_val
    names(estimate) <- ses_name

    # Set up confidence interval
    conf.int <- ci_val
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

    # Add hypothesis test components if requested
    if (alternative != "none") {
      # Compute z-statistic and p-value
      if (alternative %in% c("equivalence", "minimal.effect")) {
        z_low <- (est_val - low_bound) / se_val
        z_high <- (est_val - high_bound) / se_val

        p_low <- pnorm(z_low, lower.tail = FALSE)  # H0: effect <= low_bound
        p_high <- pnorm(z_high)  # H0: effect >= high_bound

        if (alternative == "equivalence") {
          # TOST: max of p-values, report z closest to null
          p_val <- max(p_low, p_high)
          z_stat <- if (abs(z_low) < abs(z_high)) z_low else z_high
        } else {
          # minimal.effect: min of p-values
          p_val <- min(p_low, p_high)
          z_stat <- if (abs(z_low) < abs(z_high)) z_low else z_high
        }

        null_val <- c(low_bound, high_bound)
        names(null_val) <- c("lower bound", "upper bound")
      } else {
        # Standard alternatives
        z_stat <- (est_val - null.value) / se_val

        p_val <- switch(alternative,
                        "two.sided" = 2 * pnorm(-abs(z_stat)),
                        "less" = pnorm(z_stat),
                        "greater" = pnorm(z_stat, lower.tail = FALSE))

        null_val <- null.value
        names(null_val) <- ses_name
      }

      names(z_stat) <- "z"
      rval$statistic <- z_stat
      rval$p.value <- p_val
      rval$null.value <- null_val
    }

    class(rval) <- "htest"
    return(rval)
  }
}

#' @rdname ses_calc
#' @method ses_calc formula
#' @export

ses_calc.formula = function(formula,
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
  y <- do.call("ses_calc", c(DATA, list(...)))


  # Update data.name for htest output

  if (inherits(y, "htest")) {
    y$data.name <- DNAME
  }

  y

}
