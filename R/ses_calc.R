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
#'     - "auto": (default) Automatically selects the most appropriate method based on the
#'       study design. Resolves to "score" for two-sample independent designs and "agresti"
#'       for one-sample or paired designs.
#'     - "agresti": Uses the Agresti/Lehmann placement-based variance estimation with
#'       confidence intervals computed on the log-odds scale and back-transformed. This method
#'       has better asymptotic properties and faster convergence to normality. Available for
#'       all designs (one-sample, paired, and two-sample independent).
#'     - "score": Uses the Fay-Malinovsky score-type approach based on the
#'       proportional odds model. Confidence intervals are constructed by
#'       test inversion (finding the values of the concordance probability
#'       where the score statistic equals the critical value), then
#'       transformed to the requested effect size scale. This method has
#'       better small-sample coverage than the Wald-based methods and
#'       produces confidence intervals that are compatible with the
#'       Wilcoxon-Mann-Whitney test. **Only available for two-sample
#'       independent designs.** See Fay and Malinovsky (2018) for details.
#'     - "fisher": Uses the legacy Fisher z-transformation method for confidence intervals.
#'       This method is retained for backward compatibility.
#' @param correct logical; whether to apply a continuity correction to the
#'   test statistic and confidence interval. When `se_method = "score"`,
#'   setting `correct = TRUE` produces p-values that match
#'   `wilcox.test(..., exact = FALSE, correct = TRUE)`.
#'   Default is `FALSE`. Only used with `se_method = "score"`.
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
#' These effect sizes are particularly useful to compliment non-parametric analyses or when analyzing ordinal data.
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
#' Three methods are available for computing standard errors and confidence intervals:
#'
#'   - **Agresti method** (`se_method = "agresti"`): This method computes the variance of the
#'     concordance probability \eqn{\hat{p} = \Pr(X > Y)} using the Lehmann/Agresti placement-based
#'     formula. For two independent samples with sizes \eqn{n_1} and \eqn{n_2}, the placement values
#'     \eqn{V_i = \Pr(X > Y | X = x_i)} and \eqn{W_j = \Pr(X > Y | Y = y_j)} are computed, and the
#'     variance is estimated as:
#'
#'     \deqn{\widehat{\mathrm{Var}}(\hat{p}) = \frac{\bar{V^2} - \hat{p}^2}{n_1}
#'     + \frac{\bar{W^2} - \hat{p}^2}{n_2}}
#'
#'     where \eqn{\bar{V^2} = \frac{1}{n_1}\sum V_i^2} and \eqn{\bar{W^2} = \frac{1}{n_2}\sum W_j^2}.
#'     For paired samples, the variance is derived from the Wilcoxon signed-rank statistic using
#'     its null-distribution variance with a tie correction.
#'
#'     Standard errors for other effect size scales are obtained via the delta method. Let
#'     \eqn{\mathrm{SE}_p} denote the standard error of \eqn{\hat{p}}. Then:
#'       - \eqn{\mathrm{SE}_{\mathrm{rb}} = 2 \cdot \mathrm{SE}_p} (since \eqn{r_b = 2p - 1})
#'       - \eqn{\mathrm{SE}_{\eta} = \mathrm{SE}_p / [\hat{p}(1 - \hat{p})]} for the log-odds
#'         \eqn{\eta = \log[\hat{p}/(1 - \hat{p})]}
#'       - \eqn{\mathrm{SE}_{\alpha} = \mathrm{SE}_p / (1 - \hat{p})^2} for the odds
#'         \eqn{\alpha = \hat{p}/(1 - \hat{p})}
#'
#'     Confidence intervals are constructed on the log-odds scale and back-transformed to the
#'     requested effect size scale. This ensures that intervals respect the natural bounds of each
#'     measure (e.g., \eqn{[0, 1]} for cstat, \eqn{[-1, 1]} for rb).
#'
#'   - **Score method** (`se_method = "score"`): Uses the Fay-Malinovsky approach based on the
#'     V_LAPH variance function from the proportional odds model. Confidence intervals are
#'     constructed via test inversion: finding the values of \eqn{\phi} (concordance probability)
#'     where the score statistic equals the critical value. This method has better small-sample
#'     coverage than Wald-type methods because the variance is evaluated at the candidate
#'     parameter value, not the estimate.
#'
#'     The reported standard error is descriptive (computed from V_LAPH at \eqn{\hat{\phi}}),
#'     but the confidence interval is **not** computed as estimate ± z × SE. Instead, it comes
#'     from test inversion, which gives better coverage properties.
#'
#'     When `correct = TRUE`, a continuity correction is applied that makes p-values match
#'     `wilcox.test(..., exact = FALSE, correct = TRUE)`. This method is only available for
#'     two-sample independent designs.
#'
#'   - **Fisher method** (`se_method = "fisher"`): This legacy method uses Fisher's z-transformation
#'     (arctanh) for the rank-biserial correlation. Confidence intervals for other effect sizes
#'     are obtained by simple transformation of the rank-biserial CI bounds.
#'
#' ## Boundary Case Handling (Complete Separation)
#'
#' When there is complete separation between groups (i.e., all pairwise comparisons favor
#' one group), the concordance probability \eqn{\hat{p}} equals exactly 0 or 1.
#' This leads to undefined odds and log-odds, and a degenerate (zero) placement-based
#' variance.
#'
#' For `se_method = "agresti"`, a Haldane-type shrinkage correction is applied:
#' \deqn{\tilde{p} = \frac{C + 0.5}{N_{\mathrm{pairs}} + 1}}
#' where \eqn{C} is the concordance count and \eqn{N_{\mathrm{pairs}}} is the total
#' number of pairwise comparisons (\eqn{n_1 n_2} for two-sample designs, or
#' \eqn{N(N+1)/2} for paired/one-sample designs where \eqn{N} is the number of
#' non-zero differences). This corresponds to the posterior mean under a Jeffreys
#' Beta(0.5, 0.5) prior and shrinks the estimate toward 0.5, with stronger shrinkage
#' for smaller samples.
#'
#' The Agresti placement variance is then evaluated at the corrected estimate. A message
#' is printed when this correction is applied. For more reliable inference at boundaries,
#' consider [perm_ses_test()] for permutation-based p-values and intervals, or
#' `se_method = "score"` (two-sample designs only) for score-type intervals that handle
#' boundaries naturally without correction.
#'
#' For `se_method = "score"` (two-sample only), no correction is needed. The score-type
#' CI is constructed via test inversion, where the variance function V_LAPH(phi) is
#' evaluated at candidate parameter values in the interior of (0, 1). When
#' \eqn{\hat{p} = 1}, the upper CI bound is trivially 1 and the lower bound is found
#' by root-finding. No modification of the point estimate is required.
#'
#' For two-sample designs using the Agresti method, a score-type CI is automatically
#' used as a fallback at boundaries for better interval coverage.
#'
#' ## Hypothesis Testing
#'
#' When `alternative != "none"`, a hypothesis test is performed. The approach depends on the
#' SE method:
#'
#'   - **Agresti method**: All hypothesis tests are conducted on the **log-odds scale**, regardless
#'     of which effect size is requested. The log-odds scale \eqn{\eta = \log[p / (1-p)]} has
#'     superior asymptotic properties: it is unbounded (\eqn{-\infty} to \eqn{+\infty}), converges
#'     more quickly to normality, and the null hypothesis (\eqn{\eta_0 = 0}) is interior to the
#'     parameter space.
#'
#'     The user specifies `null.value` on the scale of their chosen effect size. These values are
#'     internally transformed to the log-odds scale for testing:
#'       - rb: \eqn{\eta_0 = \log[(r_0 + 1) / (1 - r_0)]}
#'       - cstat: \eqn{\eta_0 = \log[p_0 / (1 - p_0)]}
#'       - odds: \eqn{\eta_0 = \log(\alpha_0)}
#'       - logodds: \eqn{\eta_0 = \eta_0} (no transformation needed)
#'
#'     The z-statistic is computed as \eqn{z = (\hat{\eta} - \eta_0) / \mathrm{SE}_\eta}, and the
#'     p-value is obtained from the standard normal distribution. When the user specifies custom
#'     null values (not the default), a message reports the transformation to the log-odds scale.
#'
#'   - **Fisher method**: Hypothesis tests are conducted on the scale of the requested effect size,
#'     using the Wald z-statistic \eqn{z = (\hat{\theta} - \theta_0) / \mathrm{SE}_\theta}.
#'
#' For equivalence testing (`alternative = "equivalence"`), the TOST procedure is used: two
#' one-sided tests are performed against the lower and upper bounds, and the p-value is the
#' maximum of the two one-sided p-values. For minimal effect testing
#' (`alternative = "minimal.effect"`), the p-value is the minimum.
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
#' # Example 7: Using score method for WMW-compatible CIs (two-sample only)
#' ses_calc(x = group1, y = group2, ses = "cstat", se_method = "score")
#'
#' # Example 8: Score method with continuity correction (matches wilcox.test)
#' ses_calc(x = group1, y = group2, ses = "cstat", se_method = "score",
#'          correct = TRUE, alternative = "two.sided", null.value = 0.5)
#'
#' @references
#' Agresti, A. (1980). Generalized odds ratios for ordinal data. *Biometrics*, 36, 59-67.
#'
#' Bamber, D. (1975). The area above the ordinal dominance graph and the area below the receiver
#' operating characteristic graph. *Journal of Mathematical Psychology*, 12, 387-415.
#'
#' Fay, M.P. and Malinovsky, Y. (2018). Confidence Intervals of the Mann-Whitney Parameter
#' that are Compatible with the Wilcoxon-Mann-Whitney Test. *Statistics in Medicine*,
#' 37, 3991-4006. \doi{10.1002/sim.7890}
#'
#' Kerby, D. S. (2014). The simple difference formula: An approach to teaching nonparametric
#' correlation. *Comprehensive Psychology*, 3, 11-IT.
#'
#' Lehmann, E.L. (1975). *Nonparametrics: Statistical Methods Based on Ranks*. Holden-Day.
#'
#' O'Brien, R.G. & Castelloe, J. (2006). Exploiting the link between the Wilcoxon-Mann-Whitney
#' test and a simple odds statistic. *SUGI 31 Proceedings*, Paper 209-31.
#'
#' @family effect sizes
#' @name ses_calc
#' @export ses_calc


#ses_calc <- setClass("ses_calc")
ses_calc <- function(x, ...,
                     paired = FALSE,
                     ses = "rb",
                     alpha = 0.05,
                     se_method = c("auto", "agresti", "score", "fisher"),
                     correct = FALSE,
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
                          se_method = c("auto", "agresti", "score", "fisher"),
                          correct = FALSE,
                          output = c("htest", "data.frame"),
                          null.value = NULL,
                          alternative = c("none", "two.sided", "less", "greater",
                                          "equivalence", "minimal.effect"),
                          ...) {

  ses = match.arg(ses)
  se_method = match.arg(se_method)
  output = match.arg(output)
  alternative = match.arg(alternative)

  if(!is.numeric(alpha) || alpha <=0 || alpha >=1){
    stop("alpha must be a numeric value between 0 and 1")
  }

  # Determine design type
  is_two_sample <- !is.null(y) && !paired

  # Resolve "auto" se_method based on design type
  if (se_method == "auto") {
    se_method <- if (is_two_sample) "score" else "agresti"
  }

  # Score method only available for two-sample independent
  if (se_method == "score" && !is_two_sample) {
    stop("se_method = 'score' is only available for two-sample independent designs. ",
         "For one-sample or paired designs, use se_method = 'agresti'.")
  }

  # Continuity correction only used with score method
  if (correct && se_method != "score") {
    message("Continuity correction (correct = TRUE) is only used with se_method = 'score'. Ignoring.")
  }

  # Track whether user provided null.value (for messaging about log-odds transformation)
  null.value_original <- null.value

  # Set default null.value based on effect size type (consistent with boot_ses_calc)
  if (is.null(null.value)) {
    null.value <- switch(ses,
                         "rb" = 0,
                         "cstat" = 0.5,
                         "odds" = 1,
                         "logodds" = 0)
  }

  # Determine if user specified custom bounds (for messaging)
  user_specified_bounds <- !is.null(null.value_original)

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
    # Agresti/Lehmann method with Haldane boundary correction
    est_results <- ses_compute_agresti(x = x, y = y, paired = paired, mu = mu,
                                        use_score_fallback = TRUE,
                                        conf.level = conf.level)

    if (is.null(est_results)) {
      stop("Unable to compute effect size - check that data has sufficient non-zero differences")
    }

    # Message if boundary correction was applied
    if (est_results$boundary_corrected) {
      msg <- paste0(
        "Complete separation detected (all pairwise comparisons favor one group). ",
        "A Haldane-type shrinkage correction was applied to enable confidence ",
        "interval construction on the log-odds scale."
      )

      # If two-sample and score CI fallback was used, note that
      if (est_results$boundary_used_score_ci) {
        msg <- paste0(msg, " Score-type CIs used for better boundary behavior.")
      } else if (!paired && !is.null(y)) {
        msg <- paste0(msg, " For more reliable inference at boundaries, consider ",
                      "se_method = 'score' for score-type intervals.")
      } else {
        msg <- paste0(msg, " For more reliable inference at boundaries, consider ",
                      "perm_ses_test() for permutation-based p-values and intervals.")
      }

      message(msg)
    }

    # If score CI fallback was used, use those CIs; otherwise use standard log-odds CIs
    if (est_results$boundary_used_score_ci && !is.null(est_results$score_ci_cstat)) {
      # Transform score CI bounds to all scales
      ci_cstat <- est_results$score_ci_cstat
      ci_rb <- 2 * ci_cstat - 1
      ci_odds <- ci_cstat / (1 - ci_cstat)
      ci_logodds <- log(ci_cstat / (1 - ci_cstat))

      ci_val <- switch(ses,
                       "rb" = ci_rb,
                       "cstat" = ci_cstat,
                       "odds" = ci_odds,
                       "logodds" = ci_logodds)
    } else {
      ci_results <- ses_ci_logodds(est_results, conf.level = conf.level)

      ci_val <- switch(ses,
                       "rb" = ci_results$ci_rb,
                       "cstat" = ci_results$ci_cstat,
                       "odds" = ci_results$ci_odds,
                       "logodds" = ci_results$ci_logodds)
    }

    # Extract SE for requested effect size
    se_val <- switch(ses,
                     "rb" = est_results$se_rb,
                     "cstat" = est_results$se_cstat,
                     "odds" = est_results$se_odds,
                     "logodds" = est_results$se_logodds)

  } else if (se_method == "score") {
    # Fay-Malinovsky score-type method (two-sample only)
    # Note: validation already ensured this is two-sample independent

    # Compute tie factor
    tf <- wmw_tie_factor(x - mu, y)

    # Handle degenerate case: all values in both groups identical
    if (tf == 0) {
      warning("All values in both groups are identical. ",
              "Effect size is undefined (phi = 0.5 with zero variance).")
      se_val <- NA
      ci_val <- c(0, 1)
      # p_hat is already 0.5 from rbs_calc
    } else {
      # Compute CI on cstat scale via test inversion
      ci_cstat <- score_ci_wmw(phi_hat = p_hat, tf = tf, n1 = n1, n2 = n2,
                                conf.level = conf.level, correct = correct)

      # Compute descriptive SE from V_LAPH evaluated at phi_hat
      # Use boundary-safe phi for SE computation
      p_hat_safe <- pmin(pmax(p_hat, 1e-10), 1 - 1e-10)
      se_cstat <- sqrt(v_laph(p_hat_safe, tf, n1, n2))

      # Delta method SEs for other scales (descriptive only; CIs come from
      # transforming the cstat CI bounds, not from SE +/- z)
      se_rb <- 2 * se_cstat
      se_odds <- se_cstat / (1 - p_hat_safe)^2
      se_logodds <- se_cstat / (p_hat_safe * (1 - p_hat_safe))

      # Transform CI bounds to all scales
      ci_rb <- 2 * ci_cstat - 1
      ci_odds <- ci_cstat / (1 - ci_cstat)
      ci_logodds <- log(ci_cstat / (1 - ci_cstat))

      # Extract for requested scale
      se_val <- switch(ses,
                       "rb" = se_rb,
                       "cstat" = se_cstat,
                       "odds" = se_odds,
                       "logodds" = se_logodds)

      ci_val <- switch(ses,
                       "rb" = ci_rb,
                       "cstat" = ci_cstat,
                       "odds" = ci_odds,
                       "logodds" = ci_logodds)
    }

    # Create est_results for consistency (used in hypothesis testing)
    est_results <- list(
      cstat = p_hat,
      rb = r_rbs,
      odds = alpha_hat,
      logodds = eta_hat,
      se_cstat = if(exists("se_cstat")) se_cstat else NA,
      se_rb = if(exists("se_rb")) se_rb else NA,
      se_odds = if(exists("se_odds")) se_odds else NA,
      se_logodds = if(exists("se_logodds")) se_logodds else NA,
      paired = FALSE,
      n1 = n1,
      n2 = n2,
      boundary_corrected = FALSE,
      tf = tf
    )

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

    # Create est_results for consistency
    est_results <- list(
      cstat = p_hat,
      rb = r_rbs,
      odds = alpha_hat,
      logodds = eta_hat,
      se_cstat = ci_results_fisher$se_cstat,
      se_rb = ci_results_fisher$se_rb,
      se_odds = ci_results_fisher$se_odds,
      se_logodds = ci_results_fisher$se_logodds,
      paired = paired || is.null(y),
      boundary_corrected = FALSE
    )
  }

  # Get point estimate for requested effect size
  # For Agresti method with boundary correction, use the corrected estimates
  if (se_method == "agresti" && est_results$boundary_corrected) {
    est_val <- switch(ses,
                      "rb" = est_results$rb,
                      "cstat" = est_results$cstat,
                      "odds" = est_results$odds,
                      "logodds" = est_results$logodds)
  } else {
    est_val <- switch(ses,
                      "rb" = r_rbs,
                      "cstat" = p_hat,
                      "odds" = alpha_hat,
                      "logodds" = eta_hat)
  }

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

    # Method string: "<Sample Type> <Estimate Name> <test|estimate with CI>"
    method_suffix <- if (alternative != "none") "test" else "estimate with CI"
    method_desc <- paste0(sample_type, " ", ses_name, " ", method_suffix)

    # Note: SE/CI methodology details
    if (se_method == "agresti") {
      if (exists("est_results") && isTRUE(est_results$boundary_used_score_ci)) {
        note_text <- "SE: Agresti/Lehmann placement (Haldane-corrected); CI: score-type test inversion (boundary fallback)"
      } else {
        note_text <- "SE: Agresti/Lehmann placement; CI: log-odds back-transform"
      }
      if (alternative != "none") {
        note_text <- paste0(note_text, "; hypothesis test conducted on log-odds scale")
      }
    } else if (se_method == "score") {
      note_text <- "SE: Fay-Malinovsky V_LAPH; CI: score-type test inversion"
      if (correct) note_text <- paste0(note_text, " with continuity correction")
      if (alternative != "none") {
        note_text <- paste0(note_text, "; score-type hypothesis test on cstat scale")
      }
    } else {
      note_text <- "SE: Fisher z-transform"
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
      note = note_text,
      data.name = dname
    )

    # Add hypothesis test components if requested
    if (alternative != "none" && se_method == "agresti") {
      # Conduct hypothesis tests on log-odds scale for better asymptotic properties
      # Get log-odds estimate and SE from est_results
      eta_hat <- est_results$logodds
      se_eta <- est_results$se_logodds

      if (alternative %in% c("equivalence", "minimal.effect")) {
        # Transform both bounds to log-odds
        eta_low <- to_logodds(low_bound, ses)
        eta_high <- to_logodds(high_bound, ses)

        # Message if user specified custom bounds
        if (user_specified_bounds) {
          message(
            "Hypothesis test conducted on log-odds scale. ",
            "Bounds [", round(low_bound, 4), ", ", round(high_bound, 4), "] on ", ses, " scale ",
            "correspond to [", round(eta_low, 4), ", ", round(eta_high, 4), "] on log-odds scale."
          )
        }

        # Compute z-statistics on log-odds scale
        z_low <- (eta_hat - eta_low) / se_eta
        z_high <- (eta_hat - eta_high) / se_eta

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

        # Report null values on user's requested scale
        null_val <- c(low_bound, high_bound)
        names(null_val) <- c("lower bound", "upper bound")

      } else {
        # Standard alternatives (two.sided, less, greater)
        eta_null <- to_logodds(null.value, ses)

        # Message if user specified custom null (non-default)
        if (user_specified_bounds) {
          message(
            "Hypothesis test conducted on log-odds scale. ",
            "Null value ", round(null.value, 4), " on ", ses, " scale ",
            "corresponds to ", round(eta_null, 4), " on log-odds scale."
          )
        }

        # Compute z-statistic on log-odds scale
        z_stat <- (eta_hat - eta_null) / se_eta

        p_val <- switch(alternative,
                        "two.sided" = 2 * pnorm(-abs(z_stat)),
                        "less" = pnorm(z_stat),
                        "greater" = pnorm(z_stat, lower.tail = FALSE))

        # Report null value on user's requested scale
        null_val <- null.value
        names(null_val) <- ses_name
      }

      names(z_stat) <- "z"
      rval$statistic <- z_stat
      rval$p.value <- p_val
      rval$null.value <- null_val
    }

    # For Fisher method, keep existing behavior (test on requested scale)
    if (alternative != "none" && se_method == "fisher") {
      # Compute z-statistic and p-value on the requested scale
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

    # Score method hypothesis testing
    if (alternative != "none" && se_method == "score") {
      tf <- est_results$tf

      if (alternative %in% c("equivalence", "minimal.effect")) {
        # Transform bounds to cstat scale
        phi_low <- to_cstat(low_bound, ses)
        phi_high <- to_cstat(high_bound, ses)

        # Score tests against each bound
        res_low <- score_pvalue_wmw(p_hat, phi_low, tf, n1, n2,
                                     alternative = "greater", correct = correct)
        res_high <- score_pvalue_wmw(p_hat, phi_high, tf, n1, n2,
                                      alternative = "less", correct = correct)

        if (alternative == "equivalence") {
          p_val <- max(res_low$p.value, res_high$p.value)
          z_stat <- if (abs(res_low$z.statistic) < abs(res_high$z.statistic)) {
            res_low$z.statistic
          } else {
            res_high$z.statistic
          }
        } else {
          # minimal.effect
          p_val <- min(res_low$p.value, res_high$p.value)
          z_stat <- if (abs(res_low$z.statistic) < abs(res_high$z.statistic)) {
            res_low$z.statistic
          } else {
            res_high$z.statistic
          }
        }

        null_val <- c(low_bound, high_bound)
        names(null_val) <- c("lower bound", "upper bound")

      } else {
        # Standard alternatives: two.sided, less, greater
        phi_null <- to_cstat(null.value, ses)

        res <- score_pvalue_wmw(p_hat, phi_null, tf, n1, n2,
                                 alternative = alternative, correct = correct)

        z_stat <- res$z.statistic
        p_val <- res$p.value

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
