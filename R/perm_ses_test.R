#' @title Permutation Test for Standardized Effect Sizes
#' @description
#' `r lifecycle::badge('maturing')`
#'
#' Performs permutation-based hypothesis tests for non-parametric effect sizes including
#' concordance probability (c-statistic), rank-biserial correlation, WMW odds, and WMW log-odds.
#' This function provides exact or randomization-based inference without relying on asymptotic
#' approximations.
#'
#' @section Purpose:
#' Use this function when:
#'   - You want a distribution-free hypothesis test for rank-based effect sizes
#'   - Sample sizes are small where asymptotic approximations (from [ses_calc()]) may be unreliable
#'   - You want exact p-values without relying on normal approximations
#'   - You need equivalence or minimal effect testing with permutation methods
#'
#' @inheritParams t_TOST
#' @param ses a character string specifying the effect size measure to test:
#'     - "rb": rank-biserial correlation (default)
#'     - "odds": Wilcoxon-Mann-Whitney odds
#'     - "logodds": Wilcoxon-Mann-Whitney log-odds
#'     - "cstat": concordance statistic (C-statistic, equivalent to the area under the ROC curve)
#' @param alpha significance level (default = 0.05).
#' @param mu a number or vector specifying the null hypothesis value(s) on the scale of the
#'   effect size specified by `ses`:
#'     - For standard alternatives: a single value (default = 0 for rb/logodds, 0.5 for cstat, 1 for odds)
#'     - For equivalence/minimal.effect: two values representing the lower and upper bounds
#' @param alternative a character string specifying the alternative hypothesis:
#'     - "two.sided" (default): Test whether effect differs from mu
#'     - "less": Test whether effect is less than mu
#'     - "greater": Test whether effect is greater than mu
#'     - "equivalence": Test whether effect is between specified bounds (TOST)
#'     - "minimal.effect": Test whether effect is outside specified bounds
#' @param R the number of permutations. Default is NULL, which computes all exact permutations
#'     if feasible. If R is specified and is less than the maximum number of possible permutations,
#'     randomization sampling is used instead.
#' @param symmetric a logical variable indicating whether to assume symmetry in the two-sided test.
#'     If `TRUE` (default) then the symmetric permutation p-value is computed, otherwise the
#'     equal-tail permutation p-value is computed. Only relevant for `alternative = "two.sided"`.
#' @param p_method the method for computing permutation p-values:
#'     - `NULL` (default): Automatically selects "exact" for exact permutation tests and
#'       "plusone" for randomization tests.
#'     - `"exact"`: Uses b/R where b is the count of permutation statistics at least as extreme
#'       as observed.
#'     - `"plusone"`: Uses (b+1)/(R+1), which guarantees p > 0 and provides exact Type I error
#'       control (Phipson & Smyth, 2010).
#' @param keep_perm logical. If `TRUE` (default), the permutation distribution of the effect size
#'     is stored in the output. Set to `FALSE` for large datasets to save memory.
#' @param ... further arguments to be passed to or from methods.
#'
#' @details
#' This function performs permutation-based hypothesis tests for rank-based and probability-based
#' effect size measures. It provides distribution-free inference that is exact (when all permutations
#' are enumerated) or approximate (when randomization sampling is used).
#'
#' The test uses an unstudentized approach where the test statistic is the effect size itself.
#' This avoids the need to compute standard errors for each permutation and sidesteps boundary
#' case issues that can arise with studentized approaches.
#'
#' The available effect size measures are:
#'
#'   - **Rank-biserial correlation ("rb")**: Ranges from -1 to 1. No effect corresponds to 0.
#'
#'   - **Concordance statistic ("cstat")**: Ranges from 0 to 1. No effect corresponds to 0.5.
#'
#'   - **WMW odds ("odds")**: Ranges from 0 to infinity. No effect corresponds to 1.
#'
#'   - **WMW log-odds ("logodds")**: Ranges from -infinity to infinity. No effect corresponds to 0.
#'
#' ## Permutation Procedure
#'
#' For two-sample (independent groups) tests, permutation is performed by randomly reassigning
#' observations to groups. For paired and one-sample tests, permutation is performed by randomly
#' flipping the sign of differences (equivalent to randomly swapping x and y within pairs).
#'
#' If the number of possible permutations is less than or equal to R (or R is NULL), all exact
#' permutations are computed and a message is printed.
#'
#' ## P-Value Computation
#'
#' For different alternatives, the p-values are calculated as follows:
#'   - "two.sided": Either symmetric (proportion of |T*| >= |T - null|) or equal-tail
#'       (2 * min of one-sided p-values) depending on the `symmetric` argument
#'   - "less": Proportion of T* <= T
#'   - "greater": Proportion of T* >= T
#'   - "equivalence": Maximum of two one-sided p-values, following the IU-NPC approach
#'       of Arboretti et al. (2021)
#'   - "minimal.effect": Minimum of two one-sided p-values
#'
#' ## Confidence Intervals
#'
#' Confidence intervals are computed using the percentile method on the permutation distribution.
#' For bounded scales, intervals are clamped to the valid range (e.g., [-1, 1] for rb,
#' [0, 1] for cstat).
#'
#' ## Ties
#'
#' Ties are handled using the midrank method (average ranks), which is the standard approach
#' for Wilcoxon-type statistics. No special correction is applied.
#'
#' @return A list with class `"htest"` containing:
#'   - statistic: the observed effect size (test statistic)
#'   - parameter: the number of permutations used (R)
#'   - p.value: the permutation p-value
#'   - conf.int: a permutation percentile confidence interval
#'   - estimate: the effect size estimate
#'   - null.value: the specified null hypothesis value(s)
#'   - alternative: a character string describing the alternative hypothesis
#'   - method: a character string indicating the test method
#'   - data.name: a character string giving the name(s) of the data
#'   - call: the matched call
#'   - R: the requested number of permutations
#'   - R.used: the actual number of permutations used
#'   - perm.dist: (if keep_perm = TRUE) the permutation distribution of the effect size
#'
#' @examples
#' # Example 1: Two-sample test with rank-biserial correlation
#' set.seed(123)
#' x <- c(1.2, 2.3, 3.1, 4.6, 5.2, 6.7)
#' y <- c(3.5, 4.8, 5.6, 6.9, 7.2, 8.5)
#' perm_ses_test(x = x, y = y, ses = "rb")
#'
#' # Example 2: Concordance statistic with one-sided test
#' perm_ses_test(x = x, y = y, ses = "cstat", alternative = "less")
#'
#' # Example 3: Paired samples
#' data(sleep)
#' with(sleep, perm_ses_test(x = extra[group == 1],
#'                            y = extra[group == 2],
#'                            paired = TRUE,
#'                            ses = "rb",
#'                            R = 999))
#'
#' # Example 4: Equivalence testing
#' perm_ses_test(x = x, y = y, ses = "rb",
#'               alternative = "equivalence",
#'               mu = c(-0.8, 0.8), R = 999)
#'
#' # Example 5: Formula interface
#' data(mtcars)
#' perm_ses_test(formula = mpg ~ am, data = mtcars,
#'               ses = "rb", R = 999)
#'
#' @references
#' Arboretti, R., Pesarin, F., & Salmaso, L. (2021). A unified approach to permutation testing
#' for equivalence. *Statistical Methods & Applications*, 30, 1033-1052.
#'
#' Agresti, A. (1980). Generalized odds ratios for ordinal data. *Biometrics*, 36, 59-67.
#'
#' Phipson, B., & Smyth, G. K. (2010). Permutation P-values should never be zero: calculating
#' exact P-values when permutations are randomly drawn. *Statistical Applications in Genetics
#' and Molecular Biology*, 9(1), Article 39.
#'
#' @family Robust tests
#' @name perm_ses_test
#' @export perm_ses_test

perm_ses_test <- function(x, ...,
                          paired = FALSE,
                          ses = c("rb", "cstat", "odds", "logodds"),
                          alpha = 0.05,
                          mu = NULL,
                          alternative = c("two.sided", "less", "greater",
                                          "equivalence", "minimal.effect"),
                          R = NULL,
                          symmetric = TRUE,
                          p_method = NULL,
                          keep_perm = TRUE) {
  UseMethod("perm_ses_test")
}

#' @rdname perm_ses_test
#' @importFrom stats na.omit setNames quantile terms
#' @method perm_ses_test default
#' @export

perm_ses_test.default <- function(x,
                                  y = NULL,
                                  paired = FALSE,
                                  ses = c("rb", "cstat", "odds", "logodds"),
                                  alpha = 0.05,
                                  mu = NULL,
                                  alternative = c("two.sided", "less", "greater",
                                                  "equivalence", "minimal.effect"),
                                  R = NULL,
                                  symmetric = TRUE,
                                  p_method = NULL,
                                  keep_perm = TRUE,
                                  ...) {

  ses <- match.arg(ses)
  alternative <- match.arg(alternative)

  # Input validation
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a numeric value between 0 and 1")
  }

  if (!is.null(R)) {
    if (length(R) != 1 || !is.finite(R) || R < 1) {
      stop("'R' must be NULL (for exact permutation) or a positive integer")
    }
    R <- as.integer(R)
  }

  # Set default null value based on effect size type
  if (is.null(mu)) {
    mu <- switch(ses,
                 "rb" = 0,
                 "cstat" = 0.5,
                 "odds" = 1,
                 "logodds" = 0)
  }

  # Null value for the permutation distribution (exchangeability null)
  null_under_H0 <- switch(ses,
                           "rb" = 0,
                           "cstat" = 0.5,
                           "odds" = 1,
                           "logodds" = 0)

  # Handle equivalence/minimal.effect bounds
  if (alternative %in% c("equivalence", "minimal.effect")) {
    if (length(mu) == 1) {
      stop("For equivalence or minimal.effect testing, 'mu' must specify two bounds (e.g., mu = c(-0.3, 0.3))")
    }
    if (length(mu) != 2) {
      stop("For equivalence or minimal.effect testing, 'mu' must be a vector of length 2")
    }
    mu <- sort(mu)
    low_bound <- mu[1]
    high_bound <- mu[2]
    conf.level <- 1 - alpha * 2
  } else {
    if (length(mu) != 1) {
      stop("'mu' must be a single value for this alternative")
    }
    conf.level <- 1 - alpha
  }

  # Data name extraction
  if (!is.null(y)) {
    dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
  } else {
    dname <- deparse(substitute(x))
  }

  # Handle missing values
  if (is.null(y)) {
    x <- na.omit(x)
    n1 <- length(x)
  } else if (paired) {
    data <- data.frame(x = x, y = y)
    data <- na.omit(data)
    x <- data$x
    y <- data$y
    n1 <- nrow(data)
  } else {
    x <- na.omit(x)
    y <- na.omit(y)
    n1 <- length(x)
    n2 <- length(y)
  }

  # --- Helper: compute effect size for a permuted sample ---
  compute_perm_es <- function(x_perm, y_perm, paired_flag, ses_type) {
    if (paired_flag) {
      r_rbs <- rbs_calc(x = y_perm, y = x_perm, mu = 0, paired = TRUE)
    } else {
      r_rbs <- rbs_calc(x = x_perm, y = y_perm, mu = 0, paired = FALSE)
    }
    p_hat <- (r_rbs + 1) / 2
    switch(ses_type,
           "rb" = r_rbs,
           "cstat" = p_hat,
           "odds" = p_hat / (1 - p_hat),
           "logodds" = log(p_hat / (1 - p_hat)))
  }

  # --- Compute observed effect size ---
  if (is.null(y)) {
    # One-sample: compare x to 0 using paired convention with y = rep(0, n)
    r_rbs_obs <- rbs_calc(x = rep(0, n1), y = x, mu = 0, paired = TRUE)
  } else if (paired) {
    r_rbs_obs <- rbs_calc(x = y, y = x, mu = 0, paired = TRUE)
  } else {
    r_rbs_obs <- rbs_calc(x = x, y = y, mu = 0, paired = FALSE)
  }

  p_hat_obs <- (r_rbs_obs + 1) / 2
  obs_es <- switch(ses,
                   "rb" = r_rbs_obs,
                   "cstat" = p_hat_obs,
                   "odds" = p_hat_obs / (1 - p_hat_obs),
                   "logodds" = log(p_hat_obs / (1 - p_hat_obs)))

  # --- Generate permutation distribution ---
  if (is.null(y) || paired) {
    # One-sample or paired: sign-flip permutation
    if (is.null(y)) {
      # One-sample: differences are x - 0 = x
      d <- x
    } else {
      # Paired: use original x and y for sign-flip
      d <- NULL  # We'll work with x and y directly
    }

    perm_result <- perm_signs(if (is.null(y)) n1 else n1, R)
    signs_matrix <- perm_result$signs
    R_used <- perm_result$R.used
    exact_perm <- perm_result$exact

    # Set p_method default
    if (is.null(p_method)) {
      p_method <- if (exact_perm) "exact" else "plusone"
    } else {
      p_method <- match.arg(p_method, c("exact", "plusone"))
    }

    PERM_ES <- numeric(R_used)
    for (i in seq_len(R_used)) {
      flip <- signs_matrix[i, ] < 0  # TRUE where sign is -1

      if (is.null(y)) {
        # One-sample: flip sign of x values
        x_perm <- ifelse(flip, -x, x)
        y_perm <- rep(0, n1)
      } else {
        # Paired: swap x and y where flip is TRUE
        x_perm <- ifelse(flip, y, x)
        y_perm <- ifelse(flip, x, y)
      }
      PERM_ES[i] <- compute_perm_es(x_perm, y_perm, paired_flag = TRUE, ses_type = ses)
    }

  } else {
    # Two-sample: shuffle group labels
    perm_result <- perm_groups(n1, length(y), R)
    idx_x_matrix <- perm_result$idx_x
    R_used <- perm_result$R.used
    exact_perm <- perm_result$exact

    # Set p_method default
    if (is.null(p_method)) {
      p_method <- if (exact_perm) "exact" else "plusone"
    } else {
      p_method <- match.arg(p_method, c("exact", "plusone"))
    }

    all_data <- c(x, y)
    N <- n1 + length(y)

    PERM_ES <- numeric(R_used)
    for (i in seq_len(R_used)) {
      idx_x <- idx_x_matrix[i, ]
      idx_y <- setdiff(seq_len(N), idx_x)
      x_perm <- all_data[idx_x]
      y_perm <- all_data[idx_y]
      PERM_ES[i] <- compute_perm_es(x_perm, y_perm, paired_flag = FALSE, ses_type = ses)
    }
  }

  # --- Warn about Inf values on odds/logodds scale ---
  if (ses %in% c("odds", "logodds")) {
    if (!is.finite(obs_es)) {
      warning("Complete separation detected. Effect size is ",
              ifelse(obs_es > 0, "Inf", "-Inf"),
              " on the ", ses, " scale. ",
              "Consider using 'cstat' or 'rb' scale for interpretable results.")
    }
    n_inf <- sum(!is.finite(PERM_ES))
    if (n_inf > 0) {
      message(n_inf, " of ", R_used, " permutations produced infinite values on the ",
              ses, " scale.")
    }
  }

  # --- Compute p-values ---
  if (alternative == "two.sided") {
    if (symmetric) {
      b <- sum(abs(PERM_ES - null_under_H0) >= abs(obs_es - mu))
      perm.pval <- compute_perm_pval(b, R_used, p_method)
    } else {
      b_low <- sum(PERM_ES <= obs_es)
      b_high <- sum(PERM_ES >= obs_es)
      p_low <- compute_perm_pval(b_low, R_used, p_method)
      p_high <- compute_perm_pval(b_high, R_used, p_method)
      perm.pval <- min(2 * min(p_low, p_high), 1)
    }

  } else if (alternative == "less") {
    b <- sum(PERM_ES <= obs_es)
    perm.pval <- compute_perm_pval(b, R_used, p_method)

  } else if (alternative == "greater") {
    b <- sum(PERM_ES >= obs_es)
    perm.pval <- compute_perm_pval(b, R_used, p_method)

  } else if (alternative == "equivalence") {
    # IU-NPC: shift observed to each bound, compare against perm null
    obs_shifted_low <- obs_es - low_bound
    obs_shifted_high <- obs_es - high_bound
    perm_shifted <- PERM_ES - null_under_H0

    b_low <- sum(perm_shifted >= obs_shifted_low)
    b_high <- sum(perm_shifted <= obs_shifted_high)

    p_low <- compute_perm_pval(b_low, R_used, p_method)
    p_high <- compute_perm_pval(b_high, R_used, p_method)
    perm.pval <- max(p_low, p_high)

  } else if (alternative == "minimal.effect") {
    obs_shifted_low <- obs_es - low_bound
    obs_shifted_high <- obs_es - high_bound
    perm_shifted <- PERM_ES - null_under_H0

    b_low <- sum(perm_shifted <= obs_shifted_low)
    b_high <- sum(perm_shifted >= obs_shifted_high)

    p_low <- compute_perm_pval(b_low, R_used, p_method)
    p_high <- compute_perm_pval(b_high, R_used, p_method)
    perm.pval <- min(p_low, p_high)
  }

  # --- Compute confidence intervals ---
  if (alternative == "two.sided") {
    perm.cint <- quantile(PERM_ES, c(alpha / 2, 1 - alpha / 2), names = FALSE)
  } else if (alternative == "less") {
    lower_min <- switch(ses,
                        "rb" = -1, "cstat" = 0, "odds" = 0, "logodds" = -Inf)
    perm.cint <- c(lower_min, quantile(PERM_ES, conf.level, names = FALSE))
  } else if (alternative == "greater") {
    upper_max <- switch(ses,
                        "rb" = 1, "cstat" = 1, "odds" = Inf, "logodds" = Inf)
    perm.cint <- c(quantile(PERM_ES, 1 - conf.level, names = FALSE), upper_max)
  } else if (alternative %in% c("equivalence", "minimal.effect")) {
    perm.cint <- quantile(PERM_ES, c(alpha, 1 - alpha), names = FALSE)
  }

  # Clamp to valid range
  if (ses == "rb") {
    perm.cint <- pmax(pmin(perm.cint, 1), -1)
  } else if (ses == "cstat") {
    perm.cint <- pmax(pmin(perm.cint, 1), 0)
  } else if (ses == "odds") {
    perm.cint <- pmax(perm.cint, 0)
  }

  attr(perm.cint, "conf.level") <- conf.level

  # --- Build output ---
  ses_name <- switch(ses,
                     "rb" = "Rank-Biserial Correlation",
                     "cstat" = "Concordance",
                     "odds" = "WMW Odds",
                     "logodds" = "WMW Log-Odds")

  if (alternative %in% c("equivalence", "minimal.effect")) {
    null.value <- c(low_bound, high_bound)
    names(null.value) <- c("lower bound", "upper bound")
  } else {
    null.value <- mu
    names(null.value) <- ses_name
  }

  estimate <- obs_es
  names(estimate) <- ses_name

  stat <- obs_es
  names(stat) <- "observed"

  param <- R_used
  names(param) <- "R"

  if (exact_perm) {
    method_prefix <- "Exact Permutation"
  } else {
    method_prefix <- "Randomization Permutation"
  }

  if (is.null(y)) {
    sample_type <- "(One-Sample)"
  } else if (paired) {
    sample_type <- "(Paired)"
  } else {
    sample_type <- "(Two-Sample)"
  }

  method <- paste(method_prefix, ses_name, "Test", sample_type)

  rval <- list(
    statistic = stat,
    parameter = param,
    p.value = perm.pval,
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
    rval$perm.dist <- PERM_ES
  }

  class(rval) <- "htest"
  return(rval)
}

#' @rdname perm_ses_test
#' @method perm_ses_test formula
#' @export

perm_ses_test.formula <- function(formula,
                                  data,
                                  subset,
                                  na.action, ...) {

  if (missing(formula) ||
      (length(formula) != 3L) ||
      (length(attr(terms(formula[-2L]), "term.labels")) != 1L)) {
    stop("'formula' missing or incorrect")
  }

  dots <- list(...)
  if ("paired" %in% names(dots)) {
    if (isTRUE(dots$paired)) {
      message("Using 'paired = TRUE' with the formula interface is not recommended. ",
              "Please ensure your data is sorted appropriately.")
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

  DATA <- setNames(split(mf[[response]], g), c("x", "y"))
  y <- do.call("perm_ses_test", c(DATA, list(...)))
  y$data.name <- DNAME
  y
}
