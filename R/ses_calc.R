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
#'   - You want to complement results from Wilcoxon-Mann-Whitney type test
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
#'
#' The function supports three study designs:
#'   - One-sample design: Compares a single sample to a specified value
#'   - Two-sample independent design: Compares two independent groups
#'   - Paired samples design: Compares paired observations
#'
#'
#' For detailed information on calculation methods, see `vignette("robustTOST")`.
#'
#' @return A data frame containing the following information:
#'   - estimate: The effect size estimate
#'   - lower.ci: Lower bound of the confidence interval
#'   - upper.ci: Upper bound of the confidence interval
#'   - conf.level: Confidence level (1-alpha)
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
#' @family effect sizes
#' @name ses_calc
#' @export ses_calc


#ses_calc <- setClass("ses_calc")
ses_calc <- function(x, ...,
                     paired = FALSE,
                     ses = "rb",
                     alpha = 0.05){
  UseMethod("ses_calc")
}

#' @rdname ses_calc
#' @importFrom stats sd cor na.omit setNames wilcox.test terms
#' @method ses_calc default
#' @export

# @method ses_calc default
ses_calc.default = function(x,
                          y = NULL,
                          paired = FALSE,
                          ses = c("rb","odds","logodds","cstat"),
                          alpha = 0.05,
                          mu = 0,
                          ...) {

  ses = match.arg(ses)
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


  # temporary until other effect size calculations available.

  smd_type = 'd'
  denom = "z"

  if(!is.numeric(alpha) || alpha <=0 || alpha >=1){
    stop("The alpha must be a numeric value between 0 and 1")
  }

  rbs_val = np_ses(
    x = x,
    y = y,
    paired = paired,
    mu = mu,
    conf.level = 1 - alpha,
    ses = ses
  )

  ses_name = switch(ses,
                    "rb" = "Rank-Biserial Correlation",
                    "odds" = "WMW Odds",
                    "logodds" = "WMW Log-Odds",
                    "cstat" = "Concordance")

  effsize = data.frame(
    estimate = c(rbs_val$est),
    lower.ci = c(rbs_val$conf.int[1]),
    upper.ci = c(rbs_val$conf.int[2]),
    conf.level = c((1-alpha)),
    row.names = c(ses_name)
  )

  return(effsize)

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
  
  # Check for paired argument in ... and reject it
  dots <- list(...)
  if("paired" %in% names(dots)){
    if(isTRUE(dots$paired)){
      stop("cannot use 'paired' in formula method")
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
  #y$data.name <- DNAME
  y

}

