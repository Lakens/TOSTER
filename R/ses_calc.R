#' @title SES Calculation
#' @description Standardized effect size (SES), these are the effect sizes not considered SMDs.
#' @param x a (non-empty) numeric vector of data values.
#' @param y an optional (non-empty) numeric vector of data values.
#' @param formula a formula of the form lhs ~ rhs where lhs is a numeric variable giving the data values and rhs either 1 for a one-sample or paired test or a factor with two levels giving the corresponding groups. If lhs is of class "Pair" and rhs is 1, a paired test is done.
#' @param data an optional matrix or data frame (or similar: see model.frame) containing the variables in the formula formula. By default the variables are taken from environment(formula).
#' @param paired a logical indicating whether you want to calculate a paired test.
#' @param alpha alpha level (default = 0.05)
#' @param mu  number indicating the value around which (a-)symmetry (for
#'   one-sample or paired samples) or shift (for independent samples) is to be
#'   estimated. See [stats::wilcox.test].
#' @param ses Standardized effect size. Default is "rb" for rank-biserial
#' correlation. Options also include "cstat" for concordance probability, or
#' "odds" for Wilcoxon-Mann-Whitney odds (otherwise known as Agresti's
#' generalized odds ratio).
#' @param subset an optional vector specifying a subset of observations to be used.
#' @param na.action a function which indicates what should happen when the data contain NAs. Defaults to getOption("na.action").
#' @param ...  further arguments to be passed to or from methods.
#' @details For details on the calculations in this function see vignette("robustTOST").
#' @return A data frame containing the standardized effect size.
#' @examples
#' \dontrun{
#' ses_calc(formula = extra ~ group, data = sleep, paired = TRUE, ses = "r")
#' }
#' @name ses_calc
#' @family effect sizes
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
                          ses = c("rb","odds","cstat"),
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

