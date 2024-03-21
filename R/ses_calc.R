#' @title SES Calculation
#' @description
#' `r lifecycle::badge('stable')`
#'
#' Standardized effect size (SES), these are the effect sizes not considered SMDs.
#' @inheritParams wilcox_TOST
#' @inheritParams boot_t_TOST
#' @param mu  number indicating the value around which (a-)symmetry (for
#'   one-sample or paired samples) or shift (for independent samples) is to be
#'   estimated. See [stats::wilcox.test].
#' @details For details on the calculations in this function see `vignette("robustTOST")`.
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

