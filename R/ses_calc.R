#' @title SES Calculation
#' @description
#' `r lifecycle::badge('stable')`
#'
#' Standardized effect size (SES), these are the effect sizes not considered SMDs.
#' @inheritParams wilcox_TOST
#' @inheritParams boot_t_TOST
#' @param boot_ci type of bootstrap confidence interval. Options include studentized (stud), empirical/basic (basic) and percentile (perc) confidence intervals.
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

# Bootstrap -------

#' @rdname ses_calc

#ses_calc <- setClass("ses_calc")
boot_ses_calc <- function(x, ...,
                          paired = FALSE,
                          ses = "rb",
                          alpha = 0.05,
                          boot_ci = c("basic","perc"),
                          R = 1999){
  UseMethod("boot_ses_calc")
}



# @method ses_calc default
boot_ses_calc.default = function(x,
                                 y = NULL,
                                 paired = FALSE,
                                 ses = c("rb","odds","logodds","cstat"),
                                 alpha = 0.05,
                                 boot_ci = c("basic", "perc"),
                                 R = 1999,
                                 ...) {
  boot_ci = match.arg(boot_ci)
  ses = match.arg(ses)
  if(paired == TRUE && !missing(y)){
    i1 <- x
    i2 <- y
    data <- data.frame(x = i1, y = i2)
    data <- na.omit(data)
    raw_ses = ses_calc(x = data$x,
                       y = data$y,
                       paired = paired,
                       ses = ses,
                       alpha = alpha)

    boots = c()

    for(i in 1:R){
      sampler = sample(1:nrow(data), replace = TRUE)
      res_boot = ses_calc(x = data$x[sampler],
                          y = data$y[sampler],
                          paired = paired,
                          ses = ses,
                          alpha = alpha)
      boots = c(boots, res_boot$estimate)
    }


  } else if(!missing(y)){

    i1 <- na.omit(x)
    i2 <- na.omit(y)
    data <- data.frame(values = c(i1,i2),
                       group = c(rep("x",length(i1)),
                                 rep(i2)))

    raw_ses = ses_calc(x = i1,
                       y = i2,
                       paired = paired,
                       ses = ses,
                       alpha = alpha)

    boots = c()

    for(i in 1:R){
      sampler = sample(1:nrow(data), replace = TRUE)
      boot_dat = data[sampler,]
      x_boot = subset(boot_dat,
                      group == "x")
      y_boot = subset(boot_dat,
                      group == "y")
      res_boot = ses_calc(x = x_boot,
                          y = y_boot,
                          paired = paired,
                          ses = ses,
                          alpha = alpha)
      boots = c(boots, res_boot$estimate)
    }

  } else {

    x1 = na.omit(x)
    n1 = length(x1)
    raw_ses = ses_calc(x = x1,
                       paired = paired,
                       ses = ses,
                       alpha = alpha)

    boots = c()

    for(i in 1:R){
      sampler = sample(1:nrow(x1), replace = TRUE)
      x_boot = x1[sampler]

      res_boot = ses_calc(x = x_boot,
                          paired = paired,
                          ses = ses,
                          alpha = alpha)
      boots = c(boots, res_boot$estimate)
    }

  }

  ci = switch(boot_ci,
              "perc" = perc(boots, alpha),
              "basic" = basic(boots, t0 = raw_ses$estimate, alpha))

  effsize = data.frame(
    estimate = raw_ses$estimate,
    bias = raw_ses$estimate - median(boots),
    SE = sd(boots),
    lower.ci = ci[1],
    upper.ci = ci[2],
    conf.level = c((1-alpha)),
    boot_ci = boot_ci,
    row.names = c(raw_ses$ses_label)
  )


  return(effsize)

}

#' @rdname ses_calc
#' @method ses_calc formula
#' @export

boot_ses_calc.formula = function(formula,
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
  y <- do.call("boot_ses_calc", c(DATA, list(...)))
  #y$data.name <- DNAME
  y

}
