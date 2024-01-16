#' @title Bootstrap SES Calculation
#' @description
#' `r lifecycle::badge('maturing')`
#'
#' Standardized effect size (SES), these are the effect sizes not considered SMDs.
#' @inheritParams wilcox_TOST
#' @inheritParams boot_t_TOST
#' @details For details on the calculations in this function see `vignette("robustTOST")`.
#' @return A data frame containing the standardized effect size.
#' @examples
#' \dontrun{
#' boot_ses_calc(formula = extra ~ group, data = sleep, paired = TRUE, ses = "r")
#' }

# Bootstrap -------

#' @name boot_ses_calc
#' @family effect sizes
#' @export boot_ses_calc

#ses_calc <- setClass("ses_calc")
boot_ses_calc <- function(x, ...,
                          paired = FALSE,
                          ses = "rb",
                          alpha = 0.05,
                          boot_ci = c("basic","stud","perc"),
                          R = 1999){
  UseMethod("boot_ses_calc")
}


#' @rdname boot_ses_calc
#' @method boot_ses_calc default
#' @export
# @method ses_calc default
boot_ses_calc.default = function(x,
                                 y = NULL,
                                 paired = FALSE,
                                 ses = c("rb","odds","logodds","cstat"),
                                 alpha = 0.05,
                                 boot_ci = c("basic","stud", "perc"),
                                 R = 1999,
                                 ...) {
  boot_ci = match.arg(boot_ci)
  ses = match.arg(ses)

  # paired ----
  if(paired == TRUE && !is.null(y)){
    i1 <- x
    i2 <- y
    data <- data.frame(x = i1, y = i2)
    data <- na.omit(data)
    nd = nrow(data)
    maxw <- (nd^2 + nd) / 2
    raw_SE = sqrt((2 * nd^3 + 3 * nd^2 + nd) / 6) / maxw
    raw_ses = ses_calc(x = data$x,
                       y = data$y,
                       paired = paired,
                       ses = "rb",
                       alpha = alpha)

    boots = c()
    boots_se = c()
    for(i in 1:R){
      sampler = sample(1:nrow(data), replace = TRUE)
      res_boot = ses_calc(x = data$x[sampler],
                          y = data$y[sampler],
                          paired = paired,
                          ses = "rb",
                          alpha = alpha)
      boots = c(boots, atanh(res_boot$estimate))
      nd <- nrow(data)
      maxw <- (nd^2 + nd) / 2
      rfSE <- sqrt((2 * nd^3 + 3 * nd^2 + nd) / 6) / maxw
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
    raw_SE = sqrt((n1 + n2 + 1) / (3 * n1 * n2))

    raw_ses = ses_calc(x = i1,
                       y = i2,
                       paired = paired,
                       ses = "rb",
                       alpha = alpha)

    boots = c()
    boots_se = c()
    for(i in 1:R){
      sampler = sample(1:nrow(data), replace = TRUE)
      boot_dat = data[sampler,]
      x_boot = subset(boot_dat,
                      group == "x")
      y_boot = subset(boot_dat,
                      group == "y")
      res_boot = ses_calc(x = x_boot$values,
                          y = y_boot$values,
                          paired = paired,
                          ses = "rb",
                          alpha = alpha)
      boots = c(boots, atanh(res_boot$estimate))
      n1 = nrow(x_boot)
      n2 = nrow(y_boot)
      rfSE <- sqrt((n1 + n2 + 1) / (3 * n1 * n2))
      boots_se = c(boots_se, rfSE)

    }

  } else {
    # one-sample -----
    x1 = na.omit(x)
    n1 = nd =  length(x1)
    maxw <- (nd^2 + nd) / 2
    raw_SE = sqrt((2 * nd^3 + 3 * nd^2 + nd) / 6) / maxw
    raw_ses = ses_calc(x = x1,
                       paired = paired,
                       ses = "rb",
                       alpha = alpha)

    boots = c()
    boots_se = c()
    for(i in 1:R){
      sampler = sample(1:length(x1), replace = TRUE)
      x_boot = x1[sampler]

      res_boot = ses_calc(x = x_boot,
                          paired = paired,
                          ses = "rb",
                          alpha = alpha)
      boots = c(boots, atanh(res_boot$estimate))
      nd <- length(x_boot)
      maxw <- (nd^2 + nd) / 2
      rfSE <- sqrt((2 * nd^3 + 3 * nd^2 + nd) / 6) / maxw
      boots_se = c(boots_se, rfSE)
    }

  }

  zci = switch(boot_ci,
               "stud" = stud(boots_est = boots, boots_se = boots_se,
                             se0=raw_SE, t0 = atanh(raw_ses$estimate[1L]),
                             alpha),
               "perc" = perc(boots, alpha),
               "basic" = basic(boots, t0 = raw_ses$estimate, alpha))
  rci = tanh(zci)

  rboots = tanh(boots)

  boots2 = switch(ses,
                  "rb" = rboots,
                  "cstat" = rb_to_cstat(rboots),
                  "odds" = rb_to_odds(rboots),
                  "logodds" = log(rb_to_odds(rboots)))

  ci = switch(ses,
              "rb" = rci,
              "cstat" = rb_to_cstat(rci),
              "odds" = rb_to_odds(rci),
              "logodds" = log(rb_to_odds(rci)))

  effsize = data.frame(
    estimate = raw_ses$estimate,
    bias = raw_ses$estimate - median(boots2),
    SE = sd(boots2),
    lower.ci = ci[1],
    upper.ci = ci[2],
    conf.level = c((1-alpha)),
    boot_ci = boot_ci,
    row.names = c(raw_ses$ses_label)
  )


  return(effsize)

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
