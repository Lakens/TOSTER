#' @title Bootstrapped SMD Calculation
#' @description
#' `r lifecycle::badge('maturing')`
#'
#' A function to only calculate standardized mean differences with bootstrap confidence intervals.
#' @inheritParams boot_t_TOST
#' @param mu Null value. Deviating from zero will give the x-y-mu.
#' @details For details on the calculations in this function see vignette("SMD_calcs").
#' @return A data frame containing the SMD estimates.
#' @examples
#' \dontrun{
#' boot_smd_calc(formula = extra ~ group,data = sleep, paired = TRUE, smd_ci = "nct")
#' }
#' @name boot_smd_calc
#' @family effect sizes
#' @export boot_smd_calc

# Bootstrap -------

#smd_calc <- setClass("smd_calc")
boot_smd_calc <- function(x, ...,
                          paired = FALSE,
                          var.equal = FALSE,
                          alpha = 0.05,
                          bias_correction = TRUE,
                          rm_correction = FALSE,
                          glass = NULL,
                          boot_ci = c("stud","basic","perc"),
                          R = 1999){
  UseMethod("boot_smd_calc")
}



#' @rdname boot_smd_calc
#' @method boot_smd_calc default
#' @export
boot_smd_calc.default = function(x,
                                 y = NULL,
                                 paired = FALSE,
                                 var.equal = FALSE,
                                 alpha = 0.05,
                                 mu = 0,
                                 bias_correction = TRUE,
                                 rm_correction = FALSE,
                                 glass = NULL,
                                 boot_ci = c("stud","basic","perc"),
                                 R = 1999,
                                 ...) {
  boot_ci = match.arg(boot_ci)
  if(paired == TRUE && !missing(y)){
    i1 <- x
    i2 <- y
    data <- data.frame(x = i1, y = i2)
    data <- na.omit(data)
    raw_smd = smd_calc(x = data$x,
                       y = data$y,
                       paired = paired,
                       var.equal = var.equal,
                       alpha = alpha,
                       mu = mu,
                       bias_correction = bias_correction,
                       rm_correction = rm_correction,
                       glass = glass,
                       smd_ci = "z")

    boots = c()
    boots_se = c()

    for(i in 1:R){
      sampler = sample(1:nrow(data), replace = TRUE)
      res_boot = smd_calc(x = data$x[sampler],
                          y = data$y[sampler],
                          paired = paired,
                          var.equal = var.equal,
                          alpha = alpha,
                          mu = mu,
                          bias_correction = bias_correction,
                          rm_correction = rm_correction,
                          glass = glass,
                          smd_ci = "z")
      boots = c(boots, res_boot$estimate)
      boots_se = c(boots_se, res_boot$SE)
    }


  } else if(!missing(y)){

    i1 <- na.omit(x)
    i2 <- na.omit(y)
    data <- data.frame(values = c(i1,i2),
                       group = c(rep("x",length(i1)),
                                 rep("y",length(i2))))

    raw_smd = smd_calc(x = i1,
                       y = i2,
                       paired = paired,
                       var.equal = var.equal,
                       alpha = alpha,
                       mu = mu,
                       bias_correction = bias_correction,
                       rm_correction = rm_correction,
                       glass = glass,
                       smd_ci = "z")

    boots = c()
    boots_se = c()
    for(i in 1:R){
      sampler = sample(1:nrow(data), replace = TRUE)
      boot_dat = data[sampler,]
      x_boot = subset(boot_dat,
                      group == "x")
      y_boot = subset(boot_dat,
                      group == "y")
      res_boot = smd_calc(x = x_boot$values,
                          y = y_boot$values,
                          paired = paired,
                          var.equal = var.equal,
                          alpha = alpha,
                          mu = mu,
                          bias_correction = bias_correction,
                          rm_correction = rm_correction,
                          glass = glass,
                          smd_ci = "z")
      boots = c(boots, res_boot$estimate)
      boots_se = c(boots_se, res_boot$SE)
    }

  } else {

    x1 = na.omit(x)
    n1 = length(x1)
    raw_smd = smd_calc(x = x1,
                       paired = paired,
                       var.equal = var.equal,
                       alpha = alpha,
                       mu = mu,
                       bias_correction = bias_correction,
                       rm_correction = rm_correction,
                       glass = glass,
                       smd_ci = "z")

    boots = c()
    boots_se = c()
    for(i in 1:R){
      sampler = sample(1:nrow(x1), replace = TRUE)
      x_boot = x1[sampler]

      res_boot = smd_calc(x = x_boot,
                          paired = paired,
                          var.equal = var.equal,
                          alpha = alpha,
                          mu = mu,
                          bias_correction = bias_correction,
                          rm_correction = rm_correction,
                          glass = glass,
                          smd_ci = "z")
      boots = c(boots, res_boot$estimate)
      boots_se = c(boots_se, res_boot$SE)
    }

  }

  ci = switch(boot_ci,
              "stud" = stud(boots_est = boots, boots_se = boots_se,
                            se0=raw_smd$SE[1L], t0 = raw_smd$estimate[1L],
                            alpha),
              "perc" = perc(boots, alpha),
              "basic" = basic(boots, t0 = raw_smd$estimate[1L], alpha=alpha))

  effsize = data.frame(
    estimate = raw_smd$estimate,
    bias = raw_smd$estimate - median(boots, na.rm=TRUE),
    SE = sd(boots),
    lower.ci = ci[1],
    upper.ci = ci[2],
    conf.level = c((1-alpha)),
    boot_ci = boot_ci,
    row.names = row.names(raw_smd)
  )


  return(list(effsize = effsize,
              boots = list(est = boots,
                           se = boots_se)))

}

#' @rdname boot_smd_calc
#' @method boot_smd_calc formula
#' @export

boot_smd_calc.formula = function(formula,
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
  y <- do.call("boot_smd_calc", c(DATA, list(...)))
  #y$data.name <- DNAME
  y

}
