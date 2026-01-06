#' @title Bootstrapped Standardized Mean Difference (SMD) Calculation
#' @description
#' `r lifecycle::badge('maturing')`
#'
#' Calculates standardized mean differences (SMDs) with bootstrap confidence intervals.
#' This function provides more robust confidence intervals for Cohen's d, Hedges' g,
#' and other SMD measures through resampling methods.
#'
#' @section Purpose:
#' Use this function when:
#'
#'   * You need more robust confidence intervals for standardized mean differences
#'   * You want to account for non-normality or heterogeneity in your effect size estimates
#'   * Sample sizes are small or standard error approximations may be unreliable
#'   * You prefer resampling-based confidence intervals over parametric approximations
#'   * You need to quantify uncertainty in SMD estimates more accurately
#'
#' @inheritParams boot_t_TOST
#' @param mu null value to adjust the calculation. If non-zero, the function calculates x-y-mu (default = 0).
#' @param ... further arguments to be passed to or from methods.
#'
#' @details
#' This function calculates bootstrapped confidence intervals for standardized mean differences.
#' It is an extension of the `smd_calc()` function that uses resampling to provide more robust
#' confidence intervals, especially for small sample sizes or when data violate assumptions
#' of parametric methods.
#'
#' The function implements the following bootstrap approach:
#'   * Calculate the raw SMD and its standard error using the original data
#'   * Create R bootstrap samples by resampling with replacement from the original data
#'   * Calculate the SMD and its standard error for each bootstrap sample
#'   * Calculate confidence intervals using the specified method
#'
#' Three bootstrap confidence interval methods are available:
#'   - **Studentized bootstrap ("stud")**: Accounts for the variability in standard error estimates. Usually provides the most accurate coverage probability and is set as the default.
#'   - **Basic bootstrap ("basic")**: Uses the empirical distribution of bootstrap estimates. Simple approach that works well for symmetric distributions.
#'   - **Percentile bootstrap ("perc")**: Uses percentiles of the bootstrap distribution directly. More robust to skewness in the bootstrap distribution.
#'
#' The function supports various SMD variants:
#'   * Classic standardized mean difference (bias_correction = FALSE)
#'   * Bias-corrected version (bias_correction = TRUE)
#'   * Glass's delta: Uses only one group's standard deviation as the denominator (glass = "glass1" or "glass2")
#'   * Repeated measures d: Accounts for correlation in paired designs (rm_correction = TRUE)
#'
#' The function supports three study designs:
#'   * One-sample design: Standardizes the difference between the sample mean and zero (or other specified value)
#'   * Two-sample independent design: Standardizes the difference between two group means
#'   * Paired samples design: Standardizes the mean difference between paired observations
#'
#' For detailed information on calculation methods, see `vignette("SMD_calcs")`.
#'
#' @return A data frame containing the following information:
#'   * estimate: The SMD calculated from the original data
#'   * bias: Estimated bias (difference between original estimate and median of bootstrap estimates)
#'   * SE: Standard error estimated from the bootstrap distribution
#'   * lower.ci: Lower bound of the bootstrap confidence interval
#'   * upper.ci: Upper bound of the bootstrap confidence interval
#'   * conf.level: Confidence level (1-alpha)
#'   * boot_ci: The bootstrap confidence interval method used
#'
#' @examples
#' # Example 1: Independent groups comparison with studentized bootstrap CI
#' set.seed(123)
#' group1 <- rnorm(30, mean = 100, sd = 15)
#' group2 <- rnorm(30, mean = 110, sd = 18)
#'
#' # Use fewer bootstrap replicates for a quick example
#' result <- boot_smd_calc(x = group1, y = group2,
#'                       boot_ci = "stud",
#'                       R = 999)
#'
#' # Example 2: Using formula notation with basic bootstrap and Hedges' g
#' df <- data.frame(
#'   value = c(group1, group2),
#'   group = factor(rep(c("A", "B"), each = 30))
#' )
#' result <- boot_smd_calc(formula = value ~ group,
#'                       data = df,
#'                       boot_ci = "basic",
#'                       bias_correction = TRUE,
#'                       R = 999)
#'
#' # Example 3: Paired samples with percentile bootstrap
#' set.seed(456)
#' before <- rnorm(30)
#' after <- rnorm(30)
#' result <- boot_smd_calc(x = before,
#'                       y = after,
#'                       paired = TRUE,
#'                       boot_ci = "perc",
#'                       R = 999)
#'
#' # Example 4: Glass's delta with homogeneous variances
#' set.seed(456)
#' control <- rnorm(25, mean = 50, sd = 10)
#' treatment <- rnorm(25, mean = 60, sd = 10)
#' result <- boot_smd_calc(x = control,
#'                       y = treatment,
#'                       glass = "glass1",
#'                       boot_ci = "stud",
#'                       R = 999)
#'
#' @family effect sizes
#' @name boot_smd_calc
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
      sampler = sample(1:length(x1), replace = TRUE)
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


  return(effsize = effsize)

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
  y <- do.call("boot_smd_calc", c(DATA, list(...)))
  #y$data.name <- DNAME
  y

}
