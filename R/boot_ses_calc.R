#' @title Bootstrapped Standardized Effect Size (SES) Calculation
#' @description
#' `r lifecycle::badge('maturing')`
#'
#' Calculates non-SMD standardized effect sizes with bootstrap confidence intervals.
#' This function provides more robust confidence intervals for rank-based and
#' probability-based effect size measures through resampling methods.
#'
#' @section Purpose:
#' Use this function when:
#'   - You need more robust confidence intervals for non-parametric effect sizes
#'   - You prefer resampling-based confidence intervals over asymptotic approximations
#'   - You need to quantify uncertainty in rank-based effect sizes more accurately
#'
#' @inheritParams wilcox_TOST
#' @inheritParams boot_t_TOST
#' @inheritParams boot_smd_calc
#' @param ses a character string specifying the effect size measure to calculate:
#'     - "rb": rank-biserial correlation (default)
#'     - "odds": Wilcoxon-Mann-Whitney odds
#'     - "logodds": Wilcoxon-Mann-Whitney log-odds
#'     - "cstat": concordance statistic (C-statistic/AUC)
#' @param ... further arguments to be passed to or from methods.
#'
#' @details
#' This function calculates bootstrapped confidence intervals for rank-based and probability-based
#' effect size measures. It is an extension of the `ses_calc()` function that uses resampling
#' to provide more robust confidence intervals, especially for small sample sizes.
#'
#' The function implements the following bootstrap approach:
#'   - Calculate the raw effect size using the original data
#'   - Create R bootstrap samples by resampling with replacement from the original data
#'   - Calculate the effect size for each bootstrap sample
#'   - Apply the Fisher z-transformation to stabilize variance for rank-biserial correlation values
#'   - Calculate confidence intervals using the specified method
#'   - Back-transform the confidence intervals to the original scale
#'   - Convert to the requested effect size measure (if not rank-biserial)
#'
#'
#' Three bootstrap confidence interval methods are available:
#'   - **Basic bootstrap ("basic")**: Uses the empirical distribution of bootstrap estimates
#'   - **Studentized bootstrap ("stud")**: Accounts for the variability in standard error estimates
#'   - **Percentile bootstrap ("perc")**: Uses percentiles of the bootstrap distribution directly
#'
#' The function supports three study designs:
#'   - One-sample design: Compares a single sample to a specified value
#'   - Two-sample independent design: Compares two independent groups
#'   - Paired samples design: Compares paired observations
#'
#' Note that extreme values (perfect separation between groups) can produce infinite values during
#' the bootstrapping process.
#' This happens often if the sample size is very small.
#' The function will issue a warning if this occurs, as it may affect
#' the accuracy of the confidence intervals.
#' Additionally, this affects the ability to calculate bias and SE estimates from the bootstrap samples.
#' If the number of infinite values is small (less than 10% of the bootstrap samples) then the infinite values
#' are replaced with the nearest next value (only for the SE and bias estimates, not confidence intervals).
#'
#' For detailed information on calculation methods, see `vignette("robustTOST")`.
#'
#' @return A data frame containing the following information:
#'   - estimate: The effect size estimate calculated from the original data
#'   - bias: Estimated bias (difference between original estimate and median of bootstrap estimates)
#'   - SE: Standard error estimated from the bootstrap distribution
#'   - lower.ci: Lower bound of the bootstrap confidence interval
#'   - upper.ci: Upper bound of the bootstrap confidence interval
#'   - conf.level: Confidence level (1-alpha)
#'   - boot_ci: The bootstrap confidence interval method used
#'
#' @examples
#' # Example 1: Independent groups comparison with basic bootstrap CI
#' set.seed(123)
#' group1 <- c(1.2, 2.3, 3.1, 4.6, 5.2, 6.7)
#' group2 <- c(3.5, 4.8, 5.6, 6.9, 7.2, 8.5)
#'
#' # Use fewer bootstrap replicates for a quick example
#' result <- boot_ses_calc(x = group1, y = group2,
#'                         ses = "rb",
#'                         boot_ci = "basic",
#'                         R = 99)
#'
#' # Example 2: Using formula notation to calculate concordance statistic
#' data(mtcars)
#' result <- boot_ses_calc(formula = mpg ~ am,
#'                         data = mtcars,
#'                         ses = "cstat",
#'                         boot_ci = "perc",
#'                         R = 99)
#'
#' # Example 3: Paired samples with studentized bootstrap CI
#' data(sleep)
#' with(sleep, boot_ses_calc(x = extra[group == 1],
#'                           y = extra[group == 2],
#'                           paired = TRUE,
#'                           ses = "rb",
#'                           boot_ci = "stud",
#'                           R = 99))
#'
#' # Example 4: Comparing different bootstrap CI methods
#' \dontrun{
#' # Basic bootstrap
#' basic_ci <- boot_ses_calc(x = group1, y = group2, boot_ci = "basic")
#'
#' # Percentile bootstrap
#' perc_ci <- boot_ses_calc(x = group1, y = group2, boot_ci = "perc")
#'
#' # Studentized bootstrap
#' stud_ci <- boot_ses_calc(x = group1, y = group2, boot_ci = "stud")
#'
#' # Compare the results
#' rbind(basic_ci, perc_ci, stud_ci)
#' }
#'
#' @family effect sizes
#' @name boot_ses_calc
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
  if(any(is.infinite(boots))){
    message("Bootstrapped results contain extreme results (i.e., no overlap), caution advised interpreting confidence intervals.")
  }

  # get CI on Fisher Z
  zci = switch(boot_ci,
               "stud" = stud(boots_est = boots, boots_se = boots_se,
                             se0=raw_SE, t0 = atanh(raw_ses$estimate[1L]),
                             alpha), # extreme problems with extreme
               "perc" = perc(boots, alpha),
               "basic" = basic(boots, t0 = atanh(raw_ses$estimate), alpha))
  # transform back to rcb
  rci = tanh(zci)

  rboots = tanh(boots)

  boots2 = switch(ses,
                  "rb" = rboots,
                  "cstat" = rb_to_cstat(rboots),
                  "odds" = rb_to_odds(rboots),
                  "logodds" = log(rb_to_odds(rboots)))

  # adjust for infinite observations
  if(any(is.infinite(boots2))){
    sum_inf = sum(is.infinite(boots2))
    # need to look into this more
    # might need more stringent cutoff
    if(sum_inf/R > .1){
      message("More than 10% of bootstrap estimates contain infinite values, bias and SE calculations will not be provided. Seek other robust bootstrap methods.")
    } else{
      message(paste0("A total of ",sum_inf, " bootstrapped estimates of ", R, " samples are infinite values. Bias and SE estimates are affected; proceed with caution."))
      upper_inf = max(boots2[is.finite(boots2) ])
      lower_inf = min(boots2[is.finite(boots2) ])
      boots2[boots2 == Inf ] <- upper_inf
      boots2[boots2 == -Inf] <- lower_inf
    }

  }
  # transform tp desired estimate
  ci = switch(ses,
              "rb" = rci,
              "cstat" = rb_to_cstat(rci),
              "odds" = rb_to_odds(rci),
              "logodds" = log(rb_to_odds(rci)))

  est2 = switch(ses,
                "rb" = raw_ses$estimate,
                "cstat" = rb_to_cstat(raw_ses$estimate),
                "odds" = rb_to_odds(raw_ses$estimate),
                "logodds" = log(rb_to_odds(raw_ses$estimate)))

  effsize = data.frame(
    estimate = est2,
    bias = est2 - median(boots2),
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
  y <- do.call("boot_ses_calc", c(DATA, list(...)))
  #y$data.name <- DNAME
  y

}
