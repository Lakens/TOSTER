#' @title Standardized Mean Difference (SMD) Calculation
#' @description
#' `r lifecycle::badge('stable')`
#'
#' Calculates standardized mean difference (SMD) effect sizes and their confidence intervals
#' from raw data. This function focuses solely on effect size estimation without performing
#' hypothesis tests.
#'
#' @section Purpose:
#' Use this function when:
#'   * You need to calculate standardized effect sizes (Cohen's d, Hedges' g, Glass's delta)
#'   * You want confidence intervals for your effect size estimates
#'   * You need effect sizes for meta-analysis or reporting
#'   * You want to compare effect sizes across different studies or measures
#'   * You don't need the hypothesis testing components of the TOST functions
#'
#' @inheritParams t_TOST
#' @inheritParams boot_t_TOST
#' @param ... further arguments to be passed to or from methods.
#'
#' @details
#' This function calculates standardized mean differences (SMD) for various study designs:
#'
#'   * One-sample design: Standardizes the difference between the sample mean and zero (or other specified value)
#'   * Two-sample independent design: Standardizes the difference between two group means
#'   * Paired samples design: Standardizes the mean difference between paired observations
#'
#' The function supports multiple SMD variants:
#'   * Cohen's d: Classic standardized mean difference (bias_correction = FALSE)
#'   * Hedges' g: Bias-corrected version of Cohen's d (bias_correction = TRUE)
#'   * Glass's delta: Uses only one group's standard deviation as the denominator (glass = "glass1" or "glass2")
#'   * Repeated measures d: Accounts for correlation in paired designs (rm_correction = TRUE)
#'
#' Different confidence interval calculation methods are available:
#'   * "nct": Uses the noncentral t-distribution (most accurate in most cases)
#'   * "goulet": Uses the Goulet-Pelletier method
#'   * "t": Uses the central t-distribution
#'   * "z": Uses the normal distribution
#'
#' Note that unlike the t_TOST and related functions, smd_calc only calculates effect sizes and
#' their confidence intervals without performing hypothesis tests.
#'
#' For detailed information on calculation methods, see `vignette("SMD_calcs")`.
#'
#' @return A data frame containing the following information:
#'   * estimate: The standardized mean difference estimate (Cohen's d, Hedges' g, or Glass's delta)
#'   * SE: Standard error of the estimate
#'   * lower.ci: Lower bound of the confidence interval
#'   * upper.ci: Upper bound of the confidence interval
#'   * conf.level: Confidence level (1-alpha)
#'
#' @examples
#' # Example 1: Independent groups comparison (Cohen's d)
#' set.seed(123)
#' group1 <- rnorm(30, mean = 100, sd = 15)
#' group2 <- rnorm(30, mean = 110, sd = 18)
#' smd_calc(x = group1, y = group2, bias_correction = FALSE)
#'
#' # Example 2: Independent groups with formula notation (Hedges' g)
#' df <- data.frame(
#'   value = c(group1, group2),
#'   group = factor(rep(c("A", "B"), each = 30))
#' )
#' smd_calc(formula = value ~ group, data = df)
#'
#' # Example 3: Paired samples with repeated measures correction
#' before <- c(5.1, 4.8, 6.2, 5.7, 6.0, 5.5, 4.9, 5.8)
#' after <- c(5.6, 5.2, 6.7, 6.1, 6.5, 5.8, 5.3, 6.2)
#' smd_calc(x = before, y = after, paired = TRUE, rm_correction = TRUE)
#'
#' # Example 4: Glass's delta (using only first group's SD)
#' smd_calc(x = group1, y = group2, glass = "glass1")
#'
#' @family effect sizes
#' @name smd_calc
#' @export smd_calc

#smd_calc <- setClass("smd_calc")
smd_calc <- function(x, ...,
                     paired = FALSE,
                     var.equal = FALSE,
                     alpha = 0.05,
                     bias_correction = TRUE,
                     rm_correction = FALSE,
                     glass = NULL,
                     smd_ci = c("nct", "goulet", "t", "z")){
  UseMethod("smd_calc")
}

#' @rdname smd_calc
#' @importFrom stats sd cor na.omit setNames t.test terms nlm optim optimize
#' @method smd_calc default
#' @export

# @method smd_calc default
smd_calc.default = function(x,
                            y = NULL,
                            paired = FALSE,
                            var.equal = FALSE,
                            alpha = 0.05,
                            mu = 0,
                            bias_correction = TRUE,
                            rm_correction = FALSE,
                            glass = NULL,
                            smd_ci = c("nct", "goulet", "t", "z"),
                            ...) {

  if(is.null(glass)){
    glass = "no"
  }
  smd_ci = match.arg(smd_ci)

  if(bias_correction){
    smd_type = 'g'
  } else {
    smd_type = 'd'
  }

  if(is.null(y)){
    sample_type = "One Sample"
  } else if(paired == TRUE) {
    sample_type = "Paired Sample"
  } else {
    sample_type = "Two Sample"
  }

  if(glass == "glass1" || glass == "glass2"){
    if(glass == "glass1"){
      denom = "glass1"
    }

    if(glass == "glass2"){
      denom = "glass2"
    }
  } else{
    if(sample_type != "Two Sample" ){
      if(rm_correction){
        denom = "rm"
      } else {
        denom = "z"
      }
    } else{
      denom = "d"
    }
  }

  if(!is.numeric(alpha) || alpha <=0 || alpha >=1){
    stop("The alpha must be a numeric value between 0 and 1")
  }

  if (!is.null(y)) {
    dname <- paste(deparse(substitute(x)), "and",
                   deparse(substitute(y)))
  }
  else {
    dname <- deparse(substitute(x))
  }

  if(paired == TRUE && !missing(y)){
    i1 <- x
    i2 <- y
    data <- data.frame(i1 = i1, i2 = i2)
    data <- na.omit(data)
    colnames(data) = c("i1", "i2")
    data2 =  data
    data2$diff = data2$i2 - data2$i1

    n <- nrow(data)
    i1 <- data$i1
    i2 <- data$i2
    m1 <- mean(i1)
    m2 <- mean(i2) + mu
    sd1  <- sd(i1)
    sd2  <- sd(i2)
    r12 <- cor(i1, i2)

    # Calculate Cohens d
    cohen_res = d_est_pair(
      n = n,
      m1 = m1,
      m2 = m2,
      sd1 = sd1,
      sd2 = sd2,
      r12 = r12,
      type = smd_type,
      denom = denom,
      alpha = alpha/2,
      smd_ci = smd_ci
    )

  } else if(!missing(y)){

    x1 = na.omit(x)
    y1 = na.omit(y)
    n1 = length(x1)
    n2 = length(y1)
    m1 = mean(x1)
    m2 = mean(y1) + mu
    sd1 = sd(x1)
    sd2 = sd(y1)

    cohen_res = d_est_ind(
      n1 = n1,
      n2 = n2,
      m1 = m1,
      m2 = m2,
      sd1 = sd1,
      sd2 = sd2,
      type = smd_type,
      var.equal = var.equal,
      alpha = alpha/2,
      denom = denom,
      smd_ci = smd_ci
    )

  } else {

    x1 = na.omit(x)
    n1 = length(x1)
    m1 = mean(x1) + mu
    sd1 = sd(x1)

    cohen_res = d_est_one(
      n = n1,
      mu = m1,
      sd = sd1,
      type = smd_type,
      testValue = 0,
      alpha = alpha/2,
      smd_ci = smd_ci
    )

  }

  effsize = data.frame(
    estimate = c(cohen_res$d),
    SE = c(cohen_res$d_sigma),
    lower.ci = c(cohen_res$dlow),
    upper.ci = c( cohen_res$dhigh),
    conf.level = c((1-alpha)),
    row.names = c(cohen_res$smd_label)
  )


  return(effsize)

}

#' @rdname smd_calc
#' @method smd_calc formula
#' @export

smd_calc.formula = function(formula,
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
  y <- do.call("smd_calc", c(DATA, list(...)))
  #y$data.name <- DNAME
  y

}

