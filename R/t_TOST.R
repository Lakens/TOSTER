#' @title Two One-Sided T-tests (TOST) for Equivalence Testing
#' @description
#' `r lifecycle::badge('stable')`
#'
#' Performs equivalence testing using the Two One-Sided Tests (TOST) procedure with t-tests.
#' This function supports one-sample, two-sample (independent), and paired t-tests, providing
#' a comprehensive framework for testing equivalence or minimal effects hypotheses.
#'
#' The TOST procedure is designed for situations where you want to demonstrate that an effect
#' falls within specified bounds (equivalence testing) or exceeds specified bounds (minimal effects testing).
#'
#' @section Purpose:
#' Use this function when:
#' * You want to show that two groups are practically equivalent
#' * You need to demonstrate that an effect is at least as large as a meaningful threshold
#' * You want to test if an observed effect is too small to be of practical importance
#'
#' @param x a (non-empty) numeric vector of data values.
#' @param y an optional (non-empty) numeric vector of data values.
#' @param formula a formula of the form lhs ~ rhs where lhs is a numeric variable giving the data values and rhs either 1 for a one-sample test or a factor with two levels giving the corresponding groups. For paired tests, use the default method with x and y vectors instead of the formula method.
#' @param data an optional matrix or data frame (or similar: see model.frame) containing the variables in the formula formula. By default the variables are taken from environment(formula).
#' @param paired a logical indicating whether you want a paired t-test. Cannot be used with the formula method; use x and y vectors instead for paired tests.
#' @param var.equal a logical variable indicating whether to treat the two variances as being equal. If TRUE then the pooled variance is used to estimate the variance otherwise the Welch (or Satterthwaite) approximation to the degrees of freedom is used.
#' @param eqb Equivalence bound. Can provide 1 value (symmetric bound, negative value is taken as the lower bound) or 2 specific values that represent the upper and lower equivalence bounds.
#' @param low_eqbound lower equivalence bounds (deprecated, use `eqb` instead).
#' @param high_eqbound upper equivalence bounds (deprecated, use `eqb` instead).
#' @param hypothesis 'EQU' for equivalence (default), or 'MET' for minimal effects test.
#' @param eqbound_type Type of equivalence bound. Can be 'SMD' for standardized mean difference (i.e., Cohen's d) or 'raw' for the mean difference. Default is 'raw'. Raw is strongly recommended as SMD bounds will produce biased results.
#' @param alpha alpha level (default = 0.05)
#' @param bias_correction Apply Hedges' correction for bias (default is TRUE).
#' @param rm_correction Repeated measures correction to make standardized mean difference Cohen's d(rm). This only applies to repeated/paired samples. Default is FALSE.
#' @param mu a number indicating the true value of the mean for the two-tailed test (or difference in means if you are performing a two sample test).
#' @param glass An option to calculate Glass's delta as an alternative to Cohen's d type SMD. Default is NULL to not calculate Glass's delta, 'glass1' will use the first group's SD as the denominator whereas 'glass2' will use the 2nd group's SD.
#' @param smd_ci Method for calculating SMD confidence intervals. Methods include 'goulet', 'noncentral t' (nct), 'central t' (t), and 'normal method' (z).
#' @param subset an optional vector specifying a subset of observations to be used.
#' @param na.action a function which indicates what should happen when the data contain NAs. Defaults to getOption("na.action").
#' @param ... further arguments to be passed to or from methods.
#'
#' @details
#' For details on the calculations in this function see vignette("IntroTOSTt") & vignette("SMD_calcs").
#'
#' For two-sample tests, the test is of \eqn{\bar{x} - \bar{y}} (mean of x minus mean of y).
#' For paired samples, the test is of the difference scores (z),
#' wherein \eqn{z = x - y}, and the test is of \eqn{\bar{z}} (mean of the difference scores).
#' For one-sample tests, the test is of \eqn{\bar{x}} (mean of x).
#'
#' The output combines three statistical tests:
#' 1. A traditional two-tailed t-test (null hypothesis: difference = `mu`)
#' 2. Lower bound test (one-tailed t-test against the lower equivalence bound)
#' 3. Upper bound test (one-tailed t-test against the upper equivalence bound)
#'
#' For equivalence testing (`hypothesis = "EQU"`):
#' * **Significant TOST**: Both one-sided tests are significant (p < alpha), indicating the effect is significantly within the equivalence bounds
#'
#' For minimal effects testing (`hypothesis = "MET"`):
#' * **Significant TOST**: At least one one-sided test is significant (p < alpha), indicating the effect is significantly outside at least one of the bounds
#'
#' Notes:
#' * For equivalence testing, the equivalence bounds represent the smallest effect sizes considered meaningful.
#' * When using `eqbound_type = "SMD"`, be aware that this can produce biased results and raw bounds are generally recommended.
#' * The function provides standardized effect sizes (Cohen's d or Hedges' g) along with their confidence intervals.
#' * For paired/repeated measures designs, setting `rm_correction = TRUE` adjusts the standardized effect size calculation to account for the correlation between measures.
#'
#' @return An S3 object of class `"TOSTt"` is returned containing the following slots:
#'
#' - **TOST**: A table of class `"data.frame"` containing two-tailed t-test and both one-tailed results.
#' - **eqb**: A table of class `"data.frame"` containing equivalence bound settings.
#' - **effsize**: Table of class `"data.frame"` containing effect size estimates.
#' - **hypothesis**: String stating the hypothesis being tested.
#' - **smd**: List containing the results of the standardized mean difference calculations (e.g., Cohen's d).
#'   * Items include: d (estimate), dlow (lower CI bound), dhigh (upper CI bound), d_df (degrees of freedom for SMD), d_sigma (SE), d_lambda (non-centrality), J (bias correction), smd_label (type of SMD), d_denom (denominator calculation).
#' - **alpha**: Alpha level set for the analysis.
#' - **method**: Type of t-test.
#' - **decision**: List included text regarding the decisions for statistical inference.
#'
#' @examples
#' # Example 1: Basic Two-Sample Test
#' data(mtcars)
#' # Testing if the difference in mpg between automatic and manual
#' # transmission cars falls within ±3 mpg
#' result <- t_TOST(mpg ~ am, data = mtcars, eqb = 3)
#'
#' # Example 2: Paired Sample Test with Specific Bounds
#' data(sleep)
#' result <- t_TOST(x = sleep$extra[sleep$group == 1],
#'                 y = sleep$extra[sleep$group == 2],
#'                 paired = TRUE,
#'                 eqb = c(-0.5, 2))  # Asymmetric bounds
#'
#' # Example 3: One Sample Equivalence Test
#' result <- t_TOST(x = rnorm(30, mean = 0.1, sd = 1),
#'                 eqb = 1)
#'
#' # Example 4: Minimal Effects Test
#' result <- t_TOST(mpg ~ am,
#'                 data = mtcars,
#'                 eqb = 1.5,
#'                 hypothesis = "MET")
#'
#' @family TOST
#' @name t_TOST
#' @export t_TOST
#'
#t_TOST <- setClass("t_TOST")
t_TOST <- function(x, ...,
                   hypothesis = "EQU",
                   paired = FALSE,
                   var.equal = FALSE,
                   eqb,
                   low_eqbound,
                   high_eqbound,
                   eqbound_type = "raw",
                   alpha = 0.05,
                   bias_correction = TRUE,
                   rm_correction = FALSE,
                   glass = NULL,
                   smd_ci = c("nct", "goulet", "t", "z")){
  UseMethod("t_TOST")
}

#' @rdname t_TOST
#' @importFrom stats sd cor na.omit setNames t.test terms nlm optim optimize
#' @method t_TOST default
#' @export

# @method t_TOST default
t_TOST.default = function(x,
                          y = NULL,
                          hypothesis = c("EQU","MET"),
                          paired = FALSE,
                          var.equal = FALSE,
                          eqb,
                          low_eqbound,
                          high_eqbound,
                          eqbound_type = c("raw","SMD"),
                          alpha = 0.05,
                          mu = 0,
                          bias_correction = TRUE,
                          rm_correction = FALSE,
                          glass = NULL,
                          smd_ci = c("nct", "goulet", "t", "z"),
                          ...) {
  hypothesis = match.arg(hypothesis)
  eqbound_type = match.arg(eqbound_type)
  digits = ifelse(alpha < 0.01, 3, 2)
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

  if(hypothesis == "EQU"){
    alt_low = "greater"
    alt_high = "less"
    test_hypothesis = "Hypothesis Tested: Equivalence"

  } else if(hypothesis == "MET"){
    alt_low = "less"
    alt_high = "greater"
    test_hypothesis = "Hypothesis Tested: Minimal Effect"

  } else{
    stop("hypothesis must be set to EQU or MET")
  }

  if(eqbound_type != "raw" && eqbound_type != "SMD"){
    stop("eqbound_type must be set to raw or SMD")
  }

  if(eqbound_type == "SMD"){
    message("Warning: setting bound type to SMD produces biased results!")
  }

  if(missing(eqb) && (missing(low_eqbound) ||
     missing(high_eqbound))){
    stop("Equivalence bounds missing and must be enterered")
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


  tresult = t.test(x = x,
                   y = y,
                   paired = paired,
                   var.equal = var.equal,
                   mu = mu,
                   conf.level = 1 - alpha*2,
                   alternative = "two.sided")

  if(paired == TRUE && !missing(y)){
    i1 <- x
    i2 <- y
    data <- data.frame(i1 = i1, i2 = i2)
    data <- na.omit(data)
    colnames(data) = c("i1", "i2")
    data2 =  data
    data2$diff = data2$i2 - data2$i1 - mu

    n <- nrow(data)
    i1 <- data$i1
    i2 <- data$i2
    m1 <- mean(i1)
    m2 <- mean(i2)
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
      alpha = alpha,
      smd_ci = smd_ci
    )

  } else if(!missing(y)){

    x1 = na.omit(x)
    y1 = na.omit(y)

    n1 = length(x1)
    n2 = length(y1)

    m1 = mean(x1)
    m2 = mean(y1)-mu

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
      alpha = alpha,
      denom = denom,
      smd_ci = smd_ci
    )

  } else {

    x1 = na.omit(x)
    n1 = length(x1)
    m1 = mean(x1)-mu
    sd1 = sd(x1)

    cohen_res = d_est_one(
      n = n1,
      mu = m1,
      sd = sd1,
      type = smd_type,
      testValue = 0,
      alpha = alpha,
      smd_ci = smd_ci
    )

  }

  if(!missing(eqb)){
    if(!is.numeric(eqb) || length(eqb) > 2){
      stop(
        "eqb must be a numeric of a length of 1 or 2"
      )
    }
    if(length(eqb) == 1){
      high_eqbound = abs(eqb)
      low_eqbound = -1*abs(eqb)
    } else {
      high_eqbound = max(eqb)
      low_eqbound = min(eqb)
    }


  }

  interval_no_zero = test_interval_no_zero(c(low_eqbound, high_eqbound))

  if(interval_no_zero){
    message("Equivalence interval does not include zero.")
  }

  if (eqbound_type == 'SMD') {
    low_eqbound_d <- low_eqbound
    high_eqbound_d <- high_eqbound
    low_eqbound  <- low_eqbound * cohen_res$d_denom
    high_eqbound <- high_eqbound * cohen_res$d_denom
  } else {
    low_eqbound_d <- low_eqbound / cohen_res$d_denom
    high_eqbound_d <- high_eqbound / cohen_res$d_denom
  }

  if(hypothesis == "EQU"){
    null_hyp = paste0(round(low_eqbound,2),
                      " >= (Mean1 - Mean2) or (Mean1 - Mean2) >= ",
                      round(high_eqbound,2))
    alt_hyp = paste0(round(low_eqbound,2),
                     " < (Mean1 - Mean2) < ",
                     round(high_eqbound,2))
  } else if(hypothesis == "MET"){
    null_hyp = paste0(round(low_eqbound,2),
                      " <= (Mean1 - Mean2)  <= ",
                      round(high_eqbound,2))
    alt_hyp = paste0(round(low_eqbound,2),
                     " > (Mean1 - Mean2) or (Mean1 - Mean2)  > ",
                     round(high_eqbound,2))
  }

  low_ttest <- t.test(
    y = y,
    x = x,
    paired = paired,
    var.equal = var.equal,
    alternative = alt_low,
    mu = low_eqbound,
    conf.level = 1-alpha*2
  )

  high_ttest <- t.test(
    y = y,
    x = x,
    paired = paired,
    var.equal = var.equal,
    alternative = alt_high,
    mu = high_eqbound,
    conf.level = 1-alpha*2
  )

  if(hypothesis == "EQU"){
    pTOST = max(low_ttest$p.value,
                high_ttest$p.value) # get highest p value for TOST result
    tTOST = ifelse(abs(low_ttest$statistic) < abs(high_ttest$statistic),
                   low_ttest$statistic,
                   high_ttest$statistic) #Get lowest t-value for summary TOST result
  } else {
    pTOST = min(low_ttest$p.value,
                high_ttest$p.value) # get highest p value for TOST result
    tTOST = ifelse(abs(low_ttest$statistic) > abs(high_ttest$statistic),
                   low_ttest$statistic,
                   high_ttest$statistic) #Get lowest t-value for summary TOST result




  }

  TOST = data.frame(
    t = c(tresult$statistic,
          low_ttest$statistic,
          high_ttest$statistic),
    SE = c(tresult$stderr,
           low_ttest$stderr,
           high_ttest$stderr),
    df = c(tresult$parameter,
           low_ttest$parameter,
           high_ttest$parameter),
    p.value = c(tresult$p.value,
                low_ttest$p.value,
                high_ttest$p.value),
    row.names = c("t-test","TOST Lower","TOST Upper")
  )

  eqb = data.frame(
    type = c("Raw",cohen_res$smd_label),
    low_eq = c(low_eqbound,low_eqbound_d),
    high_eq = c(high_eqbound,high_eqbound_d)
  )

  effsize = data.frame(
    estimate = c(tresult$statistic * tresult$stderr,
                 cohen_res$d),
    SE = c(tresult$stderr,cohen_res$d_sigma),
    lower.ci = c(tresult$conf.int[1], cohen_res$dlow),
    upper.ci = c(tresult$conf.int[2], cohen_res$dhigh),
    conf.level = c((1-alpha*2),(1-alpha*2)),
    row.names = c("Raw",cohen_res$smd_label)
  )
  TOSToutcome<-ifelse(pTOST<alpha,"significant","non-significant")
  testoutcome<-ifelse(tresult$p.value<alpha,"significant","non-significant")

  # Change text based on two tailed t test if mu is not zero
  if(mu == 0){
    mu_text = "zero"
  } else {
    mu_text = mu
  }

  if(hypothesis == "EQU"){
    #format(low_eqbound, digits = 3, nsmall = 3, scientific = FALSE)
    TOST_restext = paste0("The equivalence test was ",
                          TOSToutcome,", t(",rounder_stat(tresult$parameter,
                                                   digits=digits),") = ",
                          rounder_stat(tTOST, digits = digits), ", ",
                          printable_pval(pTOST, digits = digits),sep="")
  } else {
    TOST_restext = paste0("The minimal effect test was ",
                          TOSToutcome,", t(",rounder_stat(tresult$parameter,
                                                          digits=digits),") = ",
                          rounder_stat(tTOST, digits = digits), ", ",
                          printable_pval(pTOST, digits = digits),sep="")
  }

  ttest_restext = paste0("The null hypothesis test was ",
                         testoutcome,", t(",rounder_stat(tresult$parameter, digits=2),") = ",
                         rounder_stat(tresult$statistic, digits = 3), ", ",
                         printable_pval(tresult$p.value, digits = digits),sep="")
  combined_outcome = tost_decision(
    hypothesis = hypothesis,
    alpha = alpha,
    pvalue = tresult$p.value,
    pTOST = pTOST,
    mu_text = mu_text
  )


  decision = list(
    TOST = TOST_restext,
    ttest = ttest_restext,
    combined = combined_outcome
  )

  #message(cat("Based on the equivalence test and the null-hypothesis test combined, we can conclude that the observed effect is ",combined_outcome,".",sep=""))

  if(hypothesis == "MET"){
     if(!interval_no_zero){
    if(pTOST <= tresult$p.value){
      message("MET test may have higher error rates than a nil two-tailed test. Consider wider equivalence bounds.")

      }
    }
  }

  rval = list(
    TOST = TOST,
    eqb = eqb,
    alpha = alpha,
    method = tresult$method,
    hypothesis = test_hypothesis,
    effsize = effsize,
    smd = cohen_res,
    decision = decision,
    data.name = dname,
    call = match.call()
  )

  class(rval) = "TOSTt"

  return(rval)

}

#' @rdname t_TOST
#' @method t_TOST formula
#' @export

t_TOST.formula = function(formula,
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
  y <- do.call("t_TOST", c(DATA, list(...)))
  y$data.name <- DNAME
  y

}

