#' @title TOST with Wilcoxon-Mann-Whitney tests
#' @description
#' `r lifecycle::badge('stable')`
#'
#' A function for TOST using the non-parametric methods of the Wilcoxon-Mann-Whitney family of tests.
#' This function uses the normal approximation and applies a continuity correction automatically.
#'
#' @inheritParams t_TOST
#' @param mu  number indicating the value around which (a-)symmetry (for
#'   one-sample or paired samples) or shift (for independent samples) is to be
#'   estimated. See [stats::wilcox.test].
#' @param eqb Equivalence bound. Can provide one value (symmetric bound) or two specific values that represent the lower and upper equivalence bounds.
#'   Like the `mu` argument, `eqb` is specified on the **raw scale** of the original data (e.g., the scale of the location shift).
#'   This parameter is independent of the `ses` argument, which only affects the type of standardized effect size that is reported.
#' @param ses Standardized effect size. Default is "rb" for rank-biserial
#' correlation. Options also include "cstat" for concordance probability, or
#' "odds" for Wilcoxon-Mann-Whitney odds (otherwise known as Agresti's
#' generalized odds ratio). Note that `ses` only determines which effect size is calculated and does not affect the equivalence bounds (`eqb`).
#' @details
#' For details on the calculations in this function see `vignette("robustTOST")`. For details on the Wilcoxon-Mann-Whitney tests see [stats::wilcox.test].
#'
#' If only x is given, or if both x and y are given and paired is TRUE,
#' a Wilcoxon signed rank test of the null that the distribution of x (in the one sample case)
#' or of x - y (in the paired two sample case) is symmetric about mu is performed.
#'
#' Otherwise, if both x and y are given and paired is FALSE,
#' a Wilcoxon rank sum test (equivalent to the Mann-Whitney test: see the Note) is carried out.
#' In this case, the null hypothesis is that the distributions of x and y differ by a location shift.
#'
#' @return An S3 object of class
#'   `"TOSTnp"` is returned containing the following slots:
#'
#'   - "TOST": A table of class `"data.frame"` containing two-tailed wilcoxon signed rank test and both one-tailed results.
#'   - "eqb": A table of class `"data.frame"` containing equivalence bound settings.
#'   - "effsize":  table of class `"data.frame"` containing effect size estimates.
#'   - "hypothesis": String stating the hypothesis being tested.
#'   - "smd": List containing information on standardized effect size.
#'   - "alpha": Alpha level set for the analysis.
#'   - "method": Type of non-parametric test.
#'   - "decision": List included text regarding the decisions for statistical inference.
#'
#' @examples
#' data(mtcars)
#' wilcox_TOST(mpg ~ am,
#' data = mtcars,
#' eqb = 3)
#' @section References:
#' David F. Bauer (1972). Constructing confidence sets using rank statistics. Journal of the American Statistical Association 67, 687–690. doi: 10.1080/01621459.1972.10481279.
#'
#' Myles Hollander and Douglas A. Wolfe (1973). Nonparametric Statistical Methods. New York: John Wiley & Sons. Pages 27–33 (one-sample), 68–75 (two-sample). Or second edition (1999).
#' @importFrom stats wilcox.test
#' @family Robust tests
#' @family TOST
#' @name wilcox_TOST
#' @export wilcox_TOST

#wilcox_TOST <- setClass("wilcox_TOST")
wilcox_TOST <- function(x, ...,
                        hypothesis = "EQU",
                        paired = FALSE,
                        eqb,
                        low_eqbound,
                        high_eqbound,
                        ses = "rb",
                        alpha = 0.05){
  UseMethod("wilcox_TOST")
}

#' @rdname wilcox_TOST
#' @importFrom stats sd cor na.omit setNames wilcox.test terms
#' @method wilcox_TOST default
#' @export

# @method wilcox_TOST default
wilcox_TOST.default = function(x,
                          y = NULL,
                          hypothesis = "EQU",
                          paired = FALSE,
                          eqb,
                          low_eqbound,
                          high_eqbound,
                          ses = c("rb","odds", "logodds", "cstat"),
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

  if(missing(eqb) && (missing(low_eqbound) ||
                      missing(high_eqbound))){
    stop("Equivalence bounds missing and must be enterered")
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

  if(!is.numeric(alpha) || alpha <=0 || alpha >=1){
    stop("The alpha must be a numeric value between 0 and 1")
  }


  tresult = wilcox.test(x = x,
                   y = y,
                   paired = paired,
                   exact = FALSE,
                   conf.int = TRUE,
                   mu = mu,
                   conf.level = 1 - alpha*2,
                   alternative = "two.sided")

  rbs_val = ses_calc(
    x = x,
    y = y,
    paired = paired,
    mu = mu,
    alpha = alpha * 2,
    ses = ses,
    se_method = "fisher",
    output = "data.frame"
  )

  if(hypothesis == "EQU"){
    null_hyp = paste0(round(low_eqbound,2),
                      " >= (Median1 - Median2) or (Median1 - Median2) >= ",
                      round(high_eqbound,2))
    alt_hyp = paste0(round(low_eqbound,2),
                     " < (Median1 - Median2) < ",
                     round(high_eqbound,2))
  } else if(hypothesis == "MET"){
    null_hyp = paste0(round(low_eqbound,2),
                      " <= (Median1 - Median2)  <= ",
                      round(high_eqbound,2))
    alt_hyp = paste0(round(low_eqbound,2),
                     " > (Median1 - Median2) or (Median1 - Median2)  > ",
                     round(high_eqbound,2))
  }

  low_ttest <- wilcox.test(
    y = y,
    x = x,
    paired = paired,
    exact = FALSE,
    alternative = alt_low,
    mu = low_eqbound,
    conf.level = 1-alpha*2
  )

  high_ttest <- wilcox.test(
    y = y,
    x = x,
    paired = paired,
    exact = FALSE,
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

    if(!interval_no_zero){
      if(pTOST <= tresult$p.value){
        message("MET test may have higher error rates than a nil two-tailed test. Consider wider equivalence bounds.")
      }
    }
  }

  TOST = data.frame(
    statistic = c(tresult$statistic,
          low_ttest$statistic,
          high_ttest$statistic),
    p.value = c(tresult$p.value,
                low_ttest$p.value,
                high_ttest$p.value),
    row.names = c("NHST","TOST Lower","TOST Upper")
  )

  eqb = c(low_eqbound, high_eqbound)

  if(!paired && is.null(y)){
    raw_name = "pseudomedian"
  } else{
    raw_name = "Median of Differences"
  }
  ses_name = switch(ses,
                    "rb" = "Rank-Biserial Correlation",
                    "odds" = "WMW Odds",
                    "cstat" = "Concordance")

  effsize = data.frame(
    estimate = c(tresult$estimate,
                 rbs_val$estimate),
    lower.ci = c(tresult$conf.int[1], rbs_val$lower.ci),
    upper.ci = c(tresult$conf.int[2], rbs_val$upper.ci),
    conf.level = c((1-alpha*2),(1-alpha*2)),
    row.names = c(raw_name,ses_name)
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
    TOST_restext = paste0("The equivalence test was ",TOSToutcome," ",
                          names(tresult$statistic), " = ",
                          format(tTOST, digits = 3,
                                 nsmall = 3, scientific = FALSE),", p = ",
                          format(pTOST, digits = 3,
                                 nsmall = 3, scientific = TRUE),sep="")
  } else {
    TOST_restext = paste0("The minimal effect test was ",TOSToutcome," ",
                          names(tresult$statistic), " = ",
                          format(tTOST, digits = 3,
                                 nsmall = 3, scientific = FALSE),", p = ",
                          format(pTOST, digits = 3,
                                 nsmall = 3, scientific = TRUE),sep="")
  }

  ttest_restext = paste0("The null hypothesis test was ",testoutcome," ",
                         names(tresult$statistic), " = ",
                         format(tresult$statistic,
                                digits = 3,
                                nsmall = 3, scientific = FALSE),", p = ",
                         format(tresult$p.value, digits = 3,
                                nsmall = 3, scientific = TRUE),sep="")
  combined_outcome = tost_decision(hypothesis = hypothesis,
                                    alpha = alpha,
                                    pvalue = tresult$p.value,
                                    pTOST = pTOST,
                                    mu_text = mu_text)


  decision = list(
    TOST = TOST_restext,
    test = ttest_restext,
    combined = combined_outcome
  )


  rval = list(
    TOST = TOST,
    eqb = eqb,
    alpha = alpha,
    method = tresult$method,
    hypothesis = test_hypothesis,
    effsize = effsize,
    seff = rbs_val,
    decision = decision,
    data.name = dname,
    call = match.call()
  )

  class(rval) = "TOSTnp"

  return(rval)

}

#' @rdname wilcox_TOST
#' @method wilcox_TOST formula
#' @export

wilcox_TOST.formula = function(formula,
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
  y <- do.call("wilcox_TOST", c(DATA, list(...)))
  y$data.name <- DNAME
  y

}

