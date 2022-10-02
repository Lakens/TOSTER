#' @title TOST with Wilcoxon Signed Rank test
#' @description A function for TOST using the non-parametric methods of the Wilcoxon signed rank test. This function uses the normal approximation and applies continuity correction automatically.
#' @param x a (non-empty) numeric vector of data values.
#' @param y an optional (non-empty) numeric vector of data values.
#' @param formula a formula of the form lhs ~ rhs where lhs is a numeric variable giving the data values and rhs either 1 for a one-sample or paired test or a factor with two levels giving the corresponding groups. If lhs is of class "Pair" and rhs is 1, a paired test is done.
#' @param data an optional matrix or data frame (or similar: see model.frame) containing the variables in the formula formula. By default the variables are taken from environment(formula).
#' @param paired a logical indicating whether you want to calculate a paired test.
#' @param eqb Equivalence bound. Can provide 1 value (negative value is taken as the lower bound) or 2 specific values that represent the upper and lower equivalence bounds.
#' @param low_eqbound Lower equivalence bounds. Deprecated use eqb.
#' @param high_eqbound Upper equivalence bounds. Deprecated use eqb.
#' @param hypothesis 'EQU' for equivalence (default), or 'MET' for minimal effects test, the alternative hypothesis.
#' @param ses Rank-biserial (rb), odds (odds), and concordance probablity (cstat).
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
#' @return An S3 object of class
#'   \code{"TOSTnp"} is returned containing the following slots:
#' \describe{
#'   \item{\code{"TOST"}}{A table of class \code{"data.frame"} containing two-tailed wilcoxon signed rank test and both one-tailed results.}
#'   \item{\code{"eqb"}}{A table of class \code{"data.frame"} containing equivalence bound settings.}
#'   \item{\code{"effsize"}}{ table of class \code{"data.frame"} containing effect size estimates.}
#'   \item{\code{"hypothesis"}}{String stating the hypothesis being tested.}
#'   \item{\code{"smd"}}{List containing information on standardized effect size.}
#'   \item{\code{"alpha"}}{Alpha level set for the analysis.}
#'   \item{\code{"method"}}{Type of non-parametric test.}
#'   \item{\code{"decision"}}{List included text regarding the decisions for statistical inference.}
#' }
#' @section References:
#' David F. Bauer (1972). Constructing confidence sets using rank statistics. Journal of the American Statistical Association 67, 687–690. doi: 10.1080/01621459.1972.10481279.
#'
#' Myles Hollander and Douglas A. Wolfe (1973). Nonparametric Statistical Methods. New York: John Wiley & Sons. Pages 27–33 (one-sample), 68–75 (two-sample). Or second edition (1999).
#' @importFrom stats wilcox.test
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

  rbs_val = np_ses(
    x = x,
    y = y,
    paired = paired,
    mu = mu,
    conf.level = 1 - alpha * 2,
    ses = ses
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
                 rbs_val$est),
    lower.ci = c(tresult$conf.int[1], rbs_val$conf.int[1]),
    upper.ci = c(tresult$conf.int[2], rbs_val$conf.int[2]),
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
    TOST_restext = paste0("The equivalence test was ",TOSToutcome," ",names(tresult$statistic), " = ", format(tTOST, digits = 3, nsmall = 3, scientific = FALSE),", p = ",format(pTOST, digits = 3, nsmall = 3, scientific = TRUE),sep="")
  } else {
    TOST_restext = paste0("The minimal effect test was ",TOSToutcome," ",names(tresult$statistic), " = ", format(tTOST, digits = 3, nsmall = 3, scientific = FALSE),", p = ",format(pTOST, digits = 3, nsmall = 3, scientific = TRUE),sep="")
  }

  ttest_restext = paste0("The null hypothesis test was ",testoutcome," ",names(tresult$statistic), " = ", format(tresult$statistic, digits = 3, nsmall = 3, scientific = FALSE),", p = ",format(tresult$p.value, digits = 3, nsmall = 3, scientific = TRUE),sep="")
  if (hypothesis == "EQU"){
    if(tresult$p.value <= alpha && pTOST <= alpha){
      combined_outcome <- paste0("NHST: reject null significance hypothesis that the effect is equal to ",mu_text," \n",
                                 "TOST: reject null equivalence hypothesis")
    }
    if(tresult$p.value < alpha && pTOST > alpha){
      combined_outcome <- paste0("NHST: reject null significance hypothesis that the effect is equal to ",mu_text," \n",
                                 "TOST: don't reject null equivalence hypothesis")
      # paste0("statistically different from ",mu_text," and not statistically equivalent")
    }
    if(tresult$p.value > alpha && pTOST <= alpha){
      combined_outcome <- paste0("NHST: don't reject null significance hypothesis that the effect is equal to ",mu_text," \n",
                                 "TOST: reject null equivalence hypothesis")
      #paste0("statistically not different from ",mu_text," and statistically equivalent")
    }
    if(tresult$p.value > alpha && pTOST > alpha){
      combined_outcome <- paste0("NHST: don't reject null significance hypothesis that the effect is equal to ",mu_text," \n",
                                 "TOST: don't reject null equivalence hypothesis")
      #paste0("statistically not different from ",mu_text," and not statistically equivalent")
    }
  } else {
    if(tresult$p.value <= alpha && pTOST <= alpha){
      combined_outcome <- paste0("NHST: reject null significance hypothesis that the effect is equal to ",mu_text," \n",
                                 "TOST: reject null MET hypothesis")
      #paste0("statistically different from ",mu_text," and statistically greater than the minimal effect threshold")
    }
    if(tresult$p.value < alpha && pTOST > alpha){
      combined_outcome <- paste0("NHST: reject null significance hypothesis that the effect is equal to ",mu_text," \n",
                                 "TOST: don't reject null MET hypothesis")
      #paste0("statistically different from ",mu_text," but not statistically greater than the minimal effect threshold")
    }
    if(tresult$p.value > alpha && pTOST <= alpha){
      combined_outcome <- paste0("NHST: don't reject null significance hypothesis that the effect is equal to ",mu_text," \n",
                                 "TOST: reject null MET hypothesis")
      #paste0("statistically not different from ",mu_text," and statistically greater than the minimal effect threshold")
    }
    if(tresult$p.value > alpha && pTOST > alpha){
      combined_outcome <- paste0("NHST: don't reject null significance hypothesis that the effect is equal to ",mu_text," \n",
                                 "TOST: don't reject null MET hypothesis")
      #paste0("statistically not different from ",mu_text," and not statistically greater than the minimal effect threshold")
    }
  }


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

