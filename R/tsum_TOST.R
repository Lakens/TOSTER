#' @title TOSTt with Summary Statistics
#' @description A function for TOST with all types of t-tests from summary statistics.
#' @param m1 mean of group 1
#' @param m2 mean of group 2
#' @param sd1 standard deviation of group 1
#' @param sd2 standard deviation of group 2
#' @param n1 sample size in group 1
#' @param n2 sample size in group 2
#' @param r12 correlation of dependent variable between group 1 and group 2
#' @param paired a logical indicating whether you want a paired t-test.
#' @param var.equal  a logical variable indicating whether to treat the two variances as being equal. If TRUE then the pooled variance is used to estimate the variance otherwise the Welch (or Satterthwaite) approximation to the degrees of freedom is used.
#' @param low_eqbound lower equivalence bounds
#' @param high_eqbound upper equivalence bounds
#' @param hypothesis 'EQU' for equivalence (default), or 'MET' for minimal effects test, the alternative hypothesis.
#' @param eqbound_type Type of equivalence bound. Can be set to "SMD" for standardized mean difference (i.e., Cohen's d) or  "raw" for the mean difference. Default is "raw".
#' @param alpha alpha level (default = 0.05)
#' @param bias_correction Apply Hedges' correction for bias (default is TRUE).
#' @param rm_correction Repeated measures correction to make standardized mean difference Cohen's d(rm). This only applies to repeated/paired samples. Default is FALSE.
#' @param mu a number indicating the true value of the mean for the two tailed test (or difference in means if you are performing a two sample test).
#' @return An S3 object of class
#'   \code{"TOSTt"} is returned containing the following slots:
#' \describe{
#'   \item{\code{"TOST"}}{A table of class \code{"data.frame"} containing two-tailed t-test and both one-tailed results.}
#'   \item{\code{"eqb"}}{A table of class \code{"data.frame"} containing equivalence bound settings.}
#'   \item{\code{"effsize"}}{ table of class \code{"data.frame"} containing effect size estimates}
#'   \item{\code{"hypothesis"}}{String stating the hypothesis being tested}
#'   \item{\code{"smd"}}{List containing the results of the standardized mean difference calculations (e.g., Cohen's d).Items include: d (estimate), dlow (lower CI bound), dhigh (upper CI bound), d_df (degrees of freedom for SMD), d_sigma (SE), d_lambda (non-centrality), J (bias correction), smd_label (type of SMD), d_denom (denominator calculation)}
#'   \item{\code{"alpha"}}{Alpha level set for the analysis.}
#'   \item{\code{"method"}}{Type of t-test.}
#'   \item{\code{"decision"}}{List included text regarding the decisions for statistical inference.}
#' }
#' @importFrom stats na.omit setNames terms
#' @name tsum_TOST
#' @export tsum_TOST

#t_TOST <- setClass("t_TOST")
tsum_TOST <- function(m1,
                      sd1,
                      n1,
                      m2 = NULL,
                      sd2 = NULL,
                      n2 = NULL,
                      r12 = NULL,
                      hypothesis = "EQU",
                      paired = FALSE,
                      var.equal = FALSE,
                      low_eqbound,
                      high_eqbound,
                      mu = 0,
                      eqbound_type = "raw",
                      alpha = 0.05,
                      bias_correction = TRUE,
                      rm_correction = FALSE){
  if(bias_correction){
    smd_type = 'g'
  } else {
    smd_type = 'd'
  }

  if(rm_correction){
    denom = "rm"
  } else {
    denom = "z"
  }

  if(is.null(n2) || is.null(m2) || is.null()){
    sample_type = "One Sample"
  } else if(paired == TRUE && !is.null(r12)) {
    sample_type = "Paired Sample"
  } else if (paired == TRUE && is.null(r12)){
    stop("paired == TRUE but r12 not provided. Must provide correlation.")
  } else {
    sample_type = "Two Sample"
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

  if(missing(low_eqbound) ||
     missing(high_eqbound)){
    stop("Equivalence bounds missing and must be enterered")
  }

  if(!is.numeric(alpha) || alpha <=0 || alpha >=1){
    stop("The alpha must be a numeric value between 0 and 1")
  }

  tresult = tsum_test(m1 = m1, sd1 = sd1, n1 = n1,
                      m2 = m2, sd2 = sd2, n2 = n2,
                      paired = paired,
                      var.equal = var.equal,
                      mu = mu,
                      conf.level = 1 - alpha * 2,
                      alternative = "two.sided")

  if(paired == TRUE && !missing(r12)){

    cohen_res = d_est_pair(
      n = n,
      m1 = m1,
      m2 = m2,
      sd1 = sd1,
      sd2 = sd2,
      r12 = r12,
      type = smd_type,
      denom = denom,
      alpha = alpha
    )

  } else if(sample_type == "Two Sample"){

    cohen_res = d_est_ind(
      n1 = n1,
      n2 = n2,
      m1 = m1,
      m2 = m2,
      sd1 = sd1,
      sd2 = sd2,
      type = smd_type,
      var.equal = var.equal,
      alpha = alpha
    )

  } else {
    cohen_res = d_est_one(
      n = n1,
      mu = m1,
      sd = sd1,
      type = smd_type,
      testValue = mu,
      alpha = alpha
    )
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

  low_ttest <- tsum_test(
    m1 = m1, sd1 = sd1, n1 = n1,
    m2 = m2, sd2 = sd2, n2 = n2,
    paired = paired,
    var.equal = var.equal,
    alternative = alt_low,
    mu = low_eqbound,
    conf.level = 1-alpha*2
  )

  high_ttest <- tsum_test(
    m1 = m1, sd1 = sd1, n1 = n1,
    m2 = m2, sd2 = sd2, n2 = n2,
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
    TOST_restext = paste0("The equivalence test was ",TOSToutcome,", t(",round(tresult$parameter, digits=2),") = ",format(tTOST, digits = 3, nsmall = 3, scientific = FALSE),", p = ",format(pTOST, digits = 3, nsmall = 3, scientific = FALSE),sep="")
  } else {
    TOST_restext = paste0("The minimal effect test was ",TOSToutcome,", t(",round(tresult$parameter, digits=2),") = ",format(tTOST, digits = 3, nsmall = 3, scientific = FALSE),", p = ",format(pTOST, digits = 3, nsmall = 3, scientific = FALSE),sep="")
  }

  ttest_restext = paste0("The null hypothesis test was ",testoutcome,", t(",round(tresult$parameter, digits=2),") = ",format(tresult$statistic, digits = 3, nsmall = 3, scientific = FALSE),", p = ",format(tresult$p.value, digits = 3, nsmall = 3, scientific = FALSE),sep="")
  if (hypothesis == "EQU"){
    if(tresult$p.value <= alpha && pTOST <= alpha){
      combined_outcome <- paste0("statistically different from ",mu_text," but statistically equivalent")
    }
    if(tresult$p.value < alpha && pTOST > alpha){
      combined_outcome <- paste0("statistically different from ",mu_text," and not statistically equivalent")
    }
    if(tresult$p.value > alpha && pTOST <= alpha){
      combined_outcome <- paste0("statistically not different from ",mu_text," and statistically equivalent")
    }
    if(tresult$p.value > alpha && pTOST > alpha){
      combined_outcome <- paste0("statistically not different from ",mu_text," and not statistically equivalent")
    }
  } else {
    if(tresult$p.value <= alpha && pTOST <= alpha){
      combined_outcome <- paste0("statistically different from ",mu_text," and statistically greater than the minimal effect threshold")
    }
    if(tresult$p.value < alpha && pTOST > alpha){
      combined_outcome <- paste0("statistically different from ",mu_text," but not statistically greater than the minimal effect threshold")
    }
    if(tresult$p.value > alpha && pTOST <= alpha){
      combined_outcome <- paste0("statistically not different from ",mu_text," and statistically greater than the minimal effect threshold")
    }
    if(tresult$p.value > alpha && pTOST > alpha){
      combined_outcome <- paste0("statistically not different from ",mu_text," and not statistically greater than the minimal effect threshold")
    }
  }

  decision = list(
    TOST = TOST_restext,
    ttest = ttest_restext,
    combined = combined_outcome
  )

  rval = list(
    TOST = TOST,
    eqb = eqb,
    alpha = alpha,
    method = tresult$method,
    hypothesis = test_hypothesis,
    effsize = effsize,
    smd = cohen_res,
    decision = decision
  )

  class(rval) = "TOSTt"

  return(rval)
}