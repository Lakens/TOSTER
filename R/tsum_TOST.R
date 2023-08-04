#' @title TOST with t-tests from Summary Statistics
#' @description
#' `r lifecycle::badge('stable')`
#'
#'  A function for TOST with all types of t-tests from summary statistics.
#' @param m1 mean of group 1.
#' @param m2 mean of group 2.
#' @param sd1 standard deviation of group 1.
#' @param sd2 standard deviation of group 2.
#' @param n1 sample size in group 1.
#' @param n2 sample size in group 2.
#' @param r12 correlation of dependent variable between group 1 and group 2.
#' @inheritParams t_TOST
#' @details
#' For details on the calculations in this function see
#' `vignette("IntroTOSTt")` & `vignette("SMD_calcs")`.
#'
#' For two-sample tests, the test is of \eqn{m1 - m2} (mean of 1 minus mean of 2).
#' For paired samples, the test is of the difference scores (z),
#' wherein \eqn{z =  m1 - m2}, and the test is of \eqn{\bar z} (mean of the difference scores).
#' For one-sample tests, the test is of \eqn{\bar m1 } (mean of group 1).
#'
#' @return An S3 object of class
#'  `"TOSTt"` is returned containing the following slots:
#'
#'   - "TOST":  A table of class `"data.frame"` containing two-tailed t-test and both one-tailed results.
#'   - "eqb":  A table of class `"data.frame"` containing equivalence bound settings.
#'   - "effsize":   table of class `"data.frame"` containing effect size estimates.
#'   - "hypothesis":  String stating the hypothesis being tested.
#'   - "smd":  List containing the results of the standardized mean difference calculations (e.g., Cohen's d).
#'      - Items include: d (estimate), dlow (lower CI bound), dhigh (upper CI bound), d_df (degrees of freedom for SMD), d_sigma (SE), d_lambda (non-centrality), J (bias correction), smd_label (type of SMD), d_denom (denominator calculation)
#'   - "alpha":  Alpha level set for the analysis.
#'   - "method":  Type of t-test.
#'   - "decision":  List included text regarding the decisions for statistical inference.
#'
#' @examples
#' # example code
#' # One sample test
#' tsum_TOST(m1 = 0.55, n1 = 18, sd1 = 4, eqb  = 2)
#' @family TOST
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
                      hypothesis = c("EQU","MET"),
                      paired = FALSE,
                      var.equal = FALSE,
                      eqb,
                      low_eqbound,
                      high_eqbound,
                      mu = 0,
                      eqbound_type = c("raw","SMD"),
                      alpha = 0.05,
                      bias_correction = TRUE,
                      rm_correction = FALSE,
                      glass = NULL,
                      smd_ci = c("nct", "goulet", "t", "z")){

  hypothesis = match.arg(hypothesis)
  eqbound_type = match.arg(eqbound_type)

  if(is.null(glass)){
    glass = "no"
  }
  smd_ci = match.arg(smd_ci)


  if(bias_correction){
    smd_type = 'g'
  } else {
    smd_type = 'd'
  }

  if(is.null(n2) || is.null(m2) || is.null(sd2)){
    sample_type = "One Sample"
  } else if(paired == TRUE && !is.null(r12)) {
    sample_type = "Paired Sample"
  } else if (paired == TRUE && is.null(r12)){
    stop("paired == TRUE but r12 not provided. Must provide correlation.")
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

  if(!(paired)){
    r12 = NULL
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

  tresult = tsum_test(m1 = m1, sd1 = sd1, n1 = n1,
                      m2 = m2, sd2 = sd2, n2 = n2,
                      r12 = r12,
                      paired = paired,
                      var.equal = var.equal,
                      mu = mu,
                      conf.level = 1 - alpha * 2,
                      alternative = "two.sided")

  if(paired == TRUE && !missing(r12)){

    cohen_res = d_est_pair(
      n = n1,
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
      alpha = alpha,
      denom = denom,
      smd_ci = smd_ci
    )

  } else {
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
    r12 = r12,
    paired = paired,
    var.equal = var.equal,
    alternative = alt_low,
    mu = low_eqbound,
    conf.level = 1-alpha*2
  )

  high_ttest <- tsum_test(
    m1 = m1, sd1 = sd1, n1 = n1,
    m2 = m2, sd2 = sd2, n2 = n2,
    r12 = r12,
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
                          TOSToutcome,", t(",
                          round(tresult$parameter, digits=2),") = ",
                          format(tTOST, digits = 3,
                                 nsmall = 3,
                                 scientific = FALSE),", p = ",
                          format(pTOST, digits = 3,
                                 nsmall = 3, scientific = TRUE),sep="")
  } else {
    TOST_restext = paste0("The minimal effect test was ",
                          TOSToutcome,", t(",
                          round(tresult$parameter, digits=2),") = ",
                          format(tTOST, digits = 3,
                                 nsmall = 3, scientific = FALSE),", p = ",
                          format(pTOST,
                                 digits = 3, nsmall = 3,
                                 scientific = TRUE),sep="")
  }

  ttest_restext = paste0("The null hypothesis test was ",
                         testoutcome,", t(",round(tresult$parameter,
                                                  digits=2),") = ",
                         format(tresult$statistic, digits = 3,
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
