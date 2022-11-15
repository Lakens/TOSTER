log_pair = function(x,
                    hypothesis = c("EQU", "MET"),
                    #var.equal = FALSE,
                    #paired = FALSE,
                    eqb = 1.25,
                    alpha = 0.05,
                    null = 1) {
  hypothesis = match.arg(hypothesis)

    sample_type = "Paired Sample"



  if(hypothesis == "EQU"){
    alt_low = "greater"
    alt_high = "less"
    test_hypothesis = "Hypothesis Tested: Equivalence"

  } else if(hypothesis == "MET"){
    alt_low = "less"
    alt_high = "greater"
    test_hypothesis = "Hypothesis Tested: Minimal Effect"

  }


  tresult = t.test(x = x,
                   mu = log(null),
                   conf.level = 1 - alpha*2,
                   alternative = "two.sided")
  logrom = (tresult$statistic) * tresult$stderr + log(null)
  logSE = tresult$stderr

  n <- length(x)

  d_df = n - 1

  rom_res = list(
    d = exp(logrom),
    d_df = d_df,
    dlow = exp(tresult$conf.int[1]),
    dhigh = exp(tresult$conf.int[2]),
    d_sigma = logSE,
    d_lambda = NULL,
    #hn = hn,
    smd_label = "Means Ratio",
    J = NULL,
    d_denom = 1,
    ntilde = 1,
    t_stat = NULL,
    smd_ci = "t"
  )

  if(length(eqb) == 1){
    if(eqb > 1){
      high_eqbound = eqb
      low_eqbound = 1/eqb
    } else{
      high_eqbound = 1/eqb
      low_eqbound = eqb
    }

  } else {
    high_eqbound = max(eqb)
    low_eqbound = min(eqb)
  }


  if(hypothesis == "EQU"){
    null_hyp = paste0(round(low_eqbound,2),
                      " >= (Mean1/Mean2) or (Mean1/Mean2) >= ",
                      round(high_eqbound,2))
    alt_hyp = paste0(round(low_eqbound,2),
                     " < (Mean1/Mean2) < ",
                     round(high_eqbound,2))
  } else if(hypothesis == "MET"){
    null_hyp = paste0(round(low_eqbound,2),
                      " <= (Mean1/Mean2)  <= ",
                      round(high_eqbound,2))
    alt_hyp = paste0(round(low_eqbound,2),
                     " > (Mean1/Mean2) or (Mean1/Mean2)  > ",
                     round(high_eqbound,2))
  }

  low_ttest <- t.test(
    x = x,
    alternative = alt_low,
    mu = log(low_eqbound),
    conf.level = 1-alpha*2
  )

  high_ttest <- t.test(
    x = x,
    alternative = alt_high,
    mu = log(high_eqbound),
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
    type = c("log(Means Ratio)","Means Ratio"),
    low_eq = c(log(low_eqbound),low_eqbound),
    high_eq = c(log(high_eqbound),high_eqbound)
  )

  effsize = data.frame(
    estimate = c(logrom,
                 exp(logrom)),
    SE = c(logSE,NA),
    lower.ci = c(tresult$conf.int[1], exp(tresult$conf.int[1])),
    upper.ci = c(tresult$conf.int[2], exp(tresult$conf.int[2])),
    conf.level = c((1-alpha*2),(1-alpha*2)),
    row.names = c("log(Means Ratio)","Means Ratio")
  )
  TOSToutcome<-ifelse(pTOST<alpha,"significant","non-significant")
  testoutcome<-ifelse(tresult$p.value<alpha,"significant","non-significant")

  # Change text based on two tailed t test if mu is not zero
  if(null == 1){
    mu_text = "one"
  } else {
    mu_text = null
  }

  if(hypothesis == "EQU"){
    #format(low_eqbound, digits = 3, nsmall = 3, scientific = FALSE)
    TOST_restext = paste0("The equivalence test was ",TOSToutcome,", t(",round(tresult$parameter, digits=2),") = ",format(tTOST, digits = 3, nsmall = 3, scientific = FALSE),", p = ",format(pTOST, digits = 3, nsmall = 3, scientific = TRUE),sep="")
  } else {
    TOST_restext = paste0("The minimal effect test was ",TOSToutcome,", t(",round(tresult$parameter, digits=2),") = ",format(tTOST, digits = 3, nsmall = 3, scientific = FALSE),", p = ",format(pTOST, digits = 3, nsmall = 3, scientific = TRUE),sep="")
  }

  ttest_restext = paste0("The null hypothesis test was ",testoutcome,", t(",round(tresult$parameter, digits=2),") = ",format(tresult$statistic, digits = 3, nsmall = 3, scientific = FALSE),", p = ",format(tresult$p.value, digits = 3, nsmall = 3, scientific = TRUE),sep="")
  if (hypothesis == "EQU"){
    if(tresult$p.value <= alpha && pTOST <= alpha){
      combined_outcome <- paste0("NHST: reject null significance hypothesis that the effect is equal to ",mu_text," \n",
                                 "TOST: reject null equivalence hypothesis")
    }
    if(tresult$p.value < alpha && pTOST > alpha){
      combined_outcome <- paste0("NHST: reject null significance hypothesis that the effect is equal to ",mu_text," \n ",
                                 "TOST: don't reject null equivalence hypothesis")
      # paste0("statistically different from ",mu_text," and not statistically equivalent")
    }
    if(tresult$p.value > alpha && pTOST <= alpha){
      combined_outcome <- paste0("NHST: don't reject null significance hypothesis that the effect is equal to ",mu_text," \n ",
                                 "TOST: reject null equivalence hypothesis")
      #paste0("statistically not different from ",mu_text," and statistically equivalent")
    }
    if(tresult$p.value > alpha && pTOST > alpha){
      combined_outcome <- paste0("NHST: don't reject null significance hypothesis that the effect is equal to ",mu_text," \n ",
                                 "TOST: don't reject null equivalence hypothesis")
      #paste0("statistically not different from ",mu_text," and not statistically equivalent")
    }
  } else {
    if(tresult$p.value <= alpha && pTOST <= alpha){
      combined_outcome <- paste0("NHST: reject null significance hypothesis that the effect is equal to ",mu_text," \n ",
                                 "TOST: reject null MET hypothesis")
      #paste0("statistically different from ",mu_text," and statistically greater than the minimal effect threshold")
    }
    if(tresult$p.value < alpha && pTOST > alpha){
      combined_outcome <- paste0("NHST: reject null significance hypothesis that the effect is equal to ",mu_text," \n ",
                                 "TOST: don't reject null MET hypothesis")
      #paste0("statistically different from ",mu_text," but not statistically greater than the minimal effect threshold")
    }
    if(tresult$p.value > alpha && pTOST <= alpha){
      combined_outcome <- paste0("NHST: don't reject null significance hypothesis that the effect is equal to ",mu_text," \n ",
                                 "TOST: reject null MET hypothesis")
      #paste0("statistically not different from ",mu_text," and statistically greater than the minimal effect threshold")
    }
    if(tresult$p.value > alpha && pTOST > alpha){
      combined_outcome <- paste0("NHST: don't reject null significance hypothesis that the effect is equal to ",mu_text," \n ",
                                 "TOST: don't reject null MET hypothesis")
      #paste0("statistically not different from ",mu_text," and not statistically greater than the minimal effect threshold")
    }
  }

  decision = list(
    TOST = TOST_restext,
    ttest = ttest_restext,
    combined = combined_outcome
  )

  #message(cat("Based on the equivalence test and the null-hypothesis test combined, we can conclude that the observed effect is ",combined_outcome,".",sep=""))


  rval = list(
    TOST = TOST,
    eqb = eqb,
    alpha = alpha,
    method = paste0("Log-transformed ",tresult$method),
    hypothesis = test_hypothesis,
    effsize = effsize,
    smd = rom_res,
    decision = decision,
    call = match.call()
  )

  class(rval) = "TOSTt"

  return(rval)

}
