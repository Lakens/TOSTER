#' Test of Proportions between 2 Independent Groups
#'
#' @description
#' `r lifecycle::badge('maturing')`
#'
#' This is a hypothesis testing function that mimics [prop.test], but focuses only on testing differences in proportions between two groups.
#' This function utilizes a z-test to calculate the p-values (may be inaccurate with small sample sizes).
#'
#' @param p1,p2 Proportions in each respective group.
#' @param n1,n2 sample size in each respective group.
#' @param null a number indicating the null hypothesis of the difference in proportions between two groups.
#' @param effect_size the effect size estimate, and confidence intervals, to calculate. Options include the difference between both proportions ("difference"), odds ratio ("odds.ratio"), or risk ratio ("risk.ratio").
#' @inheritParams z_cor_test
#' @return An S3 object of the class `htest`.
#'
#' @details
#' The hypothesis test for differences in proportions can be made on the raw proportions scale, the odds ratio, or the risk ratio (details below).
#' This function uses the large sample size asymptotic approximations for both the p-value and confidence interval calculations.
#' There should be a good deal of caution when sample sizes are small.
#' The p-values for the differences in proportions will differ from base [prop.test] due to the use of the unpooled standard error (see below).
#'
#' ## Differences in Proportions
#'
#' Differences in proportions test is based on the following calculation:
#'
#' \deqn{d = p_1 - p_2}
#'
#' The standard error of \eqn{d} is calculated as the following:
#'
#' \deqn{se(d) = \sqrt{\frac{p_1 \cdot (1-p_1)}{n_1} + \frac{p_2 \cdot (1-p_2)}{n_2}} }
#'
#' The z-test, with \eqn{d_0} being the null value, is then calculated as the following (standard normal distribution evaluated to calculate p-value):
#'
#' \deqn{z = \frac{d - d_0}{se(d)}}
#'
#' The confidence interval can then be calculated as the following:
#'
#' \deqn{d_{lower},d_{upper} = d \pm z_{\alpha} \cdot se(d)}
#'
#' ## Risk Ratio
#'
#' The ratio between proportions test is based on the following calculation:
#'
#' \deqn{\phi = p_1/p_2}
#'
#' The standard error of \eqn{ln(\phi)} is calculated as the following:
#'
#' \deqn{se(ln(\phi)) = \sqrt{\frac{1-p_1}{n_1 \cdot p_1} + \frac{1-p_2}{n_2 \cdot p_2}} }
#'
#' The z-test, with \eqn{\phi_0} being the null value, is then calculated as the following (standard normal distribution evaluated to calculate p-value):
#'
#' \deqn{z = \frac{ln(\phi) - ln(\phi_0)}{se(ln(\phi))}}
#'
#' The confidence interval can then be calculated as the following:
#'
#'
#' \deqn{\phi_{lower} = \phi \cdot e^{-z_{\alpha} \cdot se(ln(\phi))}}
#'
#' \deqn{\phi_{upper} = \phi \cdot e^{z_{\alpha} \cdot se(ln(\phi))}}
#'
#' ## Odds Ratio
#'
#' The ratio between proportions test is based on the following calculation:
#' (p1/q1) / (p2/q2)
#'
#' \deqn{OR = \frac{p_1}{1-p_1} / \frac{p_2}{1-p_2}}
#'
#' The standard error of \eqn{ln(OR)} is calculated as the following:
#'
#'
#' \deqn{se(ln(OR)) = \sqrt{\frac{1}{n_1 \cdot p_1 + 0.5} + \frac{1}{n_1 \cdot (1-p_1) + 0.5} + \frac{1}{n_2 \cdot p_2 + 0.5} + \frac{1}{n_2 \cdot (1-p_2) + 0.5} } }
#'
#' The z-test, with \eqn{OR_0} being the null value, is then calculated as the following (standard normal distribution evaluated to calculate p-value):
#'
#' \deqn{z = \frac{ln(OR) - ln(OR_0)}{se(ln(OR))}}
#'
#' The confidence interval can then be calculated as the following:
#'
#' \deqn{OR_{lower},OR_{upper} = exp(ln(OR) \pm z_{\alpha} \cdot se(ln(OR)))}
#'
#'
#' @references
#'
#' Gart, J. J., & Nam, J. M. (1988). Approximate interval estimation of the ratio of binomial parameters: a review and corrections for skewness. Biometrics, 323-338.
#'
#' Tunes da Silva, G., Logan, B. R., & Klein, J. P. (2008). Methods for Equivalence and Noninferiority Testing. Biology of Blood Marrow Transplant, 15(1 Suppl), 120-127.
#'
#' Yin, G. (2012). Clinical Trial Design: Bayesian and Frequentist Adaptive Methods. Hoboken, New Jersey: John Wiley & Sons, Inc.
#' @export

twoprop_test = function(p1, p2,
                        n1, n2,
                        null = NULL,
                        alpha = .05,
                        alternative = c("two.sided",
                                        "less",
                                        "greater",
                                        "equivalence",
                                        "minimal.effect"),
                        effect_size = c("difference",
                                        "odds.ratio",
                                        "risk.ratio")){

  alternative = match.arg(alternative)
  effect_size = match.arg(effect_size)
  if(!is.numeric(alpha) || alpha <=0 || alpha >=1){
    stop("The alpha must be a numeric value between 0 and 1")
  }
  if (any(c(p1, p2) > 1 | c(p1, p2) < 0)
  ) {
    stop("elements of p1 or p2 must be in [0,1]")
  }
  if(any(c(n1, n2) < 10)){
    stop("elements of n1 or n2 must be greater than 9")
  }

  if(any(c(n1, n2) <= 50)){
    message("Small sample size in at least one group; proceed with caution")
  }

  # effect size ----
  res_tests = switch(effect_size,
                       difference = test_prop_dif(p1, p2, n1, n2,
                                                  null,
                                                  alternative,
                                                  alpha),
                       odds.ratio = test_odds_ratio(p1, p2, n1, n2,
                                                    null,
                                                    alternative,
                                                    alpha),
                       risk.ratio = test_risk_ratio(p1, p2, n1, n2,
                                                    null,
                                                    alternative,
                                                    alpha))



  RVAL <- list(statistic = res_tests$STATISTIC,
               #parameter = PARAMETER,
               p.value = as.numeric(res_tests$PVAL),
               estimate = res_tests$ESTIMATE,
               null.value = res_tests$NVAL,
               conf.int = res_tests$CINT,
               alternative = alternative,
               method = res_tests$METHOD)
  class(RVAL) <- "htest"
  return(RVAL)

}

test_prop_dif = function(p1,p2,n1,n2,
                         null,
                         alternative,
                         alpha) {
  if(is.null(null)){
    null = 0
  }

  prop_dif <- p1 - p2
  # proportion se
  prop_se <- sqrt((p1*(1-p1))/n1 + (p2*(1-p2))/n2)
  p_bar = (n1*p1 + n2*p2) / (n1 + n2)
  prop_se_pool = sqrt(p_bar * (1-p_bar) * (1/n1+1/n2))

  #YATES <- abs(prop_dif) / sum(1 / n1, 1 / n2)
  if (any((null <= -1) | (null >= 1))) {
    stop("elements of 'null' must be in (-1,1)")
  }

  if(alternative %in% c("equivalence","minimal.effect")){

    if (length(null) == 1) {
      if (null ==  0) {
        stop("null cannot be zero if alternative is equivalence or minimal.effect")
      }
      null = c(null,-1 * null)

    }

    lo_ztest = (prop_dif - min(null))/prop_se
    hi_ztest = (prop_dif - max(null))/prop_se

    lo_pvalue = switch(
      alternative,
      "equivalence" = p_from_z(lo_ztest,
                               alternative = "greater"),
      "minimal.effect" = p_from_z(lo_ztest,
                                  alternative = "less")
    )

    hi_pvalue = switch(
      alternative,
      "equivalence" = p_from_z(hi_ztest,
                               alternative = "less"),
      "minimal.effect" = p_from_z(hi_ztest,
                                  alternative = "greater")
    )

    PVAL = switch(
      alternative,
      "equivalence" = max(lo_pvalue, hi_pvalue),
      "minimal.effect" = min(lo_pvalue, hi_pvalue)
    )

    test_p = PVAL == c(lo_pvalue, hi_pvalue)
    test_z = c(lo_ztest, hi_ztest)
    ZTEST = test_z[test_p]

    conf_level = 1-alpha*2
    conf = 1-alpha
  } else{
    if(length(null) != 1){
      stop("null must have length of 1 if alternative is not a TOST test.")
    }

    if(null == 0){
      message("For nil-hypothesis tests (null = 0), it is recommended that prop.test be utilized.")
    }
    ZTEST = (prop_dif - null) / prop_se
    PVAL = p_from_z(ZTEST,
                    alternative = alternative)

    conf_level = switch(
      alternative,
      "two.sided" = 1-alpha,
      "less" =  1-alpha*2,
      "greater" = 1-alpha*2
    )

    conf = switch(
      alternative,
      "two.sided" = 1-alpha/2,
      "less" =  1-alpha,
      "greater" = 1-alpha
    )
  }

  z_mult = qnorm(conf)
  ESTIMATE = prop_dif
  names(ESTIMATE) = "difference in proportions"

  CINT = prop_dif + c(-1,1)*(qnorm(conf) * prop_se)
  if(alternative == "less"){
    CINT[1] = -Inf
  }

  if(alternative == "greater"){
    CINT[2] = Inf
  }
  attr(CINT, "conf.level") <- conf_level

  STATISTIC = ZTEST
  names(STATISTIC) <- "z"

  NVAL = null
  names(NVAL) = rep("difference in proportions", length(null))
  METHOD = "difference in two proportions z-test"

  list(STATISTIC = STATISTIC,
       PVAL = PVAL,
       NVAL = NVAL,
       ESTIMATE = ESTIMATE,
       CINT = CINT,
       METHOD = METHOD)
}



test_odds_ratio = function(p1, p2, n1, n2,
                           null,
                           alternative,
                           alpha){
  if(is.null(null)){
    null = 1
  }
  # Fleiss, J. L., Levin, B., Paik, M.C. 2003. Statistical Methods for Rates and Proportions. Third Edition. John Wiley & Sons. New York.
  q1 = 1-p1
  q2 = 1-p2
  a = 1/(n1*p1+.5)
  b = 1/(n1*q1+.5)
  c = 1/(n2*p2+.5)
  d = 1/(n2*q2+.5)
  m1 = n1*p1 + n2*p2
  OR = (p1/q1) / (p2/q2)
  se_logodds = sqrt(sum(a,b,c,d))

  if(alternative %in% c("equivalence","minimal.effect")){

    if (length(null) == 1) {
      if (null ==  1) {
        stop("null for odds ratio cannot be 1 if alternative is equivalence or minimal.effect")
      }
      null = c(null,null^(-1))

    }

    lo_ztest = (log(OR) - min(log(null))) / se_logodds
    hi_ztest = (log(OR) - max(log(null))) / se_logodds

    lo_pvalue = switch(
      alternative,
      "equivalence" = p_from_z(lo_ztest,
                               alternative = "greater"),
      "minimal.effect" = p_from_z(lo_ztest,
                                  alternative = "less")
    )

    hi_pvalue = switch(
      alternative,
      "equivalence" = p_from_z(hi_ztest,
                               alternative = "less"),
      "minimal.effect" = p_from_z(hi_ztest,
                                  alternative = "greater")
    )

    PVAL = switch(
      alternative,
      "equivalence" = max(lo_pvalue, hi_pvalue),
      "minimal.effect" = min(lo_pvalue, hi_pvalue)
    )

    test_p = PVAL == c(lo_pvalue, hi_pvalue)
    test_z = c(lo_ztest, hi_ztest)
    ZTEST = test_z[test_p]

    conf_level = 1-alpha*2
    conf = 1-alpha
  } else{
    if(length(null) != 1){
      stop("null must have length of 1 if alternative is not a TOST test.")
    }
    if(null == 1){
      message("For nil-hypothesis tests (null = 1), it is recommended that prop.test be utilized.")
    }
    ZTEST = (log(OR) - log(null)) / se_logodds
    PVAL = p_from_z(ZTEST,
                    alternative = alternative)

    conf_level = switch(
      alternative,
      "two.sided" = 1-alpha,
      "less" =  1-alpha*2,
      "greater" = 1-alpha*2
    )

    conf = switch(
      alternative,
      "two.sided" = 1-alpha/2,
      "less" =  1-alpha,
      "greater" = 1-alpha
    )
  }

  z_mult = qnorm(conf)


  CINT = exp(
    log(OR) + c(-1,1)*z_mult*se_logodds
    )

  if(alternative == "less"){
    CINT[1] = -Inf
  }

  if(alternative == "greater"){
    CINT[2] = Inf
  }
  attr(CINT, "conf.level") <- conf_level

  ESTIMATE = OR
  names(ESTIMATE) = "Odds Ratio"

  STATISTIC = ZTEST
  names(STATISTIC) <- "z"

  NVAL = null
  names(NVAL) = rep("Odds Ratio", length(null))
  METHOD = "approximate Odds Ratio z-test"
  list(STATISTIC = STATISTIC,
       PVAL = PVAL,
       NVAL = NVAL,
       ESTIMATE = ESTIMATE,
       CINT = CINT,
       METHOD = METHOD)

}


test_risk_ratio = function(p1, p2, n1, n2,
                           null,
                           alternative,
                           alpha){
  if(is.null(null)){
    null = 1
  }
  # Gart and Nam (1988), page 324
  # Gart, John J. and Nam, Jun-mo. 1988. 'Approximate Interval Estimation of the Ratio of Binomial Parameters: Review and Corrections for Skewness.' Biometrics, Volume 44, 323-338
  phi = p1/p2
  #
  q1 = 1-p1
  q2 = 1-p2
  se_val = sqrt(q1/(n1*p1) + q2/(n2*p2))

  if(alternative %in% c("equivalence","minimal.effect")){

    if (length(null) == 1) {
      if (null ==  1) {
        stop("null cannot be zero if alternative is equivalence or minimal.effect")
      }
      null = c(null,null^(-1))

    }

    lo_ztest = (log(phi) - min(log(null))) / se_val
    hi_ztest = (log(phi) - max(log(null))) / se_val

    lo_pvalue = switch(
      alternative,
      "equivalence" = p_from_z(lo_ztest,
                               alternative = "greater"),
      "minimal.effect" = p_from_z(lo_ztest,
                                  alternative = "less")
    )

    hi_pvalue = switch(
      alternative,
      "equivalence" = p_from_z(hi_ztest,
                               alternative = "less"),
      "minimal.effect" = p_from_z(hi_ztest,
                                  alternative = "greater")
    )

    PVAL = switch(
      alternative,
      "equivalence" = max(lo_pvalue, hi_pvalue),
      "minimal.effect" = min(lo_pvalue, hi_pvalue)
    )

    test_p = PVAL == c(lo_pvalue, hi_pvalue)
    test_z = c(lo_ztest, hi_ztest)
    ZTEST = test_z[test_p]

    conf_level = 1-alpha*2
    conf = 1-alpha
  } else{
    if(length(null) != 1){
      stop("null must have length of 1 if alternative is not a TOST test.")
    }

    if(null == 1){
      message("For nil-hypothesis tests (null = 1), it is recommended that prop.test be utilized.")
    }
    ZTEST = (log(phi) - log(null)) / se_val
    PVAL = p_from_z(ZTEST,
                    alternative = alternative)

    conf_level = switch(
      alternative,
      "two.sided" = 1-alpha,
      "less" =  1-alpha*2,
      "greater" = 1-alpha*2
    )

    conf = switch(
      alternative,
      "two.sided" = 1-alpha/2,
      "less" =  1-alpha,
      "greater" = 1-alpha
    )
  }

  z_mult = qnorm(conf)
  CINT = phi * exp(c(-1,1)*z_mult*se_val)

  if(alternative == "less"){
    CINT[1] = -Inf
  }

  if(alternative == "greater"){
    CINT[2] = Inf
  }
  attr(CINT, "conf.level") <- conf_level

  ESTIMATE = phi
  names(ESTIMATE) = "Risk Ratio"

  STATISTIC = ZTEST
  names(STATISTIC) <- "z"

  NVAL = null
  names(NVAL) = rep("Risk Ratio", length(null))

  METHOD = "approximate Risk Ratio z-test"

  list(STATISTIC = STATISTIC,
       PVAL = PVAL,
       NVAL = NVAL,
       ESTIMATE = ESTIMATE,
       CINT = CINT,
       METHOD = METHOD)
}




