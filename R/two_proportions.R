#' Testing for differences in proportions in 2 groups
#'
#' #' @description
#' `r lifecycle::badge('maturing')`
#'
#' This is a hypothesis testing function that mimics [prop.test], but focuses only on testing differences in proportions between two groups.
#' This function utilizes a z-test with unpooled variance estimates.
#'
#' @param p1,p2 Proportions in each respective group.
#' @param n1,n2 sample size in each respective group.
#' @param null a number indicating the null hypothesis of the difference in proportions between two groups.
#' @param effect_size the effect size estimate, and confidence intervals, to calculate. Options include the difference between both proportions ("difference"), odds ratio ("odds.ratio"), or risk ratio ("risk.ratio").
#' @inheritParams z_cor_test
#' @return An S3 object of the class `htest`.
#'
#' @references
#' Tunes da Silva, G., Logan, B. R., & Klein, J. P. (2008). Methods for Equivalence and Noninferiority Testing. Biology of Blood Marrow Transplant, 15(1 Suppl), 120-127.
#'
#' Yin, G. (2012). Clinical Trial Design: Bayesian and Frequentist Adaptive Methods. Hoboken, New Jersey: John Wiley & Sons, Inc.
#' @export

twoprop_test = function(p1, p2,
                        n1, n2,
                        null,
                        alpha = .05,
                        alternative = c("two.sided",
                                        "less",
                                        "greater",
                                        "equivalence",
                                        "minimal.effect"),
                        effect_size = c("difference",
                                        "odds.ratio",
                                        "risk.ratio")){

  if(alternative %in% c("equivalence","minimal.effect")){

    conf_level = 1-alpha*2
    conf = 1-alpha

  } else{

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


  METHOD = "difference in two proportions z-test"
  RVAL <- list(statistic = res_tests$STATISTIC,
               #parameter = PARAMETER,
               p.value = as.numeric(res_tests$PVAL),
               estimate = res_tests$ESTIMATE,
               null.value = res_tests$NVAL,
               conf.int = res_tests$CINT,
               alternative = alternative,
               method = METHOD)
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
  YATES <- abs(prop_dif) / sum(1 / n1, 1 / n2)
  if (any((null <= 0) | (null >= 1))) {
    stop("elements of 'null' must be in (0,1)")
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

  CINT = prop_dif - c(-1,1)*(qnorm(conf) * prop_se)
  if(alternative == "less"){
    CINT[1] = -Inf
  }

  if(alternative == "greater"){
    CINT[2] = Inf
  }
  attr(CINT, "conf.level") <- conf_level

  STATISTIC = ztest
  names(STATISTIC) <- "z"

  NVAL = null
  names(NVAL) = rep("difference in proportions", length(null))

  list(STATISTIC = STATISTIC,
       PVAL = PVAL,
       NVAL = NVAL,
       ESTIMATE = ESTIMATE,
       CINT = CINT)
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

  if(alternative %in% c("equivalence","minimal.effect")){

    if (length(null) == 1) {
      if (null ==  1) {
        stop("null for odds ratio cannot be 1 if alternative is equivalence or minimal.effect")
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
  se_logodds = sqrt(sum(a,b,c,d))

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

  list(ESTIMATE = ESTIMATE,
       CINT = CINT)

}

miettinen_nurminen_z = function(p1, p2,
                                n1, n2,
                                m1, null){

  a_parameter = n2*(null-1)
  b_parameter = n1*null+n2-m1*(null-1)
  c_parameter = -1*m1


  p_tilde2 = (-1*b_parameter + sqrt(b_parameter^2-4*a_parameter*c_parameter) ) / (2*a_parameter)

  p_tilde1 = (p_tilde2 * null)/ (1+ p_tild2*(null-1))
}

test_risk_ratio = function(p1, p2, n1, n2,
                           null,
                           alternative,
                           alpha){
  # Gart and Nam (1988), page 324
  # Gart, John J. and Nam, Jun-mo. 1988. 'Approximate Interval Estimation of the Ratio of Binomial Parameters: Review and Corrections for Skewness.' Biometrics, Volume 44, 323-338
  phi = p1/p2
  z_mult = qnorm(conf)
  q1 = 1-p1
  q2 = 1-p2

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

  se_val = sqrt(q1/(n1*p1) + q2/(n2*p2))
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

  list(ESTIMATE = ESTIMATE,
       CINT = CINT)
}




