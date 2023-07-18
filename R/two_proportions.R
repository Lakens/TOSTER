

twoprop_test = function(p1, p2,
                        n1, n2,
                        null = 0,
                        alpha = .05,
                        alternative = c("two.sided",
                                        "less",
                                        "greater",
                                        "equivalence",
                                        "minimal.effect")){

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

  # NORMAL APPROX.
  conf_int = prop_dif - c(-1,1)*(qnorm(conf) * prop_se)

}
