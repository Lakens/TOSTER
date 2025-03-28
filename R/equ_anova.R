#' @title Equivalence Test for ANOVA Results
#' @description
#' `r lifecycle::badge('stable')`
#'
#' Performs equivalence or minimal effect testing on the partial eta-squared (pes) value
#' from ANOVA results to determine if effects are practically equivalent to zero or
#' meaningfully different from zero.
#'
#' @param object An object returned by either `Anova`, `aov`, or `afex_aov`.
#' @param eqbound Equivalence bound for the partial eta-squared. This value represents
#'   the smallest effect size considered meaningful or practically significant.
#' @param MET Logical indicator to perform a minimal effect test rather than equivalence
#'   test (default is FALSE). When TRUE, the alternative hypothesis becomes that the effect
#'   is larger than the equivalence bound.
#' @param alpha Alpha level used for the test (default = 0.05).
#'
#' @details
#' This function tests whether ANOVA effects are practically equivalent to zero (when
#' `MET = FALSE`) or meaningfully different from zero (when `MET = TRUE`) using the approach
#' described by Campbell & Lakens (2021).
#'
#' The function works by:
#'
#' 1. Extracting ANOVA results from the input object
#' 2. Converting the equivalence bound for partial eta-squared to a non-centrality parameter
#' 3. Performing an equivalence test or minimal effect test for each effect in the ANOVA
#'
#' For equivalence tests (`MET = FALSE`), a significant result (p < alpha) indicates that the
#' effect is statistically equivalent to zero (smaller than the equivalence bound).
#'
#' For minimal effect tests (`MET = TRUE`), a significant result (p < alpha) indicates that
#' the effect is meaningfully different from zero (larger than the equivalence bound).
#'
#' For details on the calculations in this function see `vignette("the_ftestTOSTER")`.
#'
#' @return
#' Returns a data frame containing the ANOVA results with equivalence tests added.
#' The following columns are included in the table:
#'
#' * **effect**: Name of the effect.
#' * **df1**: Degrees of Freedom in the numerator (i.e., DF effect).
#' * **df2**: Degrees of Freedom in the denominator (i.e., DF error).
#' * **F.value**: F-value.
#' * **p.null**: p-value for the traditional null hypothesis test (probability of the data given the null hypothesis).
#' * **pes**: Partial eta-squared measure of effect size.
#' * **eqbound**: Equivalence bound used for testing.
#' * **p.equ**: p-value for the equivalence or minimal effect test.
#'
#' @examples
#' # One-way ANOVA
#' data(iris)
#' anova_result <- aov(Sepal.Length ~ Species, data = iris)
#'
#' # Equivalence test with bound of 0.1
#' equ_anova(anova_result, eqbound = 0.1)
#'
#' # Minimal effect test with bound of 0.1
#' equ_anova(anova_result, eqbound = 0.1, MET = TRUE)
#'
#' # Two-way ANOVA with lower equivalence bound
#' anova_result2 <- aov(Sepal.Length ~ Species * Petal.Width, data = iris)
#' equ_anova(anova_result2, eqbound = 0.05)
#'
#' @references
#' Campbell, H., & Lakens, D. (2021). Can we disregard the whole model? Omnibus non‐inferiority
#' testing for R2 in multi‐variable linear regression and in ANOVA. British Journal of
#' Mathematical and Statistical Psychology, 74(1), 64-89. doi: 10.1111/bmsp.12201
#'
#' @family f-test
#' @export


equ_anova <- function(object,
                      eqbound,
                      MET = FALSE,
                      alpha = 0.05){

  #message("Note: equ_anova only validated for one-way ANOVA; use with caution")

  if(inherits(object, "Anova.mlm")){
    results <- anova_summary(object)
  }
  else if(inherits(object, "anova")){
    results <- anova_summary(object)
  }
  else if(inherits(object, c("aov", "aovlist"))){
    results <- anova_summary(object)
  } else if (inherits(object, "afex_aov")){
    if(is.null(object$aov)){
      aov_res = object$Anova
    } else{
      aov_res = object$aov
    }
    results <- anova_summary(aov_res)
  } else{
    stop("Non-supported object passed: ",
         paste(class(object), collapse = ", "), ". ",
         "Object needs to be of class 'Anova.mlm', 'afex_aov', or 'anova'.")
  }

  res2 = results[c("Effect","df1","df2","F.value","p.value","pes")]
  colnames(res2) = c("effect","df1","df2","F.value","p.null","pes")
  res2$f2 = eqbound/(1 - eqbound)
  res2$lambda = (res2$f2 * (res2$df1 + res2$df2 + 1))

  res2$p.equ = pf(res2$F.value,
                  df1 = res2$df1,
                  df2 = res2$df2,
                  ncp = res2$lambda,
                  lower.tail = ifelse(MET,FALSE,TRUE))

  #res2$p.equ = suppressMessages({equ_ftest(
  #  Fstat = res2$F.value,
  #  df1 = res2$df1,
  #  df2 = res2$df2,
  #  eqbound = eqbound,
  #  MET = MET,
  #  alpha = alpha
  #)$p.value})

  res2$eqbound = eqbound

  res3 = res2[c("effect",
                "df1",
                "df2",
                "F.value",
                "p.null",
                "pes",
                "eqbound",
                "p.equ")]

  return(res3)


}
