#' @title Equivalence Test for ANOVA Results
#' @description Performs equivalence test on the partial eta-squared (pes) value from ANOVA results.
#'@param object an object of returned by either \code{Anova},
#'  \code{aov}, or \code{afex_aov}
#' @param eqbound Equivalence bound for the partial eta-squared.
#' @param MET logical indicator to perform a minimal effect test rather than equivalence test (default is FALSE).
#' @param alpha alpha used for the test (e.g., 0.05).
#' @details For details on the calculations in this function see vignette("the_ftestTOSTER").
#'@return Returns a data frame containing the ANOVA results with equivalence tests added.
#'
#'  The following abbreviations are used in the table:
#'
#'  \itemize{ \item effect: name of the effect. \item df1:	Degrees of Freedom in the numerator (i.e. DF effect).
#'  \item df2:	Degrees of Freedom in the denominator (i.e., DF error). \item F:	F-value. \item p.null:	p-value
#'  (probability of the data given the null hypothesis).  \item pes:	partial
#'  Eta-Squared measure of effect size. \item eqbound: equivalence bound. \item p.equ: p-value
#'  (probability of the data given the equivalence hypothesis) }
#'
#' @section References:
#' Campbell, H., & Lakens, D. (2021). Can we disregard the whole model? Omnibus non‐inferiority testing for R2 in multi‐variable linear regression and in ANOVA. British Journal of Mathematical and Statistical Psychology, 74(1), 64-89. doi: 10.1111/bmsp.12201
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
