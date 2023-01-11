#' Methods for TOSTnp objects
#'
#' Methods defined for objects returned from the wilcox_TOST function.
#'
#' @param x object of class \code{TOSTnp}
#' @param digits Number of digits to print for p-values
#' @param ... further arguments passed through, see description of return value
#'   for details.
#'   \code{\link{TOSTnp-methods}}.
#' @return
#' \describe{
#'   \item{\code{print}}{Prints short summary of the tests.}
#'   \item{\code{describe}}{Verbose description of results.}
#' }
#'
#' @name TOSTnp-methods


### methods for TOSTnp objects

#' @rdname TOSTnp-methods
#' @method print TOSTnp
#' @export

print.TOSTnp <- function(x,
                        digits = 4,
                        ...){
  effsize = x$effsize
  TOST = x$TOST
  TOST$p.value = ifelse(TOST$p.value < 0.001,
                        "< 0.001",
                        round(TOST$p.value, 3))
  effsize$CI = paste0("[",
                      round(effsize$lower.ci,digits),
                      ", ",
                      round(effsize$upper.ci,digits),
                      "]")
  effsize = effsize[c("estimate", "CI","conf.level")]
  #TOST = TOST[c("t","df","p.value")]
  colnames(TOST) = c("Test Statistic", "p.value")
  colnames(effsize) = c("Estimate", "C.I.", "Conf. Level")
  cat("\n")
  cat(strwrap(x$method), sep = "\n")
  #cat(x$hypothesis, "\n", sep = "")
  #cat("Equivalence Bounds (raw):",format(x$eqb[1], digits = 3, nsmall = 3, scientific = FALSE)," & ",format(x$eqb[2], digits = 3, nsmall = 3, scientific = FALSE), sep="")
  #cat("\n")
  #cat("Alpha Level:", x$alpha, sep="")
  cat("\n")
  cat(x$decision$TOST)
  cat("\n")
  cat(x$decision$test)
  cat("\n")
  cat(x$decision$combined)
  cat("\n")
  cat("\n")
  cat("TOST Results \n")
  print(TOST, digits = digits)
  cat("\n")
  cat("Effect Sizes \n")
  print(effsize, digits = digits)
  cat("\n")
  if("boot" %in% names(x)){
    cat("Note: percentile boostrap method utilized.")
  }
  cat("\n")

}

#' @rdname TOSTnp-methods
#' @method summary TOSTnp
#' @export

summary.TOSTnp <- function(x,
                          digits = 3,
                          ...){

  tosty = x
  htest = as_htest(x)

  text = describe_htest(htest = htest,
                        alpha = tosty$alpha,
                        digits = digits)

  text2 = paste0(text,
                 " Additionally, a standardized effect size (SES), ",
                 row.names(x$effsize)[2],
                 ", was estimated, SES = ",
                 rounder_stat(x$effsize$estimate[2], digits = digits),
                 " ",
                 x$effsize$conf.level[2]*100,
                 "% C.I.[",
                 rounder_stat(x$effsize$estimate[2], digits = digits),
                 ", ",
                 rounder_stat(x$effsize$estimate[2],digits = digits),
                 "].")

  return(text2)
}

#' @rdname TOSTnp-methods
#' @export

describe <- function(x) {
  UseMethod("describe")
}

#' @rdname TOSTnp-methods
#' @method describe TOSTnp
#' @export

describe.TOSTnp <- function(x,
                           digits = 3,
                           ...){

  tosty = x
  htest = as_htest(x)

  type_tost = ifelse(htest$alternative == "equivalence",
                     "equivalence",
                     "minimal effect")
  nhst_null = ifelse(is.null(tosty$call$mu),
                     0,
                     tosty$call$mu)
  alt_nhst = paste0("true ",
                    names(htest$null.value[1]),
                    " is ",
                    "not equal to",
                    " ",
                    nhst_null)
  if(htest$alternative == "equivalence"){

    alt_tost = paste0("true ",
                      names(htest$null.value)[1],
                      " is ",
                      "between",
                      " ",
                      htest$null.value[1], " and ",
                      htest$null.value[2])

  }
  if(htest$alternative == "minimal.effect"){

    alt_tost = paste0("true ",
                      names(htest$null.value)[1],
                      " is ",
                      "less than ",
                      htest$null.value[1], " or ",
                      "greater than ",
                      htest$null.value[2])


  }

  method_state = paste0("Using the ", htest$method,
                        " a null hypothesis significance test (NHST)",
                        ", and two one-sided tests (TOST) were performed",
                        " with an alpha-level of ", x$alpha, ".",
                        " The NHST tested whether the ", alt_nhst,
                        ", and the TOST tested whether the ", alt_tost,".")

  smd_name = tosty$smd$smd_label

  pTOST = htest$p.value
  sigTOST = ifelse(pTOST < tosty$alpha, TRUE,FALSE)
  pNHST = tosty$TOST$p.value[1]
  sigNHST = ifelse(pNHST < tosty$alpha, TRUE,FALSE)

  if(sigTOST){
    sig_text = paste0("The ", type_tost, " test",
                      " was significant, ",
                      "WMW"," = ",
                      rounder_stat(htest$statistic, digits = digits),
                      ", ",
                      printable_pval(pTOST, digits = digits))

    claim_text = paste0("Therefore, we can claim that ",
                        alt_tost, ".")
  } else if(sigNHST){
    sig_text = paste0("The ", type_tost, " test was not significant (",
                      printable_pval(pTOST, digits = digits),").",
                      " The NHST",
                      " was significant, ",
                      "WMW", " = ",
                      rounder_stat(tosty$TOST$t[1], digits = digits),
                      ", ",
                      printable_pval(pNHST, digits = digits))

    claim_text = paste0("Therefore, we can claim that ",
                        alt_nhst,
                        " (i.e., no ",type_tost,
                        ").")
  } else {
    sig_text = paste0("Both the ", type_tost, " test (",
                      printable_pval(pTOST, digits = digits),
                      ") the NHST (",
                      printable_pval(pNHST, digits = digits),
                      ")",
                      " were not significant")

    claim_text = paste0("Therefore, the results are ambigious (no decision regarding the effect size/direction can be made).")
  }

  stat_text = paste0(sig_text,
                     " (",
                     names(htest$null.value[1]),
                     " = ",
                     rounder_stat(htest$estimate,
                                  digits = digits),
                     " ",
                     100 * attr(htest$conf.int, "conf.level"),
                     "% C.I.[",
                     rounder_stat(min(htest$conf.int),
                                  digits = digits),
                     ", ",
                     rounder_stat(max(htest$conf.int),
                                  digits = digits),
                     "]",
                     "; ",
                     smd_name,
                     " = ",
                     rounder_stat(tosty$effsize$estimate[2],
                                  digits = digits),
                     " ",
                     100 * attr(htest$conf.int, "conf.level"),
                     "% C.I.[",
                     rounder_stat(tosty$effsize$lower.ci[2],
                                  digits = digits),
                     ", ",
                     rounder_stat(tosty$effsize$upper.ci[2],
                                  digits = digits),
                     "]",
                     ").")

  text2 = paste(method_state, stat_text,
                claim_text,
                sep = " ")

  return(text2)
}
