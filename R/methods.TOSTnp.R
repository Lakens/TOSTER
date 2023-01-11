#' Methods for TOSTnp objects
#'
#' Methods defined for objects returned from the agree functions.
#'
#' @param x object of class \code{TOSTnp} as returned from the reli_stats function
#' @param digits Number of digits to print for p-values
#' @param ... further arguments passed through, see description of return value
#'   for details.
#'   \code{\link{TOSTnp-methods}}.
#' @return
#' \describe{
#'   \item{\code{print}}{Prints short summary of the Limits of Agreement}
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
