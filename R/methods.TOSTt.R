#' Methods for TOSTt objects
#'
#' Methods defined for objects returned from the agree functions.
#'
#' @param x object of class \code{TOSTt} as returned from the reli_stats function
#' @param ... further arguments passed through, see description of return value
#'   for details.
#'   \code{\link{TOSTt}}.
#' @return
#' \describe{
#'   \item{\code{print}}{Prints short summary of the Limits of Agreement}
#'   \item{\code{plot}}{Returns a plot of the data points used in the reliability analysis}
#' }
#'
#' @name TOSTt-methods


### methods for TOSTt objects

#' @rdname TOSTt-methods
#' @method print TOSTt
#' @export

print.TOSTt <- function(x,digits = getOption("digits"), prefix = "\t",...){
  cat("\n")
  cat(strwrap(x$method, prefix = prefix), sep = "\n")
  cat("\n")
  cat(x$hypothesis$test_hypothesis, "\n", sep = "")
  cat("\n")
  cat("TOST Results \n")
  print(x$TOST)
  cat("\n")
  cat("Effect Sizes \n")
  print(x$effsize)
  cat("\n")

}

#' @rdname TOSTt-methods
#' @method plot TOSTt
#' @import ggplot2
#' @export

plot.TOSTt <- function(x,  ...){

  return(x$plot.reliability)

}
