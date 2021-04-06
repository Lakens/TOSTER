#' Methods for TOSTt objects
#'
#' Methods defined for objects returned from the agree functions.
#'
#' @param x object of class \code{TOSTt} as returned from the reli_stats function
#' @param ... further arguments passed through, see description of return value
#'   for details.
#'   \code{\link{reli_stats}}.
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

print.TOSTt <- function(x,...){
  cat("\n")
  cat("Coefficient of Variation (%): ",round(x$cv*100,2))
  cat("\n")
  cat("Standard Error of Measurement (SEM): ",round(x$SEM,4))
  cat("\n")
  cat("Standard Error of the Estimate (SEE): ",round(x$SEE,4))
  cat("\n")
  cat("Standard Error of Prediction (SEP): ",round(x$SEP,4))
  cat("\n")
  cat("\n")
  cat("Intraclass Correlation Coefficients")
  cat("\n")
  print(x$icc,digits=4)
  cat("\n")
}

#' @rdname TOSTt-methods
#' @method plot TOSTt
#' @import ggplot2
#' @export

plot.TOSTt <- function(x,  ...){

  return(x$plot.reliability)

}
