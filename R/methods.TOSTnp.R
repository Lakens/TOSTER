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
                        digits = getOption("digits"),
                        ...){
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
  print(x$TOST)
  cat("\n")
  cat("Effect Sizes \n")
  print(x$effsize)
  cat("\n")
  if("boot" %in% names(x)){
    cat("Note: percentile boostrap method utilized.")
  }
  cat("\n")

}
