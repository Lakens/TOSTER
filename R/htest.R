#' Convert to class 'htest'
#'
#' Convert a TOSTER result object of class 'TOSTt' or 'TOSTnp' to a list of class 'htest'.
#'
#' @param result.object A TOSTER result object of class 'TOSTt' or 'TOSTnp' to a list of class 'htest.
#'
#' @return Returns a list containing a list of class 'htest' for the result of each test with the following elements:
#' \item{data.name}{A character string giving the names of the data.}
#' \item{estimate}{Estimated difference in raw units.}
#' \item{method}{A character string indicating the performed test.}
#' \item{null.value}{Equivalence bound.}
#' \item{alternative}{A character string describing the alternative hypothesis.}
#' \item{parameter}{The degrees of freedom of the distribution of the test statistic.}
#' \item{statistic}{The value of the test statistic.}
#' \item{p.value}{The p-value of the test.}
#' \item{conf.int}{The confidence interval of the difference.}
#'
#' @examples
#' # To be added
#'
#'
#' # as.htest(result)
#'
#' @docType methods
#' @rdname as.htest
#' @export
setGeneric("as.htest", function(result.object) standardGeneric("as.htest"))

#' @aliases as.htest,TOSTER-method
#' @rdname as.htest
setMethod("as.htest", "TOST",
          function(result.object) {

            estimate <- switch(
              class(result.object),
              TOSTt = result.object$TOST[c("t-test"),]$t,
              TOSTnp = c(
                r.jk = result.object@r.jk,
                r.jh = result.object@r.jh,
                r.kh = result.object@r.kh
              )
            )

            conf.int <- test.object$conf.int
            if (!is.null(conf.int))
              attr(conf.int, "conf.level") <- result.object@conf.level

            statistic <- test.object$statistic
            names(statistic) <- test.object$distribution

            estimate <- switch(
              class(result.object),
              TOSTt = c(r1.jk = result.object@r1.jk, r2.hm =
                                       result.object@r2.hm),
              TOSTnp = c(
                r.jk = result.object@r.jk,
                r.jh = result.object@r.jh,
                r.kh = result.object@r.kh
              )
            )

                htest <- list(
                  statistic=statistic,
                  parameter=c(df=test.object$df),
                  p.value=test.object$p.value,
                  estimate=estimate,
                  null.value=c("difference in correlations"=result.object@null.value),
                  alternative=result.object@alternative,
                  method=all.tests[[test]],
                  data.name=data.description(result.object@data.name, result.object@var.labels),
                  conf.int=conf.int
                )
                class(htest) <- "htest"




            htest.list
          }
)
