#' Convert to class 'htest'
#'
#' Convert a TOSTER result object of class 'TOSTt' or 'TOSTnp' to a list of class 'htest'.
#'
#' @param result.object A TOSTER result object of class 'TOSTt' or 'TOSTnp'.
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
#' @export

as_htest = function(result.object) {
  if(!(class(result.object) %in% c("TOSTt", "TOSTnp"))){
    stop("Class cannot be converted to htest with this function.")
  }


  # get row

  TOSTp = result.object$TOST[2:3, ]
  TOSTp = TOSTp[which.max(TOSTp$p.value), ]

  # assign
  statistic <- switch(class(result.object),
                      TOSTt = TOSTp$t,
                      TOSTnp = TOSTp$statistic)
  names(statistic) <- switch(class(result.object),
                             TOSTt = "t",
                             TOSTnp = "WMW")

  parameter <- switch(class(result.object),
                      TOSTt = TOSTp$df,
                      TOSTnp = NULL)
  names(parameter) <- switch(class(result.object),
                             TOSTt = "df",
                             TOSTnp = NULL)

  which_row = rownames(TOSTp)
  which_eqb <- switch(
    class(result.object),
    TOSTt = ifelse(grepl("Lower", which_row), 2, 3),
    TOSTnp = ifelse(grepl("Lower", which_row), 1, 2)
  )
  eqb = result.object$eqb

  null.value <- switch(class(result.object),
                       TOSTt = eqb[1, which_eqb],
                       TOSTnp = eqb[which_eqb])
  names(null.value) <- switch(class(result.object),
                              TOSTt = "mean",
                              TOSTnp = "location shift")

  stderr <- switch(class(result.object),
                   TOSTt = TOSTp$SE,
                   TOSTnp = NULL)

  p.value <- switch(class(result.object),
                    TOSTt = TOSTp$p.value,
                    TOSTnp = TOSTp$p.value)

  method <- switch(class(result.object),
                   TOSTt = result.object$method,
                   TOSTnp = result.object$method)

  alternative <- switch(class(result.object),
                        TOSTt = "one.sided",
                        TOSTnp = "one.sided")

  alt_bound <-
    ifelse(grepl("Lower", which_row), "lower", "upper")

  if (alt_bound == "lower" &&
      grepl("equ|Equ",result.object$hypothesis)) {
    alt_text = "greater"
  } else if (alt_bound == "lower" &&
             !grepl("equ|Equ",result.object$hypothesis)) {
    alt_text = "less"
  } else if (alt_bound != "lower" &&
             !grepl("equ|Equ",result.object$hypothesis)) {
    alt_text = "greater"
  } else {
    alt_text = "less"
  }

  if (class(result.object) == "TOSTnp") {
    htest <- list(
      statistic = statistic,
      parameter = parameter,
      p.value = p.value,
      null.value = null.value,
      alternative = alt_text,
      method = method,
      data.name = TOST$data.name
    )
  } else {
    conf.int <- c(result.object$effsize$lower.ci[1],
                  result.object$effsize$upper.ci[1])
    if (!is.null(conf.int)) {
      attr(conf.int, "conf.level") <- result.object$effsize$conf.level[1]
    }
    if (grepl("One", result.object$method)) {
      estimate <- result.object$effsize$estimate[1]
      names(estimate) = "mean of x"
    } else {
      estimate <- result.object$effsize$estimate[1]
      names(estimate) = "mean difference"
    }

    htest <- list(
      statistic = statistic,
      parameter = parameter,
      p.value = p.value,
      estimate = estimate,
      null.value = null.value,
      alternative = alt_text,
      method = method,
      data.name = result.object$data.name,
      conf.int = conf.int
    )
  }

  class(htest) <- "htest"
  htest
}

