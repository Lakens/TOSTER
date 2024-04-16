#' @title Convert to class 'htest'
#'
#' @description
#' `r lifecycle::badge('maturing')`
#'
#' Convert a TOSTER result object of class 'TOSTt' or 'TOSTnp' to a list of class 'htest'.
#'
#' @param TOST A TOSTER result object of class 'TOSTt' or 'TOSTnp'.
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
#' res1 = t_TOST(formula = extra ~ group,data = sleep,eqb = .5,smd_ci = "goulet")
#' as_htest(res1)
#' @family htest
#' @export

as_htest = function(TOST) {
  if(!(class(TOST) %in% c("TOSTt", "TOSTnp"))){
    stop("Class cannot be converted to htest with this function.")
  }

  alt_text = ifelse(grepl("equ", TOST$hypothesis,
                          ignore.case = TRUE),
                    "equivalence",
                    "minimal.effect")
  # get row

  TOSTp = TOST$TOST[2:3, ]
  TOSTp = switch(alt_text,
                 "equivalence" = TOSTp[which.max(TOSTp$p.value), ],
                 "minimal.effect" = TOSTp[which.min(TOSTp$p.value), ]
  )

  # assign
  statistic <- switch(class(TOST),
                      TOSTt = TOSTp$t,
                      TOSTnp = TOSTp$statistic)
  names(statistic) <- switch(class(TOST),
                             TOSTt = "t",
                             TOSTnp = "WMW")

  parameter <- switch(class(TOST),
                      TOSTt = TOSTp$df,
                      TOSTnp = NULL)
  names(parameter) <- switch(class(TOST),
                             TOSTt = "df",
                             TOSTnp = NULL)

  which_row = rownames(TOSTp)
  which_eqb <- switch(
    class(TOST),
    TOSTt = ifelse(grepl("Lower", which_row), 2, 3),
    TOSTnp = ifelse(grepl("Lower", which_row), 1, 2)
  )
  eqb = TOST$eqb

  null.value <- switch(class(TOST),
                       TOSTt = as.numeric(eqb[1,2:3]),
                       TOSTnp = eqb)
  if(grepl("one",TOST$method, ignore.case = TRUE)){
    names(null.value) <- switch(class(TOST),
                                TOSTt = rep("mean",length(null.value)),
                                TOSTnp = rep("location",length(null.value)))
  } else {
    names(null.value) <- switch(class(TOST),
                                TOSTt = rep("mean difference",length(null.value)),
                                TOSTnp = rep("location shift",length(null.value)))
  }


  stderr <- switch(class(TOST),
                   TOSTt = TOSTp$SE,
                   TOSTnp = NULL)

  p.value <- switch(class(TOST),
                    TOSTt = TOSTp$p.value,
                    TOSTnp = TOSTp$p.value)

  method <- switch(class(TOST),
                   TOSTt = TOST$method,
                   TOSTnp = TOST$method)

  alternative <- switch(class(TOST),
                        TOSTt = "one.sided",
                        TOSTnp = "one.sided")

  alt_bound <-
    ifelse(grepl("Lower", which_row), "lower", "upper")




    conf.int <- c(TOST$effsize$lower.ci[1],
                  TOST$effsize$upper.ci[1])

    attr(conf.int, "conf.level") <- TOST$effsize$conf.level[1]


    estimate <- TOST$effsize$estimate[1]

    names(estimate) = names(null.value)[1]


    htest <- list(
      statistic = statistic,
      parameter = parameter,
      p.value = p.value,
      estimate = estimate,
      null.value = null.value,
      alternative = alt_text,
      method = method,
      data.name = TOST$data.name,
      conf.int = conf.int
    )


  class(htest) <- "htest"
  htest
}

