#' @title Convert TOSTER Results to Class 'htest'
#'
#' @description
#' `r lifecycle::badge('maturing')`
#'
#' Converts a TOSTER result object of class 'TOSTt' or 'TOSTnp' to a list of class 'htest',
#' making it compatible with standard R hypothesis testing functions and workflows.
#'
#' @param TOST A TOSTER result object of class 'TOSTt' or 'TOSTnp'.
#'
#' @details
#' This function allows you to convert the specialized TOSTER result objects to the standard
#' 'htest' class used by most R hypothesis testing functions (e.g., `t.test()`, `cor.test()`).
#' This enables:
#'
#' * Integration with other statistical functions that expect 'htest' objects
#' * Using helper functions like `df_htest()` or `describe_htest()`
#' * Consistent reporting and interpretation of results
#'
#'
#' @return
#' Returns a list of class 'htest' containing the following components:
#'
#' * **statistic**: The value of the test statistic (t for TOSTt, WMW for TOSTnp).
#' * **parameter**: The degrees of freedom of the test statistic (df for TOSTt, NULL for TOSTnp).
#' * **p.value**: The p-value of the test.
#' * **estimate**: Estimated difference in raw units.
#' * **null.value**: Equivalence bounds.
#' * **alternative**: A character string describing the alternative hypothesis ("equivalence" or "minimal.effect").
#' * **method**: A character string indicating the performed test.
#' * **data.name**: A character string giving the names of the data.
#' * **conf.int**: The confidence interval of the difference.
#'
#' @examples
#' # Example 1: Converting TOST t-test results to htest
#' res1 <- t_TOST(formula = extra ~ group, data = sleep, eqb = .5, smd_ci = "goulet")
#' htest_result <- as_htest(res1)
#' htest_result  # Print the htest object
#'
#' # Example 2: Using the converted result with htest helpers
#' describe_htest(htest_result)
#' df_htest(htest_result)
#'
#' # Example 3: Converting a non-parametric TOST result
#' res2 <- wilcox_TOST(extra ~ group, data = sleep, eqb = 2)
#' as_htest(res2)
#'
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

