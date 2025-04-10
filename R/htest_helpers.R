#' @title Helper Functions for Working with 'htest' Objects
#'
#' @description
#' A collection of utility functions designed to help interpret, display, and standardize
#' information from objects of class 'htest' (hypothesis test results). These functions
#' make it easier to extract, format, and report statistical results from various test
#' functions in R.
#'
#' @param htest An S3 object of class 'htest', such as (but not limited to) output from `t.test()`, `cor.test()`,
#'   `wilcox.test()`, or TOSTER functions converted with `as_htest()`.
#' @param test_statistics A logical variable indicating whether to display the test statistics
#'   in the output (default = TRUE).
#' @param show_ci A logical variable indicating whether to display the confidence interval
#'   in the output (default = TRUE).
#' @param extract_names A logical variable indicating whether to take the names from the S3 object
#'   (i.e., statistic for `t.test()` would be "t") (default = TRUE).
#' @param alpha The significance level to use for determining statistical significance.
#'   If NULL (default), it will be extracted from the confidence interval of the htest object
#'   or default to 0.05.
#' @param digits Integer indicating the number of decimal places to display in the output
#'   (default = 3).
#'
#' @details
#' The package provides two main helper functions:
#'
#' 1. `df_htest()`: Converts an 'htest' object to a data frame with standardized columns,
#'    making it easier to combine multiple test results or export them for further analysis.
#'
#' 2. `describe_htest()`: Generates a formatted text description of the test results,
#'    following APA style guidelines and providing a complete statistical report with
#'    test statistics, p-values, effect sizes, and confidence intervals.
#'
#' These functions work with standard R hypothesis tests (e.g., `t.test()`, `wilcox.test()`,
#' `cor.test()`) as well as TOSTER-specific tests that have been converted to 'htest' format
#' using the `as_htest()` function.
#'
#'
#' @return
#' * `df_htest()`: Returns a data frame containing the formatted test information.
#' * `describe_htest()`: Returns a character string with a formatted description of the test results.
#'
#' @examples
#' # Example 1: Working with a standard t-test
#' t_result <- t.test(extra ~ group, data = sleep)
#'
#' # Convert to data frame
#' df_htest(t_result)
#'
#' # Generate formatted description
#' describe_htest(t_result)
#'
#' # Example 2: Working with a TOST result
#' tost_result <- t_TOST(extra ~ group, data = sleep, eqb = 1)
#' htest_conv <- as_htest(tost_result)
#' describe_htest(htest_conv)
#'
#' # Example 3: Customizing output format
#' df_htest(t_result, test_statistics = TRUE, show_ci = FALSE)
#' describe_htest(t_result, alpha = 0.01, digits = 2)
#'
#' # Example 4: Working with correlation tests
#' cor_result <- cor.test(mtcars$mpg, mtcars$wt)
#' df_htest(cor_result)
#' describe_htest(cor_result)
#'
#' @name htest-helpers
#' @family htest
#' @export
df_htest = function(htest,
                    test_statistics = TRUE,
                    show_ci = TRUE,
                    extract_names = TRUE){

  if(!(class(htest) %in% c("htest"))){
    stop("htest must be of the class htest")
  }

  df1 = data.frame(method = htest$method)

  if(test_statistics) {

  # Get columns
  if(!is.null(htest$statistic)){
    if(length(htest$statistic) == 1){
      df1[ifelse(extract_names, names(htest$statistic),"statistic")] <- unname(htest$statistic)
    } else{
      for(i in 1:length(htest$statistic)){
        parm_name = ifelse(extract_names, names(htest$statistic)[i],paste0("statistic",i))
        df1[parm_name] <- unname(htest$statistic)[i]
      }
    }
  }

  if(!is.null(htest$parameter)){
    if(length(htest$parameter) == 1){
      df1[ifelse(extract_names, names(htest$parameter),"parameter")] <- unname(htest$parameter)
    } else{
      for(i in 1:length(htest$parameter)){
        parm_name = ifelse(extract_names, names(htest$parameter)[i],paste0("parameter",i))
        df1[parm_name] <- unname(htest$parameter)[i]
      }
    }
  }

  if(!is.null(htest$p.value)){
    df1$p.value <- unname(htest$p.value)
  }

  }

  if(!is.null(htest$estimate)){
    if(grepl("two sample t-test",htest$method,ignore.case = TRUE) && length(htest$estimate) > 1){
      htest$estimate = htest$estimate[1] - htest$estimate[2]
      names(htest$estimate) = c("mean difference")
    }
    if(length(htest$estimate) == 1){
      df1[ifelse(extract_names, names(htest$estimate),"estimate")] <- unname(htest$estimate)
    } else{
      for(i in 1:length(htest$estimate)){
        parm_name = ifelse(extract_names, names(htest$estimate)[i],paste0("estimate",i))
        df1[parm_name] <- unname(htest$estimate)[i]
      }
    }
  }

  if(!is.null(htest$stderr)){
    if(length(htest$stderr) == 1){
      df1["SE"] <- unname(htest$stderr)
    } else{
      for(i in 1:length(htest$stderr)){
        parm_name = paste0("SE",i)
        df1[parm_name] <- unname(htest$stderr)[i]
      }
    }
  }
  if(show_ci == TRUE){

  if(!is.null(htest$conf.int)){

    df1$lower.ci = min(htest$conf.int)
    df1$upper.ci = max(htest$conf.int)
    df1$conf.level = attr(htest$conf.int, "conf.level")

  }

}
  if(!is.null(htest$alternative)){
    df1$alternative = htest$alternative
  }

  if(!is.null(htest$null.value)){

    if(length(htest$null.value) == 1){
      df1["null"] <- unname(htest$null.value)
    } else{
      for(i in 1:length(htest$null.value)){
        parm_name = paste0("null",i)
        df1[parm_name] <- unname(htest$null.value)[i]
      }
    }
  }
  # Rename
  # names(df1)[names(df1) == 'old.var.name'] <- 'new.var.name'

  df_final = df1
  return(df_final)

}

#' @rdname htest-helpers
#' @export

describe_htest = function(htest,alpha = NULL,digits = 3){
  if(!(class(htest) %in% c("htest"))){
    stop("htest must be of the class htest")
  }

  if(is.null(htest$p.value)){
    stop("htest must have p.value")
  }
  # Set alpha
  if(is.null(alpha)){
    if(is.null(htest$conf.int)){
      alpha = 0.05
      message("No alpha level set or confidence interval provided. Defaulting to 0.05")
    } else{
      alpha = 1-attr(htest$conf.int,"conf.level")
      if(!is.null(htest$alternative)){
        if(htest$alternative %in% c("equivalence","minimal.effect")){
          alpha = alpha/2
        }
      }
    }
  }
  sig_state = ifelse(
    htest$p.value < alpha,
    "statistically significant",
    "not statistically significant"
  )

  hyp_state = ifelse(
    htest$p.value < alpha,
    "The null hypothesis can be rejected.",
    "The null hypothesis cannot be rejected."
  )

  if(!is.null(htest$parameter)){

    for (i in 1:length(htest$parameter)) {
      if (i == 1) {
        par_state = paste0(rounder_stat(unname(htest$parameter[i])))
      } else{
        par_state = paste0(par_state,", ",rounder_stat(unname(htest$parameter[i])))
      }

    }
    par_state = paste0("(",par_state,")")

  } else {
    par_state = ""
  }

  if(!is.null(htest$p.value)){
    pval_state = printable_pval(htest$p.value,
                                digits = digits)
  }else{
    pval_state = ""
  }

  if(!is.null(htest$statistic)){

    stat_name = names(htest$statistic)
    stat_est = unname(htest$statistic)

    stat_state = paste0(stat_name,
                        par_state,
                        " = ",
                        rounder_stat(stat_est),
                        ", ",
                        pval_state)

  } else{
    stat_state = pval_state
  }

  if(!is.null(htest$estimate) && !is.null(htest$conf.int)){
    for(i in 1:length(htest$estimate)){
      if(i == 1){
        est_state = paste0(names(htest$estimate[i]),
                           " = ",
                           rounder_stat(unname(htest$estimate[i]),
                                        digits = digits),
                           ", ")
      } else {
        est_state = paste0(est_state,
                           names(htest$estimate[i]),
                           " = ",
                           rounder_stat(unname(htest$estimate[i]),
                                        digits = digits),
                           ", ")
      }

    }

    conf_state = paste0(
      100 * attr(htest$conf.int, "conf.level"),
      "% C.I.[",
      rounder_stat(min(htest$conf.int),
                   digits = digits),
      ", ",
      rounder_stat(max(htest$conf.int),
                   digits = digits),
      "]"
    )
    est_state = paste0(est_state, conf_state)
    stat_state = paste0(stat_state,", ", est_state)

  }

  if(!is.null(htest$alternative)){

    if(!is.null(htest$null.value)) {
      if(length(htest$null.value) == 1) {
        alt.char <-
          switch(htest$alternative,
                 two.sided = "not equal to",
                 less = "less than",
                 greater = "greater than")
        alt_state = paste0("true ",
                           names(htest$null.value),
                           " is ",
                           alt.char,
                           " ",
                           htest$null.value)
        decision = ifelse(htest$p.value < alpha,
                          " At the desired error rate, it can be stated that the ",
                          " At the desired error rate, it cannot be stated that the ")

        print_state = paste0("The ",
                             htest$method,
                             " is ",
                             sig_state,
                             " (",
                             stat_state, ") at a ", alpha, " alpha-level. ",
                             hyp_state,
                             #" Based on the alternative hypothesis,",
                             decision,
                             alt_state, ".")
      } else {

        if(htest$alternative %in% c("equivalence","minimal.effect")){
          if(htest$alternative == "equivalence"){

            alt_state = paste0("true ",
                               names(htest$null.value)[1],
                               " is ",
                               "between",
                               " ",
                               htest$null.value[1], " and ",
                               htest$null.value[2])

          }
          if(htest$alternative == "minimal.effect"){

            alt_state = paste0("true ",
                               names(htest$null.value)[1],
                               " is ",
                               "less than ",
                               htest$null.value[1], " or ",
                               "greater than ",
                               htest$null.value[2])

          }
          decision = ifelse(htest$p.value < alpha,
                            " At the desired error rate, it can be stated that the ",
                            " At the desired error rate, it cannot be stated that the ")

          print_state = paste0("The ",
                               htest$method,
                               " is ",
                               sig_state,
                               " (",
                               stat_state, ") at a ", alpha, " alpha-level. ",
                               hyp_state,
                               #" Based on the alternative hypothesis,",
                               decision,
                               alt_state, ".")
        } else{
          print_state = paste0("The ",
                               htest$method,
                               " is ",
                               sig_state,
                               " effect, ",
                               stat_state, ", at a ", alpha, " alpha-level. ",
                               hyp_state,
                               " alternative: ", htest$alternative, " with ",
                               htest$null.value,
                               ".")
        }


      }
    } else {
      print_state = paste0("The ",
                           htest$method,
                           " is ",
                           sig_state,
                           ", ",
                           stat_state, ", at a ", alpha, " alpha-level. ",
                           hyp_state,
                           " alternative: ", htest$alternative,".")

    }


  } else{
    print_state = paste0("The ",
                         htest$method,
                         " is ",
                         sig_state,
                         ", ",
                         stat_state, ", at a ", alpha, " alpha-level. ",
                         hyp_state)
  }

  return(print_state)
}

rounder_stat = function(number,
                        digits = 3){
  if(is.na(number) || !is.numeric(number)){
    number = NA
  } else{
    cutoff = 1*10^(-1*digits)
    if(number < cutoff){
      number = signif(number, digits = digits)
    } else{
      number = round(number, digits = digits)
    }
  }

  return(number)
}

printable_pval = function(pval,
                          digits = 3){
  cutoff = 1*10^(-1*digits)
  if(pval < cutoff){
    pval = paste0("p < ",cutoff)
  } else{
    pval = paste0("p = ",round(pval, digits = digits))
  }

  return(pval)
}
