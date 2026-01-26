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

#' @title Plot Estimate from 'htest' Object
#'
#' @description
#' `r lifecycle::badge('stable')`
#'
#' Creates a simple point estimate plot with confidence interval from any 'htest' object
#' that contains an estimate and confidence interval. This provides a visual representation
#' of the effect size and its uncertainty, similar to a forest plot.
#'
#' @param htest An S3 object of class 'htest' containing at minimum an `estimate` and
#'   `conf.int` component. Examples include output from `t.test()`, `cor.test()`,
#'   or TOSTER functions converted with `as_htest()`.
#' @param alpha Significance level for determining the confidence level label.
#' @param describe Logical. If TRUE (default), includes a concise statistical description
#'   in the plot subtitle showing the test statistic, p-value, estimate, confidence interval,
#'   and the null hypothesis.
#'
#' @details
#' The function creates a horizontal point-range plot showing:
#' \itemize{
#'   \item Point estimate (black dot)
#'   \item Confidence interval (horizontal line)
#'   \item Null value(s) as dashed vertical reference line(s)
#' }
#'
#' For two-sample t-tests, R's `t.test()` returns both group means as the estimate
#' rather than their difference. This function automatically computes the difference to display
#' a single meaningful estimate with its confidence interval.
#'
#' If the 'htest' object contains equivalence bounds (two values in `null.value`),
#' both bounds are displayed as dashed vertical lines.
#'
#' When `describe = TRUE`, the plot includes a three-line subtitle:
#' \enumerate{
#'   \item Test statistic and p-value
#'   \item Point estimate and confidence interval
#'   \item Null hypothesis statement
#' }
#' The method name appears as the plot title.
#'
#' @return A `ggplot` object that can be further customized using ggplot2 functions.
#'
#' @examples
#' # Standard t-test
#' t_result <- t.test(extra ~ group, data = sleep)
#' plot_htest_est(t_result)
#'
#' # One-sample t-test
#' t_one <- t.test(sleep$extra, mu = 0)
#' plot_htest_est(t_one)
#'
#' # Correlation test
#' cor_result <- cor.test(mtcars$mpg, mtcars$wt)
#' plot_htest_est(cor_result)
#'
#' # TOST result converted to htest
#' tost_res <- t_TOST(extra ~ group, data = sleep, eqb = 1)
#' plot_htest_est(as_htest(tost_res))
#'
#' # Without description
#' plot_htest_est(t_result, describe = FALSE)
#'
#' @import ggplot2
#' @import ggdist
#' @family htest
#' @export
plot_htest_est <- function(htest, alpha = NULL, describe = TRUE) {


  if (!inherits(htest, "htest")) {
    stop("Input must be an object of class
'htest'")
  }

  if (is.null(htest$estimate)) {
    stop("Cannot create estimate plot: htest object has no estimate")
  }

  if (is.null(htest$conf.int)) {
    stop("Cannot create estimate plot: htest object has no confidence interval")
  }

  # Handle two-sample t-test case where estimate contains both group means
  estimate <- htest$estimate
  estimate_name <- names(estimate)

  if (grepl("two sample t-test", htest$method, ignore.case = TRUE) &&
      length(estimate) > 1) {
    estimate <- estimate[1] - estimate[2]
    estimate_name <- "mean difference"
  } else if (length(estimate) > 1) {
    # For other cases with multiple estimates, warn and use first
    warning("htest object has multiple estimates; using first estimate only")
    estimate <- estimate[1]
    estimate_name <- names(htest$estimate)[1]
  }

  # Get confidence interval
  ci_lower <- min(htest$conf.int)
  ci_upper <- max(htest$conf.int)
  conf_level <- attr(htest$conf.int, "conf.level")

  if (is.null(conf_level)) {
    conf_level <- 0.95
    if (is.null(alpha)) {
      message("No confidence level found in htest object. Defaulting to 95%.")
    }
  }

  # Determine label for facet
  if (is.null(estimate_name) || length(estimate_name) == 0) {
    facet_label <- "Estimate"
  } else {
    # Capitalize first letter
    facet_label <- paste0(toupper(substr(estimate_name, 1, 1)),
                          substr(estimate_name, 2, nchar(estimate_name)))
  }

  # Create data frame for plotting (include facet_label in the data)
  df_plot <- data.frame(
    estimate = unname(estimate),
    lower.ci = ci_lower,
    upper.ci = ci_upper,
    facet_label = facet_label,
    stringsAsFactors = FALSE
  )

  # Build description for subtitle if requested
  if (describe) {
    # Build concise description similar to describe_htest but shorter
    desc_parts <- c()

    # Add test statistic if available
    if (!is.null(htest$statistic)) {
      stat_name <- names(htest$statistic)
      stat_val <- rounder_stat(unname(htest$statistic), digits = 3)

      if (!is.null(htest$parameter)) {
        par_val <- rounder_stat(unname(htest$parameter), digits = 2)
        stat_str <- paste0(stat_name, "(", par_val, ") = ", stat_val)
      } else {
        stat_str <- paste0(stat_name, " = ", stat_val)
      }
      desc_parts <- c(desc_parts, stat_str)
    }

    # Add p-value if available
    if (!is.null(htest$p.value)) {
      desc_parts <- c(desc_parts, printable_pval(htest$p.value, digits = 3))
    }

    # Build first line: test statistic and p-value
    line1 <- paste(desc_parts, collapse = ", ")

    # Build second line: estimate and CI
    est_str <- paste0(estimate_name, " = ",
                      rounder_stat(unname(estimate), digits = 3))
    ci_str <- paste0(round(conf_level * 100), "% CI [",
                     rounder_stat(ci_lower, digits = 3), ", ",
                     rounder_stat(ci_upper, digits = 3), "]")
    line2 <- paste(est_str, ci_str, sep = ", ")

    # Build third line: null hypothesis
    line3 <- NULL
    if (!is.null(htest$null.value) && !is.null(htest$alternative)) {
      null_name <- names(htest$null.value)
      if (is.null(null_name) || length(null_name) == 0) {
        null_name <- estimate_name
      }

      if (length(htest$null.value) == 1) {
        # Standard hypothesis test - show null based on alternative
        null_rel <- switch(htest$alternative,
                           two.sided = "is equal to",
                           less = "is greater than or equal to",
                           greater = "is less than or equal to",
                           "is equal to")
        line3 <- paste0("null: ", null_name, " ", null_rel, " ",
                        rounder_stat(unname(htest$null.value), digits = 3))
      } else if (length(htest$null.value) == 2) {
        # Equivalence or minimal effect test
        null_vals <- sort(unname(htest$null.value))
        if (htest$alternative == "equivalence") {
          line3 <- paste0("null: ", null_name, " < ", rounder_stat(null_vals[1], digits = 3),
                          " or > ", rounder_stat(null_vals[2], digits = 3))
        } else if (htest$alternative == "minimal.effect") {
          line3 <- paste0("null: ", rounder_stat(null_vals[1], digits = 3),
                          " < ", null_name, " < ", rounder_stat(null_vals[2], digits = 3))
        } else {
          # Fallback for other cases with two bounds
          line3 <- paste0("null: ", null_name, " in [",
                          rounder_stat(null_vals[1], digits = 3), ", ",
                          rounder_stat(null_vals[2], digits = 3), "]")
        }
      }
    } else if (!is.null(htest$null.value)) {
      # No alternative specified, just show null value
      null_name <- names(htest$null.value)
      if (is.null(null_name) || length(null_name) == 0) {
        null_name <- estimate_name
      }
      if (length(htest$null.value) == 1) {
        line3 <- paste0("null: ", null_name, " = ",
                        rounder_stat(unname(htest$null.value), digits = 3))
      }
    }

    # Combine lines
    if (!is.null(line3)) {
      subtitle_text <- paste(line1, line2, line3, sep = "\n")
    } else {
      subtitle_text <- paste(line1, line2, sep = "\n")
    }
    title_text <- htest$method
  } else {
    subtitle_text <- NULL
    title_text <- htest$method
  }

  # Build the plot
  p <- ggplot(df_plot,
              aes(x = estimate,
                  y = 1,
                  xmin = lower.ci,
                  xmax = upper.ci)) +
    geom_pointrange() +
    facet_grid(~facet_label) +
    theme_tidybayes() +
    labs(caption = paste0(conf_level * 100, "% Confidence Interval"),
         title = title_text,
         subtitle = subtitle_text) +
    theme(strip.text = element_text(face = "bold", size = 10),
          plot.title = element_text(size = 11),
          plot.subtitle = element_text(size = 9),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())

  # Add null value reference line(s)
  if (!is.null(htest$null.value)) {
    null_vals <- unname(htest$null.value)

    if (length(null_vals) == 1) {
      # Single null value (standard hypothesis test)
      p <- p + geom_vline(xintercept = null_vals, linetype = "dashed")
    } else if (length(null_vals) == 2) {
      # Two null values (equivalence bounds)
      p <- p +
        geom_vline(xintercept = null_vals[1], linetype = "dashed") +
        geom_vline(xintercept = null_vals[2], linetype = "dashed") +
        scale_x_continuous(sec.axis = dup_axis(
          breaks = round(null_vals, 3),
          name = ""
        ))
    }
  }

  return(p)
}
