#' @title One, Two, and Paired Samples Hypothesis Tests with Extended Options
#'
#' @description
#' `r lifecycle::badge("maturing")`
#'
#' Performs statistical hypothesis tests with extended functionality beyond standard implementations.
#' Supports t-tests, Wilcoxon-Mann-Whitney tests, and Brunner-Munzel tests with additional
#' alternatives such as equivalence and minimal effect testing.
#'
#' @section Purpose:
#' Use this function when:
#'   - You need a unified interface for different types of hypothesis tests
#'   - You want to perform equivalence testing or minimal effect testing with non-parametric methods
#'   - You need more flexibility in hypothesis testing than provided by standard functions
#'   - You want to easily switch between parametric and non-parametric methods
#'
#' @inheritParams t_TOST
#' @inheritParams z_cor_test
#' @inheritParams brunner_munzel
#' @param test a character string specifying the type of hypothesis test to use:
#'     - "t.test": Student's t-test (parametric, default)
#'     - "wilcox.test": Wilcoxon-Mann-Whitney test (non-parametric)
#'     - "brunner_munzel": Brunner-Munzel test (non-parametric)
#'
#'   You can specify just the initial letter (e.g., "t" for "t.test").
#' @param alternative the alternative hypothesis:
#'     - "two.sided": different from mu (default)
#'     - "less": less than mu
#'     - "greater": greater than mu
#'     - "equivalence": between specified bounds
#'     - "minimal.effect": outside specified bounds
#' @param mu a number or vector specifying the null hypothesis value(s):
#'     - For standard alternatives (two.sided, less, greater): a single value (default: 0 for t-test/wilcox.test, 0.5 for brunner_munzel)
#'     - For equivalence/minimal.effect: either a single value (symmetric bounds will be created) or a vector of two values representing the lower and upper bounds
#' @param ... further arguments to be passed to or from the underlying test functions.
#'
#' @details
#' This function provides a unified interface to several common hypothesis tests with expanded
#' alternative hypotheses, particularly for equivalence testing and minimal effect testing.
#'
#' When `alternative = "equivalence"`, the test evaluates whether the effect is contained
#' within the bounds specified by `mu`. This corresponds to the alternative hypothesis that
#' the true effect is between the specified bounds. The function performs two one-sided tests and
#' returns the most conservative result (highest p-value).
#'
#' When `alternative = "minimal.effect"`, the test evaluates whether the effect is outside
#' the bounds specified by `mu`. This corresponds to the alternative hypothesis that the true
#' effect is either less than the lower bound or greater than the upper bound. The function performs
#' two one-sided tests and returns the most significant result (lowest p-value).
#'
#' For standard alternatives ("two.sided", "less", "greater"), the function behaves similarly to
#' the underlying test functions with some additional standardization in the output format.
#'
#' The interpretation of `mu` depends on the test used:
#'   - For t-test and wilcox.test: mu represents the difference in means/medians (default: 0)
#'   - For brunner_munzel: mu represents the probability that a randomly selected value from the first sample exceeds a randomly selected value from the second sample (default: 0.5)
#'
#'
#' If `mu` is a single value for equivalence or minimal effect alternatives, symmetric bounds
#' will be created automatically:
#'   - For t-test and wilcox.test: bounds become c(mu, -mu)
#'   - For brunner_munzel: bounds become c(mu, abs(mu-1))
#'
#'
#' @return A list with class `"htest"` containing the following components:
#'
#'   - statistic: the value of the test statistic.
#'   - parameter: the parameter(s) for the test statistic (e.g., degrees of freedom for t-tests).
#'   - p.value: the p-value for the test.
#'   - conf.int: a confidence interval appropriate to the specified alternative hypothesis.
#'   - estimate: the estimated effect (e.g., mean difference for t-tests, probability estimate for Brunner-Munzel).
#'   - null.value: the specified hypothesized value(s). For equivalence and minimal effect tests, this will be two values.
#'   - stderr: the standard error of the estimate (for t-tests).
#'   - alternative: a character string describing the alternative hypothesis.
#'   - method: a character string indicating what type of test was performed.
#'   - data.name: a character string giving the name(s) of the data.
#'
#' @examples
#' # Example 1: Basic t-test with equivalence alternative
#' # Testing if the difference in mpg between automatic and manual transmission cars
#' # is equivalent within ±3 units
#' data(mtcars)
#' simple_htest(mpg ~ am, data = mtcars, alternative = "equivalence", mu = 3)
#'
#' # Example 2: Using a non-parametric test with minimal effect alternative
#' # Testing if the effect of transmission type on mpg is meaningfully large
#' # (either less than -2 or greater than 2)
#' simple_htest(mpg ~ am, data = mtcars,
#'              test = "wilcox",
#'              alternative = "minimal.effect",
#'              mu = c(-2, 2))
#'
#' # Example 3: Paired samples test
#' # Using the sleep dataset to test if drug has an effect on sleep
#' data(sleep)
#' with(sleep, simple_htest(x = extra[group == 1],
#'                         y = extra[group == 2],
#'                         paired = TRUE,
#'                         alternative = "greater"))
#'
#' # Example 4: Brunner-Munzel test
#' # Testing if values in one group tend to exceed values in another
#' set.seed(123)
#' group1 <- rnorm(20, mean = 5, sd = 1)
#' group2 <- rnorm(20, mean = 6, sd = 2)
#' simple_htest(x = group1, y = group2,
#'              test = "brunner_munzel",
#'              alternative = "less")
#'
#' @family TOST
#' @family htest
#' @name simple_htest
#' @export simple_htest


#simple_htest <- setClass("simple_htest")
simple_htest <- function(x, ...,
                         paired = FALSE,
                         alternative = c("two.sided",
                                         "less",
                                         "greater",
                                         "equivalence",
                                         "minimal.effect"),
                         mu = NULL,
                         alpha = 0.05){
  UseMethod("simple_htest")
}
#' @family htest
#' @rdname simple_htest
#' @method simple_htest default
#' @export

# @method simple_htest default
simple_htest.default = function(x,
                                y = NULL,
                                test = c("t.test","wilcox.test","brunner_munzel"),
                                paired = FALSE,
                                alternative = c("two.sided",
                                                "less",
                                                "greater",
                                                "equivalence",
                                                "minimal.effect"),
                                mu = NULL,
                                alpha = 0.05,
                                ...) {
 alternative = match.arg(alternative)
 test = match.arg(test)
 if(is.null(mu)){
   if(test == "brunner_munzel"){
     mu = 0.5
     message(paste0("mu set to ", mu))
   } else{
     mu = 0
     message(paste0("mu set to ", mu))
   }
 }

 if(alternative %in% c("equivalence","minimal.effect")){

   if(length(mu) == 1){
     if(mu ==  0 && test %in% c("t.test","wilcox.test")){
       stop("mu cannot be zero if alternative is equivalence or minimal.effect")
     }

     if(mu ==  0.5 && test %in% c("brunner_munzel")){
       stop("mu cannot be zero if alternative is equivalence or minimal.effect")
     }

     if(test %in% c("t.test","wilcox.test")){
       mu = c(mu, -1*mu)
     } else {
       mu = c(mu, abs(mu-1))
     }
   }


   lo_bound = min(mu)
   hi_bound = max(mu)

   ci_test = switch(
     test,
     t.test = t.test(
       x = x,
       y = y,
       paired = paired,
       mu = 0,
       conf.level = 1 - alpha * 2,
       alternative = "two.sided",
       ...
     ),
     wilcox.test = wilcox.test(
       x = x,
       y = y,
       paired = paired,
       mu = 0,
       conf.int = TRUE,
       conf.level = 1 - alpha * 2,
       alternative = "two.sided",
       ...
     ),
     brunner_munzel = brunner_munzel(
       x = x,
       y = y,
       paired = paired,
       mu = 0.5,
       #conf.int = TRUE,
       alpha = alpha*2,
       alternative = "two.sided",
       ...
     )
   )

   if(alternative == "equivalence"){

     lo_test = switch(
       test,
       t.test = t.test(
       x = x,
       y = y,
       paired = paired,
       mu = lo_bound,
       alternative = "greater",
       ...
     ),
     wilcox.test = wilcox.test(
       x = x,
       y = y,
       paired = paired,
       mu = lo_bound,
       alternative = "greater",
       ...),
     brunner_munzel = brunner_munzel(
       x = x,
       y = y,
       paired = paired,
       mu = lo_bound,
       alternative = "greater",
       ...))
     hi_test = switch(
       test,
       t.test = t.test(
         x = x,
         y = y,
         paired = paired,
         mu = hi_bound,
         alternative = "less",
         ...
       ),
       wilcox.test = wilcox.test(
         x = x,
         y = y,
         paired = paired,
         mu = hi_bound,
         alternative = "less",
         ...),
       brunner_munzel = brunner_munzel(
         x = x,
         y = y,
         paired = paired,
         mu = hi_bound,
         alternative = "less",
         ...))


     if(hi_test$p.value >= lo_test$p.value){
       rval = hi_test
     } else {
       rval = lo_test
     }

     name_val = names(ci_test$null.value)
     rval$conf.int = ci_test$conf.int
     rval$alternative = alternative
     rval$null.value = c(lo_bound, hi_bound)
     names(rval$null.value) = rep(name_val,2)

   }

   if(alternative == "minimal.effect"){
     two_test = switch(
       test,
       t.test = t.test(
         x = x,
         y = y,
         mu = 0,
         paired = paired,
         alternative = "two.sided",
         ...
       ),
       wilcox.test = wilcox.test(
         x = x,
         y = y,
         paired = paired,
         mu = 0,
         alternative = "two.sided",
         ...),
       brunner_munzel = brunner_munzel(
         x = x,
         y = y,
         mu = 0.5,
         paired = paired,
         alternative = "two.sided",
         ...))

     lo_test = switch(
       test,
       t.test = t.test(
         x = x,
         y = y,
         paired = paired,
         mu = lo_bound,
         alternative = "less",
         ...
       ),
       wilcox.test = wilcox.test(
         x = x,
         y = y,
         paired = paired,
         mu = lo_bound,
         alternative = "less",
         ...),
       brunner_munzel = brunner_munzel(
         x = x,
         y = y,
         paired = paired,
         mu = lo_bound,
         alternative = "less",
         ...))

     hi_test = switch(
       test,
       t.test = t.test(
         x = x,
         y = y,
         paired = paired,
         mu = hi_bound,
         alternative = "greater",
         ...
       ),
       wilcox.test = wilcox.test(
         x = x,
         y = y,
         paired = paired,
         mu = hi_bound,
         alternative = "greater",
         ...),
       brunner_munzel = brunner_munzel(
         x = x,
         y = y,
         paired = paired,
         mu = hi_bound,
         alternative = "greater",
         ...))

     if(hi_test$p.value <= lo_test$p.value){
       rval = hi_test
     } else {
       rval = lo_test
     }

       if(two_test$p.value >= rval$p.value){
         message("MET test may have higher error rates than a nil two-tailed test. Consider wider equivalence bounds.")
       }

     name_val = names(ci_test$null.value)
     rval$conf.int = ci_test$conf.int
     rval$alternative = alternative
     rval$null.value = c(lo_bound, hi_bound)
     names(rval$null.value) = rep(name_val,2)

   }

 } else {

   rval = switch(
     test,
     t.test = t.test(
       x = x,
       y = y,
       paired = paired,
       mu = mu,
       conf.level = 1 - alpha,
       alternative = alternative,
       ...
     ),
     wilcox.test = wilcox.test(
       x = x,
       y = y,
       paired = paired,
       conf.int = TRUE,
       conf.level = 1 - alpha,
       alternative = alternative,
       mu = mu,
       ...
     ),
     brunner_munzel = brunner_munzel(
       x = x,
       y = y,
       paired = paired,
       alpha = alpha,
       alternative = alternative,
       mu = mu,
       ...
     )
   )

 }

  return(rval)

}

#' @rdname simple_htest
#' @method simple_htest formula
#' @export

simple_htest.formula = function(formula,
                                data,
                                subset,
                                na.action, ...) {

  if(missing(formula)
     || (length(formula) != 3L)
     || (length(attr(terms(formula[-2L]), "term.labels")) != 1L))
    stop("'formula' missing or incorrect")
  
  # Check for paired argument in ... and reject it
  dots <- list(...)
  if("paired" %in% names(dots))
    stop("cannot use 'paired' in formula method")
  
  m <- match.call(expand.dots = FALSE)
  if(is.matrix(eval(m$data, parent.frame())))
    m$data <- as.data.frame(data)
  ## need stats:: for non-standard evaluation
  m[[1L]] <- quote(stats::model.frame)
  m$... <- NULL
  mf <- eval(m, parent.frame())
  DNAME <- paste(names(mf), collapse = " by ")
  names(mf) <- NULL
  response <- attr(attr(mf, "terms"), "response")
  g <- factor(mf[[-response]])
  if(nlevels(g) != 2L)
    stop("grouping factor must have exactly 2 levels")
  DATA <- setNames(split(mf[[response]], g), c("x", "y"))
  y <- do.call("simple_htest", c(DATA, list(...)))
  y$data.name <- DNAME
  y

}

