#' @title One, two, and paired samples hypothesis tests
#'
#' @description
#' `r lifecycle::badge("maturing")`
#'
#' Performs one or two sample t-tests or Wilcoxon-Mann-Whitney rank-based tests with expanded options compared to [t.test], [brunner_munzel], or [wilcox.test].
#' @param test a character string specifying what type of hypothesis test to use. Options are limited to "wilcox.test", "t.test", or "brunner_munzel". You can specify just the initial letter.
#' @inheritParams t_TOST
#' @inheritParams z_cor_test
#' @inheritParams brunner_munzel
#' @param 	mu a number specifying an optional parameter used to form the null hypothesis. See ‘Details’.
#' @param ...  further arguments to be passed to or from methods.
#' @details
#' The type of test, t-test/Wilcoxon-Mann-Whitney/Brunner-Munzel, can be selected with the `"test"` argument.
#'
#'
#' \code{alternative = "greater"} is the alternative that x is larger than y (on average).
#' If \code{alternative = "equivalence"} then the alternative is that the difference between x and y is between the two null values `mu`..
#' If \code{alternative = "minimal.effect"} then the alternative is that the difference between x and y is less than the lowest null value or greater than the highest.
#'
#' For more details on each possible test ([brunner_munzel], [stats::t.test], or [stats::wilcox.test]), please read their individual documentation.
#'
#' @return A list with class `"htest"` containing the following components:
#'
#'   - statistic: the value of the t-statistic.
#'   - parameter: the degrees of freedom for the t-statistic.
#'   - p.value: the p-value for the test.
#'   - conf.int: a confidence interval for the mean appropriate to the specified alternative hypothesis.
#'   - estimate: the estimated mean or difference in means depending on whether it was a one-sample test or a two-sample test.
#'   - null.value: the specified hypothesized value of the mean or mean difference. May be 2 values.
#'   - stderr: the standard error of the mean (difference), used as denominator in the t-statistic formula.
#'   - alternative: a character string describing the alternative hypothesis.
#'   - method: a character string indicating what type of t-test was performed.
#'   - data.name: a character string giving the name(s) of the data.
#'
#' @examples
#' data(mtcars)
#' simple_htest(mpg ~ am,
#' data = mtcars,
#' alternative = "e",
#' mu = 3)
#' @family TOST
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

