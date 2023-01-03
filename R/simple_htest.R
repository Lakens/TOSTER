#' @title One, two, and paired samples hypothesis tests
#' @description Performs one or two sample t-tests or Wilcoxon-Mann-Whitney rank-based tests with expanded options compared to \code{t.test} or \code{wilcox.test}.
#' @inheritParams t_TOST
#' @inheritParams z_cor_test
#' @param ...  further arguments to be passed to or from methods.
#' @details \code{alternative = "greater"} is the alternative that x has a larger mean than y.
#' For the one-sample case: that the mean is positive.
#'
#' If paired is TRUE then both x and y must be specified and they must be the same length.
#' Missing values are silently removed (in pairs if paired is TRUE).
#' If var.equal is TRUE then the pooled estimate of the variance is used.
#'  By default, if \code{var.equal} is FALSE then the variance is estimated separately for both groups and the Welch modification to the degrees of freedom is used.
#'
#' If the input data are effectively constant (compared to the larger of the two means) an error is generated.
#' @return A list with class \code{"htest"} containing the following components:
#' \describe{
#'   \item{\code{statistic}}{the value of the t-statistic.}
#'   \item{\code{parameter}}{the degrees of freedom for the t-statistic.}
#'   \item{\code{p.value}}{the p-value for the test.}
#'   \item{\code{conf.int}}{a confidence interval for the mean appropriate to the specified alternative hypothesis.}
#'   \item{\code{estimate}}{the estimated mean or difference in means depending on whether it was a one-sample test or a two-sample test.}
#'   \item{\code{null.value}}{the specified hypothesized value of the mean or mean difference. May be 2 values.}
#'   \item{\code{stderr}}{the standard error of the mean (difference), used as denominator in the t-statistic formula.}
#'   \item{\code{alternative}}{a character string describing the alternative hypothesis.}
#'   \item{\code{method}}{a character string indicating what type of t-test was performed.}
#'   \item{\code{data.name}}{a character string giving the name(s) of the data..}
#' }
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
                         mu = 0,
                         alpha = 0.05){
  UseMethod("simple_htest")
}

#' @rdname simple_htest
#' @method simple_htest default
#' @export

# @method simple_htest default
simple_htest.default = function(x,
                                y = NULL,
                                test = c("t.test","wilcox.test"),
                                paired = FALSE,
                                alternative = c("two.sided",
                                                "less",
                                                "greater",
                                                "equivalence",
                                                "minimal.effect"),
                                mu = 0,
                                alpha = 0.05,
                                ...) {
 alternative = match.arg(alternative)
 test = match.arg(test)

 if(alternative %in% c("equivalence","minimal.effect")){

   if(length(mu) == 1){
     if(mu ==  0){
       stop("mu cannot be zero if alternative is equivalence or minimal.effect")
     }
     mu = c(mu, -1*mu)
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
       mu = 0,
       conf.level = 1 - alpha,
       alternative = alternative,
       ...
     ),
     wilcox.test = wilcox.test(
       x = x,
       y = y,
       paired = paired,
       mu = 0,
       conf.int = TRUE,
       conf.level = 1 - alpha,
       alternative = alternative,
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

