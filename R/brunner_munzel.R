#' @title Brunner-Munzel Test
#' @description This is a generic function that performs a generalized asymptotic Brunner-Munzel test in a fashion similar to \code{"t.test"}
#' @param x a (non-empty) numeric vector of data values.
#' @param y an optional (non-empty) numeric vector of data values.
#' @param formula a formula of the form lhs ~ rhs where lhs is a numeric variable giving the data values and rhs either 1 for a one-sample or paired test or a factor with two levels giving the corresponding groups. If lhs is of class "Pair" and rhs is 1, a paired test is done.
#' @param data an optional matrix or data frame (or similar: see model.frame) containing the variables in the formula formula. By default the variables are taken from environment(formula).
#' @param paired a logical indicating whether you want a paired test.
#' @param p_null 	a number specifying an optional parameter used to form the null hypothesis (Default = 0.5). This can be thought of as the null in terms of probability of exchangeability, p = P (X < Y ) + 0.5 * P (X = Y); See ‘Details’.
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less". You can specify just the initial letter.
#' @param subset an optional vector specifying a subset of observations to be used.
#' @param na.action a function which indicates what should happen when the data contain NAs. Defaults to getOption("na.action").
#' @inheritParams t_TOST
#' @param ...  further arguments to be passed to or from methods.
#' @details
#'
#' The estimate of probability of exchangeability refers to the following:
#'
#'  p(X<Y) + .5*P(X=Y)
#'
#' @return A list with class \code{"htest"} containing the following components:
#' \describe{
#'   \item{\code{statistic}}{the value of the test statistic.}
#'   \item{\code{parameter}}{the degrees of freedom for the test statistic.}
#'   \item{\code{p.value}}{the p-value for the test.}
#'   \item{\code{conf.int}}{a confidence interval for the mean appropriate to the specified alternative hypothesis.}
#'   \item{\code{estimate}}{the estimated mean or difference in means depending on whether it was a one-sample test or a two-sample test.}
#'   \item{\code{null.value}}{the specified hypothesized value of the mean or mean difference. May be 2 values.}
#'   \item{\code{stderr}}{the standard error.}
#'   \item{\code{alternative}}{a character string describing the alternative hypothesis.}
#'   \item{\code{method}}{a character string indicating what type of test was performed.}
#'   \item{\code{data.name}}{a character string giving the name(s) of the data.}
#' }
#' @examples
#' data(mtcars)
#' #brunner_munzel(mpg ~ am,
#' #data = mtcars)
#' @references reference
#' Brunner, E., Munzel, U. (2000). The Nonparametric Behrens-Fisher Problem: Asymptotic Theory and a Small Sample Approximation. Biometrical Journal 42, 17 -25.
#'
#' Neubert, K., Brunner, E., (2006). A Studentized Permutation Test for the Nonparametric Behrens-Fisher Problem. Computational Statistics and Data Analysis.
#'
#' Munzel, U., Brunner, E. (2002). An Exact Paired Rank Test. Biometrical Journal 44, 584-593.
#' @family TOST
#' @name brunner_munzel
#' @export brunner_munzel

#brunner_munzel <- setClass("brunner_munzel")
brunner_munzel <- function(x,
                           alternative = c("two.sided",
                                           "less",
                                           "greater"),
                           p_null = 0.5,
                           alpha = 0.05,
                           ...) {

  UseMethod("brunner_munzel")
}

#' @rdname brunner_munzel
#' @importFrom stats sd cor na.omit setNames t.test terms nlm optim optimize
#' @method brunner_munzel default
#' @export

# @method brunner_munzel default
brunner_munzel.default = function(x,
                                  y = NULL,
                                  paired = FALSE,
                                  alternative = c("two.sided",
                                                  "less",
                                                  "greater"),
                                  p_null = 0.5,
                                  alpha = 0.05,
                                  ...
) {
  alternative = match.arg(alternative)
  if(!missing(p_null) &&
     ((length(p_null) > 1L) || !is.finite(p_null)) &&
     p_null > 0 && p_null < 1)
    stop("'p_null' must be a single number between 0 and 1.")

  if(!is.numeric(x)) stop("'x' must be numeric")
  if(!is.null(y)) {
    if(!is.numeric(y)) stop("'y' must be numeric")
    DNAME <- paste(deparse(substitute(x)), "and",
                   deparse(substitute(y)))
    if(paired) {
      if(length(x) != length(y))
        stop("'x' and 'y' must have the same length")
      ok <- complete.cases(x,y)
      x = x[ok]
      y = y[ok]
    }
    else {
      x <- x[is.finite(x)]
      y <- y[is.finite(y)]
    }
  } else {

    stop("'y' is missing. One sample tests currently not supported.")

  }
  # Paired -----
  if(paired){

    METHOD = "Exact Paired Brunner-Munzel test"
    n = length(x)

    df.sw = n - 1
    all_data <- c(y, x)
    N = length(all_data)

    xinverse <- c(x, y)
    x1 <- y
    x2 <- x
    rx <- rank(all_data)
    rxinverse <- rank(xinverse)
    rx1 <- rx[1:n]
    rx2 <- rx[(n+1):N]
    rix1 <- rank(x1)
    rix2 <- rank(x2)
    BM1 <- 1 / n * (rx1 - rix1)
    BM2 <- 1 / n * (rx2 - rix2)
    BM3 <- BM1 - BM2
    BM4 <- 1 / (2 * n) * (rx1 - rx2)
    pd <- mean(BM2)

    m <- mean(BM3)
    v <- (sum(BM3 ^ 2) - n * m ^ 2) / (n - 1)
    v0 <- (v == 0)
    v[v0] <- 1 / n
    test_stat <- sqrt(n) * (pd - p_null) / sqrt(v)

    p.value = switch(alternative,
                     "two.sided" = (2*min(pt(test_stat,n-1),
                                          1-pt(test_stat,n-1))),
                     "less" = pt(test_stat,
                                 df=df.sw),
                     "greater" =  1-pt(test_stat,
                                       df=df.sw))

    pd.lower = switch(alternative,
                      "two.sided" = pd-qt(1-alpha/2,df.sw)*sqrt(v/n),
                      "less" = 0,
                      "greater" =  pd-qt(1-alpha,df.sw)*sqrt(v/n))

    pd.upper = switch(alternative,
                      "two.sided" = pd+qt(1-alpha/2,df.sw)*sqrt(v/n),
                      "less" = pd+qt(1-alpha,df.sw)*sqrt(v/n),
                      "greater" =  1)

  } else{

  # Two-sample ------
    rxy <- rank(c(x, y))
    rx <- rank(x)
    ry <- rank(y)
    n.x <- as.double(length(x))
    n.y <- as.double(length(y))
    N = n.x + n.y

    pl2 <- 1/n.y*(rxy[1:n.x]-rx)
    pl1 <- 1/n.x*(rxy[(n.x+1):N]-ry)
    pd <- mean(pl2)
    pd1 <- (pd == 1)
    pd0 <- (pd == 0)
    pd[pd1] <- 0.9999
    pd[pd0] <- 0.0001
    s1 <- var(pl2)/n.x
    s2 <- var(pl1)/n.y

    V <- N*(s1 +s2)
    singular.bf <- (V == 0)
    V[singular.bf] <- N/(2 * n.x * n.y)


    ## nparcomp ------
    test_stat <- sqrt(N)*(pd - p_null)/sqrt(V)
    df.sw <- (s1 + s2)^2/(s1^2/(n.x - 1) + s2^2/(n.y - 1))
    df.sw[is.nan(df.sw)] <- 1000
    METHOD = "Two Sample Brunner-Munzel test"

    p.value = switch(alternative,
                     "two.sided" = min(c(2 - 2 * pt(test_stat,
                                                  df=df.sw),
                                       2 * pt(test_stat,
                                              df=df.sw))),
                     "less" = pt(test_stat, df=df.sw),
                     "greater" =  1-pt(test_stat, df=df.sw))

    pd.lower = switch(alternative,
                      "two.sided" = pd - qt(1-alpha/2, df=df.sw)/sqrt(N)*sqrt(V),
                      "less" = 0,
                      "greater" =  pd - qt(1-alpha, df=df.sw)/sqrt(N)*sqrt(V))
    pd.lower = ifelse(pd.lower < 0, 0, pd.lower)

    pd.upper = switch(alternative,
                      "two.sided" = pd + qt(1-alpha/2, df=df.sw)/sqrt(N)*sqrt(V),
                      "less" = pd + qt(1-alpha, df=df.sw)/sqrt(N)*sqrt(V),
                      "greater" =  1)
    pd.upper = ifelse(pd.upper > 1, 1, pd.upper)
  }

  names(p_null) <- "probability of exchangeability"
  names(test_stat) = "t"
  names(df.sw) = "df"
  cint = c(pd.lower,pd.upper)
  attr(cint,"conf.level") = ifelse(alternative == "two.sided",
                                   1-alpha,
                                   1-2*alpha)
  estimate = pd
  names(estimate) = "p(X<Y) + .5*P(X=Y)"

  rval <- list(statistic = test_stat,
               parameter = df.sw,
               p.value = as.numeric(p.value),
               estimate = estimate,
               conf.int = cint,
               null.value = p_null,
               alternative = alternative,
               method = METHOD,
               data.name = DNAME)
  class(rval) <- "htest"
  return(rval)

}

#' @rdname brunner_munzel
#' @method brunner_munzel formula
#' @export

brunner_munzel.formula = function(formula,
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
  y <- do.call("brunner_munzel", c(DATA, list(...)))
  y$data.name <- DNAME
  y

}

