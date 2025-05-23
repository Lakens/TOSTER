% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/brunner_munzel.R
\name{brunner_munzel}
\alias{brunner_munzel}
\alias{brunner_munzel.default}
\alias{brunner_munzel.formula}
\title{Brunner-Munzel Test}
\usage{
brunner_munzel(
  x,
  ...,
  paired = FALSE,
  alternative = c("two.sided", "less", "greater"),
  mu = 0.5,
  alpha = 0.05,
  perm = FALSE,
  max_n_perm = 10000
)

\method{brunner_munzel}{default}(
  x,
  y,
  paired = FALSE,
  alternative = c("two.sided", "less", "greater"),
  mu = 0.5,
  alpha = 0.05,
  perm = FALSE,
  max_n_perm = 10000,
  ...
)

\method{brunner_munzel}{formula}(formula, data, subset, na.action, ...)
}
\arguments{
\item{x}{a (non-empty) numeric vector of data values.}

\item{...}{further arguments to be passed to or from methods.}

\item{paired}{a logical indicating whether you want a paired test.}

\item{alternative}{a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less". You can specify just the initial letter.}

\item{mu}{a number specifying an optional parameter used to form the null hypothesis (Default = 0.5). This can be thought of as the null in terms of the relative effect, p = P (X < Y ) + 0.5 * P (X = Y); See ‘Details’.}

\item{alpha}{alpha level (default = 0.05)}

\item{perm}{a logical indicating whether or not to perform a permutation test over approximate t-distribution based test (default is FALSE). Highly recommend to set perm = TRUE when sample size per condition is less than 15.}

\item{max_n_perm}{the maximum number of permutations (default is 10000).}

\item{y}{an optional (non-empty) numeric vector of data values.}

\item{formula}{a formula of the form lhs ~ rhs where lhs is a numeric variable giving the data values and rhs either 1 for a one-sample or paired test or a factor with two levels giving the corresponding groups. If lhs is of class "Pair" and rhs is 1, a paired test is done.}

\item{data}{an optional matrix or data frame (or similar: see model.frame) containing the variables in the formula formula. By default the variables are taken from environment(formula).}

\item{subset}{an optional vector specifying a subset of observations to be used.}

\item{na.action}{a function which indicates what should happen when the data contain NAs. Defaults to getOption("na.action").}
}
\value{
A list with class \code{"htest"} containing the following components:
\itemize{
\item "statistic": the value of the test statistic.
\item "parameter": the degrees of freedom for the test statistic.
\item "p.value": the p-value for the test.
\item "conf.int": a confidence interval for the relative effect appropriate to the specified alternative hypothesis.
\item "estimate": the estimated relative effect.
\item "null.value": the specified hypothesized value of the relative effect.
\item "stderr": the standard error of the relative effect.
\item "alternative": a character string describing the alternative hypothesis.
\item "method": a character string indicating what type of test was performed.
\item "data.name": a character string giving the name(s) of the data.
}
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#maturing}{\figure{lifecycle-maturing.svg}{options: alt='[Maturing]'}}}{\strong{[Maturing]}}

This is a generic function that performs a generalized asymptotic Brunner-Munzel test in a fashion similar to \link{t.test}.
}
\details{
This function is made to provide a test of stochastic equality between two samples (paired or independent), and is referred to as the Brunner-Munzel test.

This tests the hypothesis that the relative effect, discussed below, is equal to the null value (default is \code{mu = 0.5}).

The estimate of the relative effect, which can be considered as value similar to the probability of superiority, refers to the following:

\deqn{\hat p = p(X>Y) + \frac{1}{2} \cdot P(X=Y)}

Note, for paired samples, this does \emph{not} refer to the probability of an increase/decrease in paired sample but rather the probability that a randomly sampled value of X.
This is also referred to as the "relative" effect in the literature. Therefore, the results will differ from the concordance probability provided by the ses_calc function.

The brunner_munzel function is based on the \code{npar.t.test} and \code{npar.t.test.paired} functions within the \code{nparcomp} package (Konietschke et al. 2015).
}
\examples{
data(mtcars)
brunner_munzel(mpg ~ am, data = mtcars)
}
\references{
Brunner, E., Munzel, U. (2000). The Nonparametric Behrens-Fisher Problem: Asymptotic Theory and a Small Sample Approximation. Biometrical Journal 42, 17 -25.

Neubert, K., Brunner, E., (2006). A Studentized Permutation Test for the Nonparametric Behrens-Fisher Problem. Computational Statistics and Data Analysis.

Munzel, U., Brunner, E. (2002). An Exact Paired Rank Test. Biometrical Journal 44, 584-593.

Konietschke, F., Placzek, M., Schaarschmidt, F., & Hothorn, L. A. (2015). nparcomp: an R software package for nonparametric multiple comparisons and simultaneous confidence intervals. Journal of Statistical Software 64 (2015), Nr. 9, 64(9), 1-17. http://www.jstatsoft.org/v64/i09/
}
\seealso{
Other Robust tests: 
\code{\link{boot_log_TOST}()},
\code{\link{boot_t_TOST}()},
\code{\link{boot_t_test}()},
\code{\link{log_TOST}()},
\code{\link{wilcox_TOST}()}
}
\concept{Robust tests}
