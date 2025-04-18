% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/power_t_TOST.R
\name{power_t_TOST}
\alias{power_t_TOST}
\title{Power calculations for TOST with t-tests}
\usage{
power_t_TOST(
  n = NULL,
  delta = 0,
  sd = 1,
  eqb,
  low_eqbound = NULL,
  high_eqbound = NULL,
  alpha = NULL,
  power = NULL,
  type = "two.sample"
)
}
\arguments{
\item{n}{number of observations per group. 2 sample sizes, in a vector, can be provided for the two sample case.}

\item{delta}{true difference in means (default is 0).}

\item{sd}{population standard deviation. Standard deviation of the differences for paired samples.}

\item{eqb}{Equivalence bound. Can provide 1 value (negative value is taken as the lower bound) or 2 specific values that represent the upper and lower equivalence bounds.}

\item{low_eqbound}{Lower equivalence bounds. Deprecated use eqb.}

\item{high_eqbound}{Upper equivalence bounds. Deprecated use eqb.}

\item{alpha}{a priori alpha-level (i.e., significance level).}

\item{power}{power of the TOST procedure (1-beta).}

\item{type}{string specifying the type of t-test.}
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}

Calculates the exact power of two one sided t-tests (TOST) for one, two, and paired samples.
}
\details{
The exact calculations of power are based on Owen’s Q-function or by direct integration of the bivariate non-central t-distribution (inspired by the PowerTOST package).
Approximate power is implemented via the non-central t-distribution or the ‘shifted’ central t-distribution.
}
\note{
The power function in this package is limited. Please see the PowerTOST R package for more options.
}
\references{
Phillips KF. Power of the Two One-Sided Tests Procedure in Bioequivalence. J Pharmacokin Biopharm. 1990;18(2):137–44. doi: 10.1007/BF01063556

Diletti D, Hauschke D, Steinijans VW. Sample Size Determination for Bioequivalence Assessment by Means of Confidence Intervals. Int J Clin Pharmacol Ther Toxicol. 1991;29(1):1–8.
}
\seealso{
Other power: 
\code{\link{power_eq_f}()},
\code{\link{power_z_cor}()}
}
\concept{power}
