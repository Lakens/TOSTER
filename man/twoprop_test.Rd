% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/two_proportions.R
\name{twoprop_test}
\alias{twoprop_test}
\title{Test of Proportions between 2 Independent Groups}
\usage{
twoprop_test(
  p1,
  p2,
  n1,
  n2,
  null = NULL,
  alpha = 0.05,
  alternative = c("two.sided", "less", "greater", "equivalence", "minimal.effect"),
  effect_size = c("difference", "odds.ratio", "risk.ratio")
)
}
\arguments{
\item{p1, p2}{Proportions in each respective group.}

\item{n1, n2}{sample size in each respective group.}

\item{null}{a number indicating the null hypothesis of the difference in proportions between two groups.}

\item{alpha}{alpha level (default = 0.05)}

\item{alternative}{a character string specifying the alternative hypothesis:
\itemize{
\item "two.sided": correlation is not equal to null (default)
\item "greater": correlation is greater than null
\item "less": correlation is less than null
\item "equivalence": correlation is within the equivalence bounds (TOST)
\item "minimal.effect": correlation is outside the equivalence bounds (TOST)
}

You can specify just the initial letter.}

\item{effect_size}{the effect size estimate, and confidence intervals, to calculate. Options include the difference between both proportions ("difference"), odds ratio ("odds.ratio"), or risk ratio ("risk.ratio").}
}
\value{
An S3 object of the class \code{htest}.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#maturing}{\figure{lifecycle-maturing.svg}{options: alt='[Maturing]'}}}{\strong{[Maturing]}}

This is a hypothesis testing function that mimics \link{prop.test}, but focuses only on testing differences in proportions between two groups.
This function utilizes a z-test to calculate the p-values (may be inaccurate with small sample sizes).
}
\details{
The hypothesis test for differences in proportions can be made on the raw proportions scale, the odds ratio, or the risk ratio (details below).
This function uses the large sample size asymptotic approximations for both the p-value and confidence interval calculations.
There should be a good deal of caution when sample sizes are small.
The p-values for the differences in proportions will differ from base \link{prop.test} due to the use of the unpooled standard error (see below).
\subsection{Differences in Proportions}{

Differences in proportions test is based on the following calculation:

\deqn{d = p_1 - p_2}

The standard error of \eqn{d} is calculated as the following:

\deqn{se(d) = \sqrt{\frac{p_1 \cdot (1-p_1)}{n_1} + \frac{p_2 \cdot (1-p_2)}{n_2}} }

The z-test, with \eqn{d_0} being the null value, is then calculated as the following (standard normal distribution evaluated to calculate p-value):

\deqn{z = \frac{d - d_0}{se(d)}}

The confidence interval can then be calculated as the following:

\deqn{d_{lower},d_{upper} = d \pm z_{\alpha} \cdot se(d)}
}

\subsection{Risk Ratio}{

The ratio between proportions test is based on the following calculation:

\deqn{\phi = p_1/p_2}

The standard error of \eqn{ln(\phi)} is calculated as the following:

\deqn{se(ln(\phi)) = \sqrt{\frac{1-p_1}{n_1 \cdot p_1} + \frac{1-p_2}{n_2 \cdot p_2}} }

The z-test, with \eqn{\phi_0} being the null value, is then calculated as the following (standard normal distribution evaluated to calculate p-value):

\deqn{z = \frac{ln(\phi) - ln(\phi_0)}{se(ln(\phi))}}

The confidence interval can then be calculated as the following:

\deqn{\phi_{lower} = \phi \cdot e^{-z_{\alpha} \cdot se(ln(\phi))}}

\deqn{\phi_{upper} = \phi \cdot e^{z_{\alpha} \cdot se(ln(\phi))}}
}

\subsection{Odds Ratio}{

The ratio between proportions test is based on the following calculation:
(p1/q1) / (p2/q2)

\deqn{OR = \frac{p_1}{1-p_1} / \frac{p_2}{1-p_2}}

The standard error of \eqn{ln(OR)} is calculated as the following:

\deqn{se(ln(OR)) = \sqrt{\frac{1}{n_1 \cdot p_1 + 0.5} + \frac{1}{n_1 \cdot (1-p_1) + 0.5} + \frac{1}{n_2 \cdot p_2 + 0.5} + \frac{1}{n_2 \cdot (1-p_2) + 0.5} } }

The z-test, with \eqn{OR_0} being the null value, is then calculated as the following (standard normal distribution evaluated to calculate p-value):

\deqn{z = \frac{ln(OR) - ln(OR_0)}{se(ln(OR))}}

The confidence interval can then be calculated as the following:

\deqn{OR_{lower},OR_{upper} = exp(ln(OR) \pm z_{\alpha} \cdot se(ln(OR)))}
}
}
\references{
Gart, J. J., & Nam, J. M. (1988). Approximate interval estimation of the ratio of binomial parameters: a review and corrections for skewness. Biometrics, 323-338.

Tunes da Silva, G., Logan, B. R., & Klein, J. P. (2008). Methods for Equivalence and Noninferiority Testing. Biology of Blood Marrow Transplant, 15(1 Suppl), 120-127.

Yin, G. (2012). Clinical Trial Design: Bayesian and Frequentist Adaptive Methods. Hoboken, New Jersey: John Wiley & Sons, Inc.
}
