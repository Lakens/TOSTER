% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tsum_TOST.R
\name{tsum_TOST}
\alias{tsum_TOST}
\title{TOST with t-tests from Summary Statistics}
\usage{
tsum_TOST(
  m1,
  sd1,
  n1,
  m2 = NULL,
  sd2 = NULL,
  n2 = NULL,
  r12 = NULL,
  hypothesis = c("EQU", "MET"),
  paired = FALSE,
  var.equal = FALSE,
  eqb,
  low_eqbound,
  high_eqbound,
  mu = 0,
  eqbound_type = c("raw", "SMD"),
  alpha = 0.05,
  bias_correction = TRUE,
  rm_correction = FALSE,
  glass = NULL,
  smd_ci = c("nct", "goulet", "t", "z")
)
}
\arguments{
\item{m1}{mean of group 1.}

\item{sd1}{standard deviation of group 1.}

\item{n1}{sample size in group 1.}

\item{m2}{mean of group 2 (not required for one-sample tests).}

\item{sd2}{standard deviation of group 2 (not required for one-sample tests).}

\item{n2}{sample size in group 2 (not required for one-sample tests).}

\item{r12}{correlation between measurements for paired designs. Required when paired = TRUE.}

\item{hypothesis}{'EQU' for equivalence (default), or 'MET' for minimal effects test.}

\item{paired}{a logical indicating whether you want a paired t-test.}

\item{var.equal}{a logical variable indicating whether to treat the two variances as being equal. If TRUE then the pooled variance is used to estimate the variance otherwise the Welch (or Satterthwaite) approximation to the degrees of freedom is used.}

\item{eqb}{Equivalence bound. Can provide 1 value (symmetric bound, negative value is taken as the lower bound) or 2 specific values that represent the upper and lower equivalence bounds.}

\item{low_eqbound}{lower equivalence bounds (deprecated, use \code{eqb} instead).}

\item{high_eqbound}{upper equivalence bounds (deprecated, use \code{eqb} instead).}

\item{mu}{a number indicating the true value of the mean for the two-tailed test (or difference in means if you are performing a two sample test).}

\item{eqbound_type}{Type of equivalence bound. Can be 'SMD' for standardized mean difference (i.e., Cohen's d) or 'raw' for the mean difference. Default is 'raw'. Raw is strongly recommended as SMD bounds will produce biased results.}

\item{alpha}{alpha level (default = 0.05)}

\item{bias_correction}{Apply Hedges' correction for bias (default is TRUE).}

\item{rm_correction}{Repeated measures correction to make standardized mean difference Cohen's d(rm). This only applies to repeated/paired samples. Default is FALSE.}

\item{glass}{An option to calculate Glass's delta as an alternative to Cohen's d type SMD. Default is NULL to not calculate Glass's delta, 'glass1' will use the first group's SD as the denominator whereas 'glass2' will use the 2nd group's SD.}

\item{smd_ci}{Method for calculating SMD confidence intervals. Methods include 'goulet', 'noncentral t' (nct), 'central t' (t), and 'normal method' (z).}
}
\value{
An S3 object of class \code{"TOSTt"} is returned containing the following slots:
\itemize{
\item \strong{TOST}: A table of class \code{"data.frame"} containing two-tailed t-test and both one-tailed results.
\item \strong{eqb}: A table of class \code{"data.frame"} containing equivalence bound settings.
\item \strong{effsize}: Table of class \code{"data.frame"} containing effect size estimates.
\item \strong{hypothesis}: String stating the hypothesis being tested.
\item \strong{smd}: List containing the results of the standardized mean difference calculations (e.g., Cohen's d).
\itemize{
\item Items include: d (estimate), dlow (lower CI bound), dhigh (upper CI bound), d_df (degrees of freedom for SMD), d_sigma (SE), d_lambda (non-centrality), J (bias correction), smd_label (type of SMD), d_denom (denominator calculation).
}
\item \strong{alpha}: Alpha level set for the analysis.
\item \strong{method}: Type of t-test.
\item \strong{decision}: List included text regarding the decisions for statistical inference.
}
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}

Performs equivalence testing using the Two One-Sided Tests (TOST) procedure with t-tests
based on summary statistics rather than raw data. This function allows TOST analysis when
only descriptive statistics are available from published studies or reports.
}
\details{
This function performs TOST equivalence testing using summary statistics instead of raw data.
It is particularly useful when analyzing published results or conducting meta-analyses where
only summary statistics are available.

The function supports three types of tests:
\itemize{
\item One-sample test: Provide m1, sd1, and n1 only
\item Two-sample independent test: Provide all parameters except r12, with paired = FALSE
\item Paired samples test: Provide all parameters including r12, with paired = TRUE
}

For two-sample tests, the test is of \eqn{m1 - m2} (mean of group 1 minus mean of group 2).
For paired samples, the test is of the difference scores, wherein \eqn{z = m1 - m2}, and the test is of \eqn{\bar{z}} (mean of the difference scores).
For one-sample tests, the test is of \eqn{\bar{m1}} (mean of group 1).

The function calculates both raw mean differences and standardized effect sizes (Cohen's d or Hedges' g),
along with their confidence intervals.

For details on the calculations in this function see
\code{vignette("IntroTOSTt")} & \code{vignette("SMD_calcs")}.
}
\section{Purpose}{

Use this function when:
\itemize{
\item You only have access to summary statistics (means, standard deviations, sample sizes)
\item You want to perform meta-analyses using published results
\item You're conducting power analyses based on previous studies
\item You need to reanalyze published results within an equivalence testing framework
}
}

\examples{
# Example 1: One-sample test
# Testing if a sample with mean 0.55 and SD 4 (n=18) is equivalent to zero within ±2 units
tsum_TOST(m1 = 0.55, n1 = 18, sd1 = 4, eqb = 2)

# Example 2: Two-sample independent test
# Testing if two groups with different means are equivalent within ±3 units
tsum_TOST(m1 = 15.2, sd1 = 5.3, n1 = 30,
         m2 = 13.8, sd2 = 4.9, n2 = 28,
         eqb = 3)

# Example 3: Paired samples test
# Testing if pre-post difference is equivalent to zero within ±2.5 units
# with correlation between measurements of 0.7
tsum_TOST(m1 = 24.5, sd1 = 6.2, n1 = 25,
         m2 = 26.1, sd2 = 5.8, n2 = 25,
         r12 = 0.7, paired = TRUE,
         eqb = 2.5)

# Example 4: Two-sample test using standardized effect size bounds
# Testing if the standardized mean difference is within ±0.5 SD
tsum_TOST(m1 = 100, sd1 = 15, n1 = 40,
         m2 = 104, sd2 = 16, n2 = 42,
         eqb = 0.5, eqbound_type = "SMD")

}
\seealso{
Other TOST: 
\code{\link{boot_log_TOST}()},
\code{\link{boot_t_TOST}()},
\code{\link{simple_htest}()},
\code{\link{t_TOST}()},
\code{\link{wilcox_TOST}()}
}
\concept{TOST}
