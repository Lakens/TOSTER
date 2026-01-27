#' @title Non-parametric standardized effect sizes (replicates of ses_calc)
#' @description
#' `r lifecycle::badge('superseded')`
#'
#' Effect sizes for simple (one or two sample) non-parametric tests. Suggested to use [ses_calc] function instead.
#' @param x a (non-empty) numeric vector of data values.
#' @param y an optional (non-empty) numeric vector of data values.
#' @param mu a number indicating the value around which (a-)symmetry (for
#'   one-sample or paired samples) or shift (for independent samples) is to be
#'   estimated. See [stats::wilcox.test].
#' @param paired a logical indicating whether you want to calculate a paired test.
#' @param ses Rank-biserial (rb), odds (odds), and concordance probability (cstat).
#' @param conf.level confidence level of the interval.
#'
#' @details
#' This method was adapted from the effectsize R package.
#' The rank-biserial correlation is appropriate for non-parametric tests of
#' differences - both for the one sample or paired samples case, that would
#' normally be tested with Wilcoxon's Signed Rank Test (giving the
#' **matched-pairs** rank-biserial correlation) and for two independent samples
#' case, that would normally be tested with Mann-Whitney's *U* Test (giving
#' **Glass'** rank-biserial correlation). See [stats::wilcox.test]. In both
#' cases, the correlation represents the difference between the proportion of
#' favorable and unfavorable pairs / signed ranks (Kerby, 2014). Values range
#' from `-1` indicating that all values of the second sample are smaller than
#' the first sample, to `+1` indicating that all values of the second sample are
#' larger than the first sample.
#'
#' In addition, the rank-biserial correlation can be transformed into a
#' concordance probability (i.e., probability of superiority) or into a
#' generalized odds (WMW odds or Agresti's generalized odds ratio).
#'
#' ## Ties
#' When tied values occur, they are each given the average of the ranks that
#' would have been given had no ties occurred. No other corrections have been
#' implemented yet.
#'
#' # Confidence Intervals
#' Confidence intervals for the standardized effect sizes
#' are estimated using the normal approximation (via Fisher's transformation).
#'
#' @return Returns a list of results including the rank biserial correlation, logical indicator if it was a paired method, setting for mu, and confidence interval.
#'
#' @references
#' - Cureton, E. E. (1956). Rank-biserial correlation. Psychometrika, 21(3),
#' 287-290.
#'
#' - Glass, G. V. (1965). A ranking variable analogue of biserial correlation:
#' Implications for short-cut item analysis. Journal of Educational Measurement,
#' 2(1), 91-95.
#'
#' - Kendall, M.G. (1948) Rank correlation methods. London: Griffin.
#'
#' - Kerby, D. S. (2014). The simple difference formula: An approach to teaching
#' nonparametric correlation. Comprehensive Psychology, 3, 11-IT.
#'
#' - King, B. M., & Minium, E. W. (2008). Statistical reasoning in the
#' behavioral sciences. John Wiley & Sons Inc.
#'
#' - Cliff, N. (1993). Dominance statistics: Ordinal analyses to answer ordinal
#' questions. Psychological bulletin, 114(3), 494.
#'
#' - Tomczak, M., & Tomczak, E. (2014). The need to report effect size estimates
#' revisited. An overview of some recommended measures of effect size.
#'
#' @export
#' @importFrom stats na.omit complete.cases qnorm

rbs <- function(x,
                y = NULL,
                mu = 0,
                conf.level = 0.95,
                paired = FALSE) {

 lifecycle::deprecate_warn(
    when = "0.9.0",
    what = "rbs()",
    with = "ses_calc()",
    details = "rbs() is deprecated. Please use ses_calc() which provides improved SE estimation methods and additional output options."
  )

  # adpated from
  # Ben-Shachar M, Lüdecke D, Makowski D (2020). effectsize: Estimation of Effect Size
  #   Indices and Standardized Parameters. Journal of Open Source Software, 5(56), 2815. doi:
  #   10.21105/joss.02815

  if (is.null(y)) {
    y <- rep(0, length.out = length(x))
    paired <- TRUE
  }

  if (paired) {
    i1 <- y
    i2 <- x
    data = data.frame(i1 = i1,
                      i2 = i2)
    data <- na.omit(data)
    colnames(data) = c("i1", "i2")
    x = data$i1
    y = data$i2

  } else {
    x <- na.omit(x)
    y <- na.omit(y)
  }

  ## Compute
  r_rbs <- rbs_calc(x, y,
                    mu = mu,
                    paired = paired)
  #out <- data.frame(r_rank_biserial = r_rbs)

  ## CI
  ci_method <- NULL
  if (is.numeric(conf.level)) {

    # Parametric method
    stopifnot(length(conf.level) == 1, conf.level < 1, conf.level > 0)
    alpha <- 1 - conf.level
    #out$CI <- conf.level
    rf <- atanh(r_rbs)
    if (paired) {
      nd <- length(x)
      maxw <- (nd^2 + nd) / 2

      # From: https://en.wikipedia.org/wiki/Wilcoxon_signed-rank_test#Historical_T_statistic
      # wSE <- sqrt((n * (n + 1) * (2 * n + 1)) / 24)
      # Delta method for f(x) = w * 2 / (maxw) - 1
      # r_rbsSE <- wSE * sqrt(4 / (maxw)^2)
      # Delta method for z: z_rbsSE <- r_rbsSE / (1 - r_rbs^2)
      #   But simulations suggest that z_rbsSE is positively biased
      #   more than r_rbsSE is negatively biased, especially when r_rbs is large,
      #   so we use r_rbsSE instead
      rfSE <- sqrt((2 * nd^3 + 3 * nd^2 + nd) / 6) / maxw
    } else {
      n1 <- length(x)
      n2 <- length(y)

      # From: https://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test#Normal_approximation_and_tie_correction
      # wSE <- sqrt((n1 * n2 * (n1 + n2 + 1)) / 12)
      # Delta method for f(x) = 1 - 2 * w / (n1 * n2) * sign(diff)
      # r_rbsSE <- wSE * sqrt(4 / (n1 * n2)^2)
      # Delta method for z: z_rbsSE <- r_rbsSE / (1 - r_rbs^2)
      #   But simulations suggest that z_rbsSE is positively biased
      #   more than r_rbsSE is negatively biased, especially when r_rbs is large,
      #   so we use r_rbsSE instead
      rfSE <- sqrt((n1 + n2 + 1) / (3 * n1 * n2))
    }

    confint <- tanh(rf + c(-1, 1) * qnorm(1 - alpha / 2) * rfSE)

  }

  rval = list(rbs = r_rbs,
              conf.int = confint,
              paired = paired,
              mu = mu)
  return(rval)
}

#' @rdname rbs
#' @export

np_ses <- function(x,
                   y = NULL,
                   mu = 0,
                   conf.level = 0.95,
                   paired = FALSE,
                   ses = c("rb","odds","logodds","cstat")) {

  lifecycle::deprecate_warn(
    when = "0.9.0",
    what = "np_ses()",
    with = "ses_calc()",
    details = "np_ses() is deprecated. Please use ses_calc() which provides improved SE estimation methods and additional output options."
  )

  ses = match.arg(ses)
  rb <- rbs(x=x,
            y = y,
            mu = mu,
            conf.level = conf.level,
            paired = paired)

  rb2 = switch(ses,
              "rb" = rb$rbs,
              "cstat" = rb_to_cstat(rb$rbs),
              "odds" = rb_to_odds(rb$rbs),
              "logodds" = log(rb_to_odds(rb$rbs)))
  confint = switch(ses,
              "rb" = rb$conf.int,
              "cstat" = rb_to_cstat(rb$conf.int),
              "odds" = rb_to_odds(rb$conf.int),
              "logodds" = log(rb_to_odds(rb$conf.int)))
  type = switch(ses,
                "rb" = "rb",
                "cstat" = "cstat",
                "odds" = "odds",
                "logodds" = "logodds")

  rval = list(
    type = ses,
    est = rb2,
    conf.int = confint,
    paired = paired,
    mu = mu
  )
  return(rval)
}
