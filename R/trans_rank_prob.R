#' @title Rescale a Probability-Scale Effect Size
#' @description
#' `r lifecycle::badge('stable')`
#'
#' Transforms a probability-scale effect size (and, optionally, its standard
#' error, confidence interval, and null value) between four scales:
#' `"probability"`, `"difference"`, `"logodds"`, and `"odds"`.
#'
#' This function serves as both a standalone utility and the internal engine
#' for `brunner_munzel(scale = ...)`.
#'
#' @param estimate numeric; the point estimate to transform.
#' @param se numeric or `NULL`; standard error of `estimate`.
#'   Transformed via the delta method.
#' @param ci numeric vector of length 2 or `NULL`; confidence interval
#'   endpoints.
#'   Because every scale conversion is monotonic, CI endpoints are transformed
#'   directly (coverage is preserved without the delta method).
#' @param null numeric (scalar or vector) or `NULL`; the null-hypothesis
#'   value(s) to transform (e.g., a single null, or two equivalence bounds).
#' @param from character; the scale `estimate` is currently on.
#'   One of `"probability"`, `"difference"`, `"logodds"`, `"odds"`.
#' @param to character; the target scale.
#'   One of `"probability"`, `"difference"`, `"logodds"`, `"odds"`.
#'
#' @details
#' The four scales and their relationship to a probability \eqn{p} are:
#'
#' | Scale | Domain | Formula | Null at stochastic equality |
#' |-------|--------|---------|-----------------------------|
#' | probability | \eqn{(0, 1)} | \eqn{p} | 0.5 |
#' | difference  | \eqn{(-1, 1)} | \eqn{2p - 1} | 0 |
#' | logodds     | \eqn{(-\infty, \infty)} | \eqn{\log[p / (1-p)]} | 0 |
#' | odds        | \eqn{(0, \infty)} | \eqn{p / (1-p)} | 1 |
#'
#' All conversions are routed through the probability scale internally.
#'
#' **Standard errors** are transformed via the delta method:
#'
#' \deqn{\mathrm{SE}_{\mathrm{target}} =
#'   \mathrm{SE}_{\mathrm{original}} \times
#'   \left|\frac{dp}{dx}\right| \times
#'   \left|\frac{dy}{dp}\right|}
#'
#' where \eqn{x} is the original scale and \eqn{y} is the target scale.
#'
#' **Confidence intervals** are transformed by applying the monotonic mapping
#' directly to each endpoint, which preserves coverage exactly.
#'
#' At the boundaries (\eqn{p = 0} or \eqn{p = 1}), transformations to the
#' logodds scale return \eqn{\pm\infty}, and the delta-method SE is infinite.
#' This is mathematically correct behaviour; no clamping or warning is applied.
#'
#' @return A list with components:
#'   \describe{
#'     \item{estimate}{transformed point estimate}
#'     \item{se}{transformed standard error (or `NULL`)}
#'     \item{ci}{transformed CI endpoints (or `NULL`)}
#'     \item{null}{transformed null value(s) (or `NULL`)}
#'     \item{from}{the `from` scale (echoed back)}
#'     \item{to}{the `to` scale (echoed back)}
#'   }
#'
#' @examples
#' # Probability to difference (rank-biserial)
#' trans_rank_prob(0.7, se = 0.05, ci = c(0.6, 0.8),
#'                null = 0.5, from = "probability", to = "difference")
#'
#' # Probability to odds
#' trans_rank_prob(0.7, from = "probability", to = "odds")
#'
#' # Round-trip: logodds -> probability -> logodds
#' lo <- trans_rank_prob(0.8473, from = "logodds", to = "probability")
#' trans_rank_prob(lo$estimate, from = "probability", to = "logodds")
#'
#' # Apply to brunner_munzel output
#' res <- brunner_munzel(mpg ~ am, data = mtcars)
#' trans_rank_prob(as.numeric(res$estimate),
#'                se = res$stderr,
#'                ci = as.numeric(res$conf.int),
#'                null = as.numeric(res$null.value),
#'                from = "probability", to = "logodds")
#'
#' @seealso [brunner_munzel()], [ses_calc()]
#' @family effect sizes
#' @export
trans_rank_prob <- function(estimate,
                            se = NULL,
                            ci = NULL,
                            null = NULL,
                            from = c("probability", "difference",
                                     "logodds", "odds"),
                            to = c("probability", "difference",
                                   "logodds", "odds")) {

  from <- match.arg(from)
  to <- match.arg(to)

  # Identity: return early
  if (from == to) {
    return(list(estimate = estimate, se = se, ci = ci,
                null = null, from = from, to = to))
  }

  # Step 1: Convert to probability scale
  p      <- to_probability(estimate, from)
  p_ci   <- if (!is.null(ci))   to_probability(ci, from)   else NULL
  p_null <- if (!is.null(null)) to_probability(null, from) else NULL
  p_se   <- if (!is.null(se))   se * abs(d_to_probability(estimate, from)) else NULL

  # Step 2: Convert from probability to target scale
  est_out  <- from_probability(p, to)
  ci_out   <- if (!is.null(p_ci))   from_probability(p_ci, to)   else NULL
  null_out <- if (!is.null(p_null)) from_probability(p_null, to) else NULL
  se_out   <- if (!is.null(p_se))   p_se * abs(d_from_probability(p, to)) else NULL

  list(estimate = est_out, se = se_out, ci = ci_out,
       null = null_out, from = from, to = to)
}

# Internal helpers --------
# These are non-exported internal functions.

# Transform a value on `scale` to the probability scale
to_probability <- function(x, scale) {
  switch(scale,
    "probability" = x,
    "difference"  = (x + 1) / 2,
    "logodds"     = plogis(x),
    "odds"        = x / (1 + x)
  )
}

# Transform a value on the probability scale to `scale`
from_probability <- function(p, scale) {
  switch(scale,
    "probability" = p,
    "difference"  = 2 * p - 1,
    "logodds"     = qlogis(p),
    "odds"        = p / (1 - p)
  )
}

# Derivative of to_probability(x, scale) w.r.t. x
# se_p = se_x * |d/dx to_probability(x)|
d_to_probability <- function(x, scale) {
  switch(scale,
    "probability" = 1,
    "difference"  = 0.5,
    "logodds"     = dlogis(x),
    "odds"        = 1 / (1 + x)^2
  )
}

# Derivative of from_probability(p, scale) w.r.t. p
# se_target = se_p * |d/dp from_probability(p)|
d_from_probability <- function(p, scale) {
  switch(scale,
    "probability" = 1,
    "difference"  = 2,
    "logodds"     = 1 / (p * (1 - p)),
    "odds"        = 1 / (1 - p)^2
  )
}

# Build probability-notation estimate labels --------
#
# Constructs the probability-notation label for an estimate on a given scale.
# Used by brunner_munzel() and ses_calc() to produce consistent labels.
#
# @param scale character; one of "probability", "difference", "logodds", "odds"
# @param xname character; quoted name for the first group/variable
# @param yname character; quoted name for the second group/variable (or NULL for one-sample)
# @param paired logical; whether the comparison is paired
# @return character string label
prob_notation_label <- function(scale, xname, yname = NULL, paired = FALSE,
                                paired_style = c("direct", "difference")) {
  paired_style <- match.arg(paired_style)

  # Quote names only when they are numeric
  quote_if_numeric <- function(nm) {
    if (grepl("^-?\\d*\\.?\\d+$", nm)) paste0("'", nm, "'") else nm
  }

  xq <- quote_if_numeric(xname)

  if (is.null(yname)) {
    # One-sample: compare against 0
    prob_label <- paste0("P(", xq, ">0) + .5*P(", xq, "=0)")
    diff_label <- paste0("P(", xq, ">0) - P(", xq, "<0)")
  } else if (paired && paired_style == "difference") {
    # Paired with difference-score notation: P(X - Y>0)
    yq <- quote_if_numeric(yname)
    d_expr <- paste0(xq, " - ", yq)
    prob_label <- paste0("P(", d_expr, ">0) + .5*P(", d_expr, "=0)")
    diff_label <- paste0("P(", d_expr, ">0) - P(", d_expr, "<0)")
  } else {
    # Two-sample independent (or paired with direct notation): P(X>Y)
    yq <- quote_if_numeric(yname)
    prob_label <- paste0("P(", xq, ">", yq, ") + .5*P(", xq, "=", yq, ")")
    diff_label <- paste0("P(", xq, ">", yq, ") - P(", xq, "<", yq, ")")
  }

  switch(scale,
    "probability" = prob_label,
    "difference"  = diff_label,
    "logodds"     = paste0("logodds(", prob_label, ")"),
    "odds"        = paste0("odds(", prob_label, ")")
  )
}
