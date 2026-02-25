#' @title Rank Difference Transformation for Paired Data
#' @description
#' `r lifecycle::badge('stable')`
#'
#' Applies the Kornbrot (1990) rank difference transformation to paired data.
#' All 2n observations are jointly ranked using midranks for ties, and the
#' ranks corresponding to each condition are returned. The transformed data
#' can then be passed to [ses_calc()] or [boot_ses_calc()]
#' for effect size estimation that is invariant under monotone transformations
#' of the original scale.
#'
#' @param x numeric vector of observations from condition 1.
#' @param y numeric vector of observations from condition 2, same length as x.
#'   Pairs are defined positionally: \code{x[i]} is paired with \code{y[i]}.
#' @param names optional character vector of length 2 giving column names for
#'   the returned data frame. Default is \code{c("x", "y")}.
#'
#' @details
#' The standard Wilcoxon signed-rank procedure for paired data computes
#' differences \eqn{d_i = x_i - y_i}, then ranks the absolute values
#' \eqn{|d_i|}. This is meaningful only when the differences themselves are
#' on an interval scale (i.e., when it makes sense to say that one difference
#' is "larger" than another).
#'
#' For purely ordinal data, the differences may not be rankable. The Kornbrot
#' (1990) rank difference procedure addresses this by:
#' \enumerate{
#'   \item Pooling all 2n observations from both conditions into a single vector.
#'   \item Ranking the pooled vector using standard midranks for ties.
#'   \item Returning the ranks corresponding to each condition.
#' }
#'
#' The resulting rank differences \eqn{R(x_i) - R(y_i)} are then suitable
#' for paired signed-rank effect size computation. Because the transformation
#' uses only the ordinal information in the data, the effect size is invariant
#' under any monotone (order-preserving) transformation of the original scale.
#'
#' ## Usage with ses_calc
#'
#' Pass the transformed columns directly to \code{ses_calc(..., paired = TRUE)}:
#'
#' \preformatted{
#'   rd <- rank_diff(x, y)
#'   ses_calc(x = rd$x, y = rd$y, paired = TRUE, ses = "rb")
#' }
#'
#' Because \code{mu} has no meaningful interpretation on the joint-rank scale,
#' always use \code{mu = 0} (the default) when analysing rank-difference data.
#'
#' @return A data frame with two columns (named by \code{names}) containing
#'   the joint ranks for condition 1 and condition 2, respectively. The number
#'   of rows equals \code{length(x)}. Missing-value pairs (where either
#'   \code{x[i]} or \code{y[i]} is \code{NA}) are removed before ranking,
#'   and a message is printed if any pairs are dropped.
#'
#' @examples
#' # Kornbrot (1990) Tables 1-2: time vs rate give different
#' # standard Wilcoxon results but identical rank difference results
#' time_plac <- c(4.6, 4.3, 6.7, 5.8, 5.0, 4.2, 6.0,
#'                2.0, 2.6, 10.0, 3.4, 7.1, 8.6)
#' time_drug <- c(2.9, 2.8, 12.0, 3.8, 5.9, 6.5, 3.3,
#'                2.3, 2.1, 14.3, 2.4, 14.0, 4.9)
#'
#' # Standard approach: different results for time vs rate
#' ses_calc(time_plac, time_drug, paired = TRUE, ses = "rb")
#' ses_calc(60 / time_plac, 60 / time_drug, paired = TRUE, ses = "rb")
#'
#' # Rank difference approach: identical results
#' rd_time <- rank_diff(time_plac, time_drug)
#' rd_rate <- rank_diff(60 / time_plac, 60 / time_drug)
#' ses_calc(rd_time$x, rd_time$y, paired = TRUE, ses = "rb")
#' ses_calc(rd_rate$x, rd_rate$y, paired = TRUE, ses = "rb")
#'
#' @references
#' Kornbrot, D. E. (1990). The rank difference test: A new and meaningful
#' alternative to the Wilcoxon signed ranks test for ordinal data.
#' *British Journal of Mathematical and Statistical Psychology*, 43, 241-264.
#'
#' @family effect sizes
#' @export
rank_diff <- function(x, y, names = c("x", "y")) {

  if (!is.numeric(x) || !is.numeric(y)) {
    stop("x and y must be numeric vectors.")
  }

  if (length(x) != length(y)) {
    stop("x and y must have the same length (paired data).")
  }

  if (length(names) != 2) {
    stop("names must be a character vector of length 2.")
  }

  # Remove pairs with any NA
  complete <- complete.cases(x, y)
  n_dropped <- sum(!complete)
  if (n_dropped > 0) {
    message(n_dropped, " pairs with missing values removed before ranking.")
    x <- x[complete]
    y <- y[complete]
  }

  if (length(x) == 0) {
    stop("No complete pairs remaining after removing missing values.")
  }

  # Pool all 2n observations and rank with midranks for ties
  n_pairs <- length(x)
  all_vals <- c(x, y)
  all_ranks <- rank(all_vals)

  # Return data frame with ranks for each condition
  out <- data.frame(
    all_ranks[1:n_pairs],
    all_ranks[(n_pairs + 1):(2 * n_pairs)]
  )
  colnames(out) <- names
  out
}
