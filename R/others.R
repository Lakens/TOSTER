

tost_decision = function(hypothesis = "EQU",
                         mu_text,
                         pvalue,
                         pTOST,
                         alpha){
  if (hypothesis == "EQU"){
    if(pvalue <= alpha && pTOST <= alpha){
      combined_outcome <- paste0("NHST: reject null significance hypothesis that the effect is equal to ",mu_text," \n",
                                 "TOST: reject null equivalence hypothesis")
    }
    if(pvalue < alpha && pTOST > alpha){
      combined_outcome <- paste0("NHST: reject null significance hypothesis that the effect is equal to ",mu_text," \n",
                                 "TOST: don't reject null equivalence hypothesis")
      # paste0("statistically different from ",mu_text," and not statistically equivalent")
    }
    if(pvalue > alpha && pTOST <= alpha){
      combined_outcome <- paste0("NHST: don't reject null significance hypothesis that the effect is equal to ",mu_text," \n",
                                 "TOST: reject null equivalence hypothesis")
      #paste0("statistically not different from ",mu_text," and statistically equivalent")
    }
    if(pvalue > alpha && pTOST > alpha){
      combined_outcome <- paste0("NHST: don't reject null significance hypothesis that the effect is equal to ",mu_text," \n",
                                 "TOST: don't reject null equivalence hypothesis")
      #paste0("statistically not different from ",mu_text," and not statistically equivalent")
    }
  } else {
    if(pvalue <= alpha && pTOST <= alpha){
      combined_outcome <- paste0("NHST: reject null significance hypothesis that the effect is equal to ",mu_text," \n",
                                 "TOST: reject null MET hypothesis")
      #paste0("statistically different from ",mu_text," and statistically greater than the minimal effect threshold")
    }
    if(pvalue < alpha && pTOST > alpha){
      combined_outcome <- paste0("NHST: reject null significance hypothesis that the effect is equal to ",mu_text," \n",
                                 "TOST: don't reject null MET hypothesis")
      #paste0("statistically different from ",mu_text," but not statistically greater than the minimal effect threshold")
    }
    if(pvalue > alpha && pTOST <= alpha){
      combined_outcome <- paste0("NHST: don't reject null significance hypothesis that the effect is equal to ",mu_text," \n",
                                 "TOST: reject null MET hypothesis")
      #paste0("statistically not different from ",mu_text," and statistically greater than the minimal effect threshold")
    }
    if(pvalue > alpha && pTOST > alpha){
      combined_outcome <- paste0("NHST: don't reject null significance hypothesis that the effect is equal to ",mu_text," \n",
                                 "TOST: don't reject null MET hypothesis")
      #paste0("statistically not different from ",mu_text," and not statistically greater than the minimal effect threshold")
    }
  }
  return(combined_outcome)
}

# Bootstrap CI functions ------

#' BCa bootstrap confidence interval (internal)
#'
#' Computes a bias-corrected and accelerated (BCa) bootstrap confidence interval.
#' The BCa method provides second-order accuracy by correcting for both bias and
#' skewness in the bootstrap distribution, using jackknife estimates to compute
#' the acceleration factor.
#'
#' @param boots_est Numeric vector of bootstrap estimates (on the working scale)
#' @param t0 Original estimate (on the same working scale as boots_est)
#' @param jack_est Numeric vector of jackknife (leave-one-out) estimates (on the same working scale)
#' @param alpha Significance level (e.g., 0.05 for 95% CI)
#' @return Numeric vector of length 2: c(lower, upper)
#' @keywords internal
bca_ci <- function(boots_est, t0, jack_est, alpha) {
  # Bias correction
  z0 <- qnorm(mean(boots_est < t0))

  # Check for infinite z0 (all boots on one side of t0)
  if (!is.finite(z0)) {
    stop(
      "BCa bias correction is infinite (all bootstrap estimates are on one side of ",
      "the original estimate). This may indicate a degenerate bootstrap distribution. ",
      "Consider using boot_ci = 'perc' instead.",
      call. = FALSE
    )
  }

  # Acceleration via jackknife
  L <- mean(jack_est) - jack_est
  denom <- 6 * sum(L^2)^(3/2)

  if (denom == 0) {
    stop(
      "BCa acceleration factor is degenerate (all jackknife estimates are identical). ",
      "This typically occurs with constant data or extreme boundary cases. ",
      "Consider using boot_ci = 'perc' instead.",
      call. = FALSE
    )
  }
  a <- sum(L^3) / denom

  # Adjusted quantiles
  z_alpha <- qnorm(c(alpha / 2, 1 - alpha / 2))
  numer <- z0 + z_alpha
  adj <- pnorm(z0 + numer / (1 - a * numer))

  if (any(adj <= 0) || any(adj >= 1)) {
    stop(
      "BCa adjusted quantiles are outside (0, 1), indicating extreme bias or skewness ",
      "in the bootstrap distribution. Consider using boot_ci = 'perc' instead.",
      call. = FALSE
    )
  }

  c(quantile(boots_est, adj[1], names = FALSE),
    quantile(boots_est, adj[2], names = FALSE))
}

basic <- function(boots_est, t0, alpha){
  conf = 1-alpha
  qq <- norm.inter(boots_est, (1 + c(conf, -conf))/2)
  c((2 *  t0 - qq[, 2L]))
}

perc <- function(boots_est, alpha = 0.05){
  conf.level = 1-alpha

  low <- (1 - conf.level)/2
  high <- 1 - low

  lower <- quantile(boots_est, low, names=FALSE)
  upper <- quantile(boots_est, high, names=FALSE)
  return(c(lower, upper))
}

stud <- function(boots_est, boots_se, se0, t0, alpha){
  conf = 1-alpha
  z <- (boots_est - t0)/(boots_se)
  qq <- norm.inter(z, (1 + c(conf, -conf))/2)
  c( ((t0 - (se0) * qq[, 2L])))

}

norm.inter = function (t, alpha) {
  t <- t[is.finite(t)]
  R <- length(t)
  rk <- (R + 1) * alpha
  if (!all(rk > 1 & rk < R))
    warning("extreme order statistics used as endpoints")
  k <- trunc(rk)
  inds <- seq_along(k)
  out <- inds
  kvs <- k[k > 0 & k < R]
  tstar <- sort(t, partial = sort(union(c(1, R), c(kvs, kvs +
                                                     1))))
  ints <- (k == rk)
  if (any(ints))
    out[inds[ints]] <- tstar[k[inds[ints]]]
  out[k == 0] <- tstar[1L]
  out[k == R] <- tstar[R]
  not <- function(v) xor(rep(TRUE, length(v)), v)
  temp <- inds[not(ints) & k != 0 & k != R]
  temp1 <- qnorm(alpha[temp])
  temp2 <- qnorm(k[temp]/(R + 1))
  temp3 <- qnorm((k[temp] + 1)/(R + 1))
  tk <- tstar[k[temp]]
  tk1 <- tstar[k[temp] + 1L]
  out[temp] <- tk + (temp1 - temp2)/(temp3 - temp2) * (tk1 -
                                                         tk)
  cbind(round(rk, 2), out)
}


# Function to test if an interval (defined by two numbers) does not contain zero
test_interval_no_zero <- function(vec) {
  # Check if vector has exactly 2 elements
  if (length(vec) != 2) {
    stop("Input must be a vector of exactly 2 numbers")
  }

  # Get the interval bounds
  lower <- min(vec)
  upper <- max(vec)

  # Check if zero is NOT in the interval [lower, upper]
  return(!(0 >= lower && 0 <= upper))
}

