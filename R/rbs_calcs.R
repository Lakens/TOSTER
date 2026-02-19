rbs_calc = function (x, y,
                     mu, paired) {
  # adapted from effectsize R package
  if (paired) {
    z <- (x-y) - mu
    z <- z[z != 0]
    if (length(z) == 0) return(0)
    abs_z = abs(z)
    RR = -1 * rank(abs_z) * sign(z)
    Rplus = sum(RR[RR > 0])
    Rminus = sum(abs(RR[RR < 0]))
    Tee = min(Rplus, Rminus)
    n = length(RR)
    if (Rplus >= Rminus) {
      rho = -4 * abs((Tee - (Rplus + Rminus)/2)/n/(n + 1))
    }
    if (Rplus < Rminus) {
      rho = 4 * abs((Tee - (Rplus + Rminus)/2)/n/(n + 1))
    }
    return(rho)
  }
  else {
    Ry <- ranktransform(c(x - mu, y))
    n1 <- length(x)
    n2 <- length(y)
    S <- (n1 * n2)
    U1 <- sum(Ry[seq_along(x)]) - n1 * (n1 + 1)/2
    U2 <- sum(Ry[-seq_along(x)]) - n2 * (n2 + 1)/2
    u_ <- U1/S
    f_ <- U2/S
    return(u_ - f_)
  }

}

ranktransform <- function(x,
                          sign = FALSE,
                          method = "average") {


  # Warning if all NaNs
  if (all(is.na(x))) {
    return(x)
  }

  # Warning if only one value
  if (length(unique(x)) == 1) {
    if (is.null(names(x))) {
      name <- deparse(substitute(x))
    } else {
      name <- names(x)
    }

    return(x)
  }


  # Warning if logical vector
  if (length(unique(x)) == 2) {
    if (is.null(names(x))) {
      name <- deparse(substitute(x))
    } else {
      name <- names(x)
    }

  }


  if (sign) {
    ZEROES <- x == 0
    #if (any(ZEROES) && verbose) warning("Zeros detected. These cannot be sign-rank transformed.")
    out <- rep(NA, length(x))
    out[!ZEROES] <- sign(x[!ZEROES]) * rank(abs(x[!ZEROES]),
                                            ties.method = method,
                                            na.last = "keep")
  } else {
    out <- rank(x, ties.method = method, na.last = "keep")
  }

  return(out)
}

odds_to_pr <- function(x, log = FALSE) {
  if (log) {
    stats::plogis(x)
  } else {
    stats::plogis(log(x))
  }
}

pr_to_odds <- function(x, log = FALSE) {
  if (log) {
    stats::qlogis(x)
  } else {
    exp(stats::qlogis(x))
  }
}

rb_to_odds <- function(x) {
  pr_to_odds(rb_to_cstat(x))
}

rb_to_cstat <- function(x) {
  (x + 1) / 2
}

cstat_to_rb <- function(x){
  2*x-1
}

z_to_rho <- function(x){

  tanh(x)
}

rho_to_z <- function(x){
  atanh(x)
}

# === Transformation helper functions ===

#' Transform effect size values to log-odds scale
#' @param value numeric value on the specified scale
#' @param scale character string specifying the scale ("rb", "cstat", "odds", "logodds")
#' @return value on the log-odds scale
#' @noRd
to_logodds <- function(value, scale) {

  switch(scale,
    "rb" = {
      p <- (value + 1) / 2
      log(p / (1 - p))
    },
    "cstat" = {
      log(value / (1 - value))
    },
    "odds" = {
      log(value)
    },
    "logodds" = {
      value
    }
  )
}

#' Convert effect size value to cstat (concordance probability) scale
#' @param value the value to convert
#' @param from_scale character: "rb", "cstat", "odds", or "logodds"
#' @return value on cstat scale
#' @noRd
to_cstat <- function(value, from_scale) {
  switch(from_scale,
    "cstat" = value,
    "rb" = (value + 1) / 2,
    "odds" = value / (1 + value),
    "logodds" = plogis(value)
  )
}

# === Score-type functions for paired/one-sample designs ===

#' Extract paired rank information for score-based inference
#'
#' Computes the key quantities from paired/one-sample data needed for
#' score-based inference: N (non-zero differences), S, Q, T+, p_hat.
#'
#' @param x numeric vector (first sample or one-sample data)
#' @param y numeric vector (second sample), or NULL for one-sample
#' @param mu hypothesized difference (default 0)
#' @return list with n_eff, S, Q, T_plus, p_hat
#' @noRd
paired_rank_info <- function(x, y = NULL, mu = 0) {
  if (is.null(y)) {
    d <- x - mu
  } else {
    d <- x - y - mu
  }

  # Remove zeros (ties with mu)
  d <- d[d != 0]
  N <- length(d)

  if (N == 0) {
    return(list(n_eff = 0, S = 0, Q = 0, T_plus = 0, p_hat = 0.5))
  }

  # Rank absolute differences (average ties)
  r <- rank(abs(d))

  # T+ = sum of ranks where d > 0
  T_plus <- sum(r[d > 0])

  S <- N * (N + 1) / 2
  Q <- sum(r^2)
  p_hat <- T_plus / S

  list(
    n_eff = N,
    S     = S,
    Q     = Q,
    T_plus = T_plus,
    p_hat  = p_hat
  )
}

#' Score CI for paired designs (Wilson-type)
#'
#' Computes a Wilson-score-type confidence interval for the concordance
#' probability (cstat) from paired/one-sample signed-rank data.
#'
#' The CI inverts the score test, yielding a quadratic in pi0 with
#' closed-form solution (no root-finding needed).
#'
#' @param p_hat concordance estimate (T_plus / S)
#' @param n_eff number of non-zero differences
#' @param Q sum of squared (average) ranks
#' @param conf.level confidence level
#' @param correct logical; apply continuity correction?
#' @return numeric vector of length 2: (lower, upper) on cstat scale
#' @noRd
score_ci_paired <- function(p_hat, n_eff, Q,
                            conf.level = 0.95,
                            correct = FALSE) {

  S <- n_eff * (n_eff + 1) / 2

  if (n_eff < 1 || Q == 0) {
    return(c(0, 1))
  }

  alpha <- 1 - conf.level
  z_crit <- qnorm(1 - alpha / 2)
  z2 <- z_crit^2

  # Ratio that appears in Wilson formula
  c_val <- z2 * Q / S^2

  # Helper: solve the Wilson quadratic for a given "center" a
  # Returns (lower_root, upper_root)
  wilson_roots <- function(a) {
    A <- 1 + c_val
    B <- -(2 * a + c_val)
    C <- a^2

    disc <- B^2 - 4 * A * C
    if (disc < 0) disc <- 0

    lo <- (-B - sqrt(disc)) / (2 * A)
    hi <- (-B + sqrt(disc)) / (2 * A)
    c(lo, hi)
  }

  if (!correct) {
    roots <- wilson_roots(p_hat)
    lower <- roots[1]
    upper <- roots[2]
  } else {
    # Continuity correction: shift p_hat by h = 0.5/S toward pi0
    # Lower bound: use a' = min(p_hat + h, 1), take smaller root
    # Upper bound: use a'' = max(p_hat - h, 0), take larger root
    h <- 0.5 / S

    roots_low  <- wilson_roots(min(p_hat + h, 1))
    roots_high <- wilson_roots(max(p_hat - h, 0))

    lower <- roots_low[1]
    upper <- roots_high[2]
  }

  # Clamp to [0, 1]
  lower <- max(0, lower)
  upper <- min(1, upper)

  c(lower, upper)
}

#' Score p-value for paired designs
#'
#' Computes a score-type p-value for testing H0: pi = pi0,
#' where pi is the concordance probability.
#' At pi0 = 0.5 this matches wilcox.test(..., exact = FALSE).
#'
#' @param p_hat concordance estimate
#' @param pi0 null hypothesis value on cstat scale
#' @param n_eff number of non-zero differences
#' @param Q sum of squared (average) ranks
#' @param alternative "two.sided", "less", or "greater"
#' @param correct apply continuity correction?
#' @return list with z.statistic and p.value
#' @noRd
score_pvalue_paired <- function(p_hat, pi0, n_eff, Q,
                                alternative = "two.sided",
                                correct = FALSE) {

  S <- n_eff * (n_eff + 1) / 2

  if (n_eff < 1 || Q == 0) {
    return(list(z.statistic = NA_real_, p.value = NA_real_))
  }

  # Variance under H0
  var_null <- pi0 * (1 - pi0) * Q / S^2
  se_null  <- sqrt(var_null)

  if (se_null < .Machine$double.eps) {
    return(list(z.statistic = NA_real_, p.value = NA_real_))
  }

  if (correct) {
    h <- 0.5 / S
    numer <- max(abs(p_hat - pi0) - h, 0)
    z <- numer / se_null * sign(p_hat - pi0)
    # Edge case: if p_hat == pi0 after correction, z = 0
    if (p_hat == pi0) z <- 0
  } else {
    z <- (p_hat - pi0) / se_null
  }

  p_val <- switch(alternative,
    "two.sided" = 2 * pnorm(-abs(z)),
    "less"      = pnorm(z),
    "greater"   = pnorm(z, lower.tail = FALSE)
  )

  list(z.statistic = z, p.value = p_val)
}

#' Descriptive SE for paired score method
#'
#' Computes descriptive standard errors evaluated at p_hat.
#' These are for reporting; the CI comes from test inversion, not SE +/- z.
#'
#' @param p_hat concordance estimate
#' @param n_eff number of non-zero differences
#' @param Q sum of squared (average) ranks
#' @return list with se_cstat, se_rb, se_logodds, se_odds
#' @noRd
score_se_paired <- function(p_hat, n_eff, Q) {

  S <- n_eff * (n_eff + 1) / 2

  # Boundary-safe p_hat for SE computation
  p_safe <- pmin(pmax(p_hat, 1e-10), 1 - 1e-10)

  se_cstat   <- sqrt(p_safe * (1 - p_safe) * Q / S^2)
  se_rb      <- 2 * se_cstat
  se_logodds <- se_cstat / (p_safe * (1 - p_safe))
  se_odds    <- se_cstat / (1 - p_safe)^2

  list(
    se_cstat   = se_cstat,
    se_rb      = se_rb,
    se_logodds = se_logodds,
    se_odds    = se_odds
  )
}

# === Score-type CI functions (Fay & Malinovsky 2018) ===

#' Compute tie adjustment factor for WMW variance
#'
#' Computes the tie correction factor used in the variance of the
#' Wilcoxon-Mann-Whitney statistic.
#'
#' @param x numeric vector, group 1
#' @param y numeric vector, group 2
#' @return scalar tie factor in [0, 1]; equals 1 when no ties
#' @noRd
wmw_tie_factor <- function(x, y) {
  r <- rank(c(x, y))
  N <- length(r)
  NTIES <- table(r)
  tf <- 1 - sum(NTIES^3 - NTIES) / (N * (N + 1) * (N - 1))
  tf
}

#' Score-type variance of the Mann-Whitney parameter (V_LAPH)
#'
#' Computes V_LAPH(phi) from Fay & Malinovsky (2018), which is the
#' variance of the Mann-Whitney estimator under the proportional odds
#' model, evaluated at a candidate parameter value phi.
#'
#' @param phi candidate value of the Mann-Whitney parameter (scalar in (0,1))
#' @param tf tie adjustment factor from wmw_tie_factor()
#' @param n1 sample size for group 1 (x)
#' @param n2 sample size for group 2 (y)
#' @return scalar variance
#' @references
#' Fay, M.P. and Malinovsky, Y. (2018). Confidence Intervals of the
#' Mann-Whitney Parameter that are Compatible with the Wilcoxon-Mann-Whitney
#' Test. Statistics in Medicine, 37, 3991-4006.
#' @noRd
v_laph <- function(phi, tf, n1, n2) {
  tf * (phi * (1 - phi) / (n1 * n2)) *
    (1 + ((n1 + n2 - 2) / 2) *
       ((1 - phi) / (2 - phi) + phi / (1 + phi)))
}

#' Score-type p-value for the Mann-Whitney parameter
#'
#' Computes a p-value testing H0: phi = phi_null using the score statistic
#' with optional continuity correction.
#'
#' @param phi_hat observed Mann-Whitney estimate (concordance probability)
#' @param phi_null null hypothesis value on the cstat scale
#' @param tf tie adjustment factor
#' @param n1 sample size for group 1
#' @param n2 sample size for group 2
#' @param alternative character: "two.sided", "less", or "greater"
#' @param correct logical: apply continuity correction?
#' @return list with elements: p.value, z.statistic
#' @noRd
score_pvalue_wmw <- function(phi_hat, phi_null, tf, n1, n2,
                             alternative = "two.sided", correct = FALSE) {

  # Continuity corrections (separate for each direction)
  corr_less <- 0
  corr_greater <- 0
  if (correct) {
    corr_less <- -0.5 / (n1 * n2)
    corr_greater <- 0.5 / (n1 * n2)
  }

  V_null <- v_laph(phi_null, tf, n1, n2)

  # Handle degenerate case (all values equal -> tf = 0 -> V = 0)
  if (V_null <= 0) {
    return(list(p.value = 1, z.statistic = 0))
  }

  se_null <- sqrt(V_null)

  z_less <- (phi_hat - phi_null - corr_less) / se_null
  z_greater <- (phi_hat - phi_null - corr_greater) / se_null

  p_val <- switch(alternative,
    "two.sided" = 2 * min(pnorm(z_less), pnorm(z_greater, lower.tail = FALSE)),
    "less" = pnorm(z_less),
    "greater" = pnorm(z_greater, lower.tail = FALSE)
  )

  # For reporting, use the z closer to the observed direction
  z_stat <- switch(alternative,
    "two.sided" = if (phi_hat >= phi_null) z_greater else z_less,
    "less" = z_less,
    "greater" = z_greater
  )

  list(p.value = min(1, p_val), z.statistic = z_stat)
}

#' Score-type confidence interval for the Mann-Whitney parameter
#'
#' Constructs a CI by inverting the score test: finds the values of phi
#' where the score statistic equals the critical value.
#'
#' @param phi_hat observed Mann-Whitney estimate (concordance probability)
#' @param tf tie adjustment factor
#' @param n1 sample size for group 1
#' @param n2 sample size for group 2
#' @param conf.level confidence level (e.g. 0.95)
#' @param correct logical: apply continuity correction?
#' @param epsilon small number for uniroot bounds (default 1e-8)
#' @return numeric vector of length 2: c(lower, upper) on cstat scale
#' @noRd
score_ci_wmw <- function(phi_hat, tf, n1, n2, conf.level = 0.95,
                         correct = FALSE, epsilon = 1e-8) {

  alpha <- 1 - conf.level

  # Continuity corrections
  corr_less <- 0
  corr_greater <- 0
  if (correct) {
    corr_less <- -0.5 / (n1 * n2)
    corr_greater <- 0.5 / (n1 * n2)
  }

  # The score function: at the true phi_null, this should equal zq
  wfunc <- function(phi_null, zq, correction = 0) {
    V <- v_laph(phi_null, tf, n1, n2)
    if (V <= 0) return(-Inf)
    (phi_hat - phi_null - correction) / sqrt(V) - zq
  }

  phimin <- epsilon
  phimax <- 1 - epsilon

  root <- function(zq, corr) {
    f.lower <- wfunc(phimin, zq, corr)
    if (f.lower <= 0) return(phimin)
    f.upper <- wfunc(phimax, zq, corr)
    if (f.upper >= 0) return(phimax)
    uniroot(wfunc, c(phimin, phimax),
            f.lower = f.lower, f.upper = f.upper,
            tol = epsilon, zq = zq, correction = corr)$root
  }

  # Two-sided CI
  # Lower bound: find phi_L where score stat = z_{1-alpha/2}
  # Upper bound: find phi_U where score stat = z_{alpha/2}

  lower <- if (phi_hat == 0) {
    0
  } else {
    root(zq = qnorm(alpha / 2, lower.tail = FALSE), corr = corr_greater)
  }

  upper <- if (phi_hat == 1) {
    1
  } else {
    root(zq = qnorm(alpha / 2), corr = corr_less)
  }

  c(lower, upper)
}

# === New variance estimation functions (Agresti/Lehmann method) ===

#' Compute placement values for two-sample comparison
#' @param x numeric vector for group 1
#' @param y numeric vector for group 2
#' @return list with V (placements for x), W (placements for y), p_hat, n1, n2
#' @details
#' Note: p_hat is computed as Pr(X > Y) to match the sign convention used by rbs_calc,
#' where positive rank-biserial indicates x values tend to be larger than y values.
#' @noRd
compute_placements <- function(x, y) {
  n1 <- length(x)
  n2 <- length(y)

  # For each x: proportion of y's that are less (x dominates) + 0.5 * ties
  # This gives Pr(X > Y | X = xi)
  V <- sapply(x, function(xi) {
    mean(y < xi) + 0.5 * mean(y == xi)
  })

  # For each y: proportion of x's that are greater + 0.5 * ties
  # This gives Pr(X > Y | Y = yj)
  W <- sapply(y, function(yj) {
    mean(x > yj) + 0.5 * mean(x == yj)
  })

  # p_hat = Pr(X > Y) - this matches the sign convention of rbs_calc
  # where positive rb means x tends to be larger than y
  p_hat <- mean(V)

  list(V = V, W = W, p_hat = p_hat, n1 = n1, n2 = n2)
}

#' Compute variance of concordance probability (two-sample) using Lehmann/Agresti formula
#' @param placements output from compute_placements
#' @return variance of p_hat
#' @references
#' Agresti, A. (1980). Generalized odds ratios for ordinal data. Biometrics, 36, 59-67.
#' Lehmann, E.L. (1975). Nonparametrics: Statistical Methods Based on Ranks. Holden-Day.
#' @noRd
var_concordance_twosample <- function(placements) {
  with(placements, {
    # Variance components from Lehmann (1975, eq. 2.21)
    # With our convention where p_hat = Pr(X > Y):
    # V[i] = Pr(X > Y | X = x_i), so E[V^2] estimates Pr(X1 > Y, X2 > Y)
    # W[j] = Pr(X > Y | Y = y_j), so E[W^2] estimates Pr(X > Y1, X > Y2)
    P_v <- mean(V^2)  # Pr(X1 > Y, X2 > Y) for two X's vs one Y
    P_w <- mean(W^2)  # Pr(X > Y1, X > Y2) for one X vs two Y's

    # Variance of p_hat = Pr(X > Y)
    # Var(p_hat) = (1/n1) * (P_v - p^2) + (1/n2) * (P_w - p^2)
    var_p <- (P_v - p_hat^2) / n1 + (P_w - p_hat^2) / n2

    max(var_p, 0)  # Ensure non-negative
  })
}

#' Compute variance of concordance probability (paired/one-sample)
#' @param d vector of differences (y - x - mu)
#' @return variance of p_hat for paired design
#' @noRd
var_concordance_paired <- function(d) {
  # Remove zeros
  d_nonzero <- d[d != 0]
  n <- length(d_nonzero)

  if (n < 2) return(NA)

  # Tie correction
  tie_counts <- table(abs(d_nonzero))
  tie_correction <- sum(tie_counts * (tie_counts - 1) * (tie_counts + 1) / 2)

  # Variance of W (Wilcoxon statistic)
  var_W <- (n * (n + 1) * (2 * n + 1) - tie_correction) / 24

  # Total rank sum
  S <- n * (n + 1) / 2

  # Variance of p_hat = W / S
  var_p <- var_W / S^2

  max(var_p, 0)
}

#' Compute SEs for all non-parametric effect sizes using Agresti method
#' @param x numeric vector (group 1 or pre-treatment)
#' @param y numeric vector (group 2 or post-treatment), NULL for one-sample
#' @param paired logical, TRUE for paired samples
#' @param mu hypothesized difference (default 0)
#' @param use_score_fallback logical, if TRUE and two-sample boundary is hit, use score CI
#' @param conf.level confidence level for score CI fallback (default 0.95)
#' @return list with point estimates and SEs for all effect sizes, including boundary_corrected flag
#' @noRd
ses_compute_agresti <- function(x, y = NULL, paired = FALSE, mu = 0,
                                 use_score_fallback = TRUE, conf.level = 0.95) {

  if (is.null(y)) {
    # One-sample: compare x to mu (treat as paired with y = mu)
    y <- rep(mu, length(x))
    paired <- TRUE
    mu <- 0  # Already accounted for
  }

  # Track if continuity correction was applied
  boundary_corrected <- FALSE
  boundary_used_score_ci <- FALSE
  score_ci_cstat <- NULL

  # Store original p_hat before any correction (needed for score fallback)
  p_hat_original <- NULL

  if (paired) {
    # === PAIRED SAMPLES ===
    # Use rbs_calc to get the point estimate (ensures consistency)
    # Natural order: P(X - Y > 0)
    r_hat <- rbs_calc(x = x, y = y, mu = mu, paired = TRUE)
    p_hat <- rb_to_cstat(r_hat)
    p_hat_original <- p_hat

    # Compute variance using Agresti method
    d <- x - y - mu
    d_nonzero <- d[d != 0]
    n_nonzero <- length(d_nonzero)

    if (n_nonzero < 2) {
      warning("Fewer than 2 non-zero differences")
      return(NULL)
    }

    # Haldane-type boundary correction for paired/one-sample
    # N_pairs = N*(N+1)/2 is the maximum possible Wilcoxon signed-rank statistic
    N_pairs <- n_nonzero * (n_nonzero + 1) / 2

    if (p_hat <= 0 || p_hat >= 1) {
      boundary_corrected <- TRUE
      # Concordance count: C = p_hat * N_pairs (will be 0 or N_pairs at boundary)
      C <- p_hat * N_pairs
      # Haldane-type shrinkage correction
      p_hat <- (C + 0.5) / (N_pairs + 1)
      r_hat <- 2 * p_hat - 1
    }

    # Variance using Agresti method (based on non-zero differences)
    # At boundary, recompute variance with corrected p_hat
    var_p <- var_concordance_paired(d)

    # If variance is degenerate at boundary, use a fallback based on corrected p_hat
    if (boundary_corrected && (is.na(var_p) || var_p <= 0)) {
      # Use approximate variance based on (1 - p^2) / n structure
      var_p <- (1 - p_hat^2) * (2 * n_nonzero + 1) / (6 * N_pairs^2)
    }

    se_p <- sqrt(max(var_p, 0))
    n1 <- n_nonzero
    n2 <- NULL

  } else {
    # === TWO INDEPENDENT SAMPLES ===
    x <- na.omit(x)
    y <- na.omit(y)

    n1 <- length(x)
    n2 <- length(y)

    # Compute placements and variance
    placements <- compute_placements(x - mu, y)
    p_hat <- placements$p_hat
    p_hat_original <- p_hat
    N_pairs <- n1 * n2

    if (p_hat <= 0 || p_hat >= 1) {
      boundary_corrected <- TRUE
      # Concordance count: C = p_hat * N_pairs (will be 0 or N_pairs at boundary)
      C <- p_hat * N_pairs
      # Haldane-type shrinkage correction
      p_hat <- (C + 0.5) / (N_pairs + 1)

      # For two-sample, optionally fall back to score CI for better boundary behavior
      if (use_score_fallback) {
        tf <- wmw_tie_factor(x - mu, y)
        score_ci_cstat <- score_ci_wmw(phi_hat = p_hat_original, tf = tf,
                                        n1 = n1, n2 = n2,
                                        conf.level = conf.level, correct = FALSE)
        boundary_used_score_ci <- TRUE
      }
    }

    r_hat <- 2 * p_hat - 1

    # Recompute variance with corrected p_hat
    # At boundary, placement values are all identical, so need to use corrected p_hat
    if (boundary_corrected) {
      # At complete separation, the placement-based variance formula fails
      # (all V_i and W_j are 0 or 1, leading to negative variance).
      # Use the V_LAPH variance formula evaluated at the corrected p_hat.
      tf <- wmw_tie_factor(x - mu, y)
      var_p <- v_laph(p_hat, tf, n1, n2)
    } else {
      var_p <- var_concordance_twosample(placements)
    }

    se_p <- sqrt(max(var_p, 0))
  }

  # Safeguard for extreme values
  p_hat_safe <- pmin(pmax(p_hat, 1e-10), 1 - 1e-10)

  # Point estimates for all scales
  alpha_hat <- p_hat / (1 - p_hat)
  eta_hat <- log(alpha_hat)  # = qlogis(p_hat)

  # Standard errors (delta method)
  se_cstat <- se_p
  se_rb <- 2 * se_p
  se_odds <- se_p / (1 - p_hat_safe)^2
  se_logodds <- se_p / (p_hat_safe * (1 - p_hat_safe))

  list(
    # Point estimates
    cstat = p_hat,
    rb = r_hat,
    odds = alpha_hat,
    logodds = eta_hat,
    # Standard errors
    se_cstat = se_cstat,
    se_rb = se_rb,
    se_odds = se_odds,
    se_logodds = se_logodds,
    # Design info
    paired = paired,
    n1 = n1,
    n2 = n2,
    # Boundary correction flags
    boundary_corrected = boundary_corrected,
    boundary_used_score_ci = boundary_used_score_ci,
    score_ci_cstat = score_ci_cstat,
    # Original uncorrected p_hat (for reference)
    p_hat_original = p_hat_original
  )
}

#' Compute CIs for non-parametric effect sizes via log-odds transformation
#' @param ses_results output from ses_compute_agresti
#' @param conf.level confidence level (default 0.95)
#' @return list with CIs for all effect sizes
#' @noRd
ses_ci_logodds <- function(ses_results, conf.level = 0.95) {

  alpha <- 1 - conf.level
  z <- qnorm(1 - alpha / 2)

  with(ses_results, {
    # Safeguard for extreme values
    p_safe <- pmin(pmax(cstat, 1e-10), 1 - 1e-10)
    eta_safe <- log(p_safe / (1 - p_safe))
    se_eta_safe <- se_cstat / (p_safe * (1 - p_safe))

    # CI on log-odds scale (best asymptotic properties)
    ci_logodds <- eta_safe + c(-1, 1) * z * se_eta_safe

    # Back-transform to all other scales
    ci_odds <- exp(ci_logodds)
    ci_cstat <- plogis(ci_logodds)
    ci_rb <- 2 * ci_cstat - 1

    list(
      ci_cstat = ci_cstat,
      ci_rb = ci_rb,
      ci_odds = ci_odds,
      ci_logodds = ci_logodds,
      conf.level = conf.level
    )
  })
}

#' Compute SE for rank-biserial using Fisher z-transform (legacy method)
#' @param n1 sample size for group 1 (or total n for paired)
#' @param n2 sample size for group 2 (NULL for paired)
#' @param paired logical
#' @return SE on Fisher z scale
#' @noRd
se_fisher_z <- function(n1, n2 = NULL, paired = FALSE) {
  if (paired || is.null(n2)) {
    # Paired/one-sample
    nd <- n1
    maxw <- (nd^2 + nd) / 2
    rfSE <- sqrt((2 * nd^3 + 3 * nd^2 + nd) / 6) / maxw
  } else {
    # Two-sample
    rfSE <- sqrt((n1 + n2 + 1) / (3 * n1 * n2))
  }
  rfSE
}

#' Compute CIs using Fisher z-transform (legacy method)
#' @param r_rbs rank-biserial correlation estimate
#' @param n1 sample size for group 1
#' @param n2 sample size for group 2 (NULL for paired)
#' @param paired logical
#' @param conf.level confidence level
#' @return list with CIs and SE for all effect sizes
#' @noRd
ses_ci_fisher <- function(r_rbs, n1, n2 = NULL, paired = FALSE, conf.level = 0.95) {
  alpha <- 1 - conf.level

  # Fisher z-transform
  rf <- atanh(r_rbs)
  rfSE <- se_fisher_z(n1, n2, paired)

  # CI on Fisher z scale, back-transform to rb
  ci_rb <- tanh(rf + c(-1, 1) * qnorm(1 - alpha / 2) * rfSE)

  # Transform to other scales
  ci_cstat <- rb_to_cstat(ci_rb)
  ci_odds <- rb_to_odds(ci_rb)
  ci_logodds <- log(ci_odds)

  # SE for rb (on original scale, using delta method from Fisher z)
  # d(tanh(z))/dz = sech^2(z) = 1 - tanh^2(z)
  se_rb_fisher <- rfSE * (1 - r_rbs^2)

  # SE for other scales (delta method from rb)
  p_hat <- rb_to_cstat(r_rbs)
  p_hat_safe <- pmin(pmax(p_hat, 1e-10), 1 - 1e-10)
  se_cstat_fisher <- se_rb_fisher / 2
  se_odds_fisher <- se_cstat_fisher / (1 - p_hat_safe)^2
  se_logodds_fisher <- se_cstat_fisher / (p_hat_safe * (1 - p_hat_safe))

  list(
    ci_rb = ci_rb,
    ci_cstat = ci_cstat,
    ci_odds = ci_odds,
    ci_logodds = ci_logodds,
    se_rb = se_rb_fisher,
    se_cstat = se_cstat_fisher,
    se_odds = se_odds_fisher,
    se_logodds = se_logodds_fisher,
    conf.level = conf.level
  )
}
