rbs_calc = function (x, y,
                     mu, paired) {
  # adapted from effectsize R package
  if (paired) {
    z <- (x-y) - mu
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
#' @return list with point estimates and SEs for all effect sizes, including boundary_corrected flag
#' @noRd
ses_compute_agresti <- function(x, y = NULL, paired = FALSE, mu = 0) {

  if (is.null(y)) {
    # One-sample: compare x to mu (treat as paired with y = mu)
    y <- rep(mu, length(x))
    paired <- TRUE
    mu <- 0  # Already accounted for
  }

  # Track if continuity correction was applied

  boundary_corrected <- FALSE

  if (paired) {
    # === PAIRED SAMPLES ===
    # Use rbs_calc to get the point estimate (ensures consistency)
    # Note: rbs_calc expects x and y swapped for paired (it computes x - y)
    r_hat <- rbs_calc(x = y, y = x, mu = mu, paired = TRUE)
    p_hat <- rb_to_cstat(r_hat)

    # Compute variance using Agresti method
    d <- y - x - mu
    d_nonzero <- d[d != 0]
    n_nonzero <- length(d_nonzero)

    if (n_nonzero < 2) {
      warning("Fewer than 2 non-zero differences")
      return(NULL)
    }

    # Continuity correction for boundary cases (paired/one-sample)
    # S = total rank sum, analogous to n1 * n2 in the two-sample case
    S_total <- n_nonzero * (n_nonzero + 1) / 2
    if (p_hat <= 0) {
      p_hat <- 0.5 / S_total
      boundary_corrected <- TRUE
    } else if (p_hat >= 1) {
      p_hat <- 1 - 0.5 / S_total
      boundary_corrected <- TRUE
    }

    # Recompute r_hat from corrected p_hat if needed
    if (boundary_corrected) {
      r_hat <- 2 * p_hat - 1
    }

    # Variance using Agresti method (based on non-zero differences)
    var_p <- var_concordance_paired(d)
    se_p <- sqrt(var_p)

  } else {
    # === TWO INDEPENDENT SAMPLES ===
    x <- na.omit(x)
    y <- na.omit(y)

    # Compute placements and variance
    placements <- compute_placements(x - mu, y)
    p_hat <- placements$p_hat
    n_pairs <- placements$n1 * placements$n2

    # Continuity correction for boundary cases (two-sample)
    if (p_hat <= 0) {
      p_hat <- 0.5 / n_pairs
      boundary_corrected <- TRUE
    } else if (p_hat >= 1) {
      p_hat <- 1 - 0.5 / n_pairs
      boundary_corrected <- TRUE
    }

    r_hat <- 2 * p_hat - 1

    var_p <- var_concordance_twosample(placements)
    se_p <- sqrt(var_p)
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
    # Boundary correction flag
    boundary_corrected = boundary_corrected
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
