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
