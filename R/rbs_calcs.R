rbs_calc = function (x, y,
                     mu, paired) {
  # adapted from effectsize R package
  if (paired) {
    Ry <- ranktransform((x - y) - mu, sign = TRUE)
    Ry <- na.omit(Ry)
    n <- length(Ry)
    S <- (n * (n + 1)/2)
    U1 <- sum(Ry[Ry > 0], na.rm = TRUE)
    U2 <- -sum(Ry[Ry < 0], na.rm = TRUE)
  }
  else {
    Ry <- ranktransform(c(x - mu, y))
    n1 <- length(x)
    n2 <- length(y)
    S <- (n1 * n2)
    U1 <- sum(Ry[seq_along(x)]) - n1 * (n1 + 1)/2
    U2 <- sum(Ry[-seq_along(x)]) - n2 * (n2 + 1)/2
  }
  u_ <- U1/S
  f_ <- U2/S
  return(u_ - f_)
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

