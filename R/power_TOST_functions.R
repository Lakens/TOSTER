pow_tTOST = function (alpha = 0.05,
                      theta1,
                      theta2,
                      theta0,
                      sd = 1,
                      n,
                      type = "two.sample")
{
  if (missing(sd))
    stop("sd must be given!")
  if (missing(n))
    stop("Number of subjects (n) must be given!")
  if (length(alpha) != 1)
    stop("alpha must be a scalar!")

  if (length(n) != 1) {
    if(length(n) > 2){
      stop("Only 2 sample sizes can be provided for a two sample case")
    }

    if(type != "two.sample"){
      stop("length of sample sizes can only be 1 for paired or one sample cases")
    }
  }

  nc <- sum(1/n)
  if (type == "two.sample"){
    df = ifelse(length(n) > 1 , sum(n)-2, 2*n-2)
  } else {
    df = n-1
  }
  se.fac <- ifelse(type == "two.sample" , sqrt(1 * nc), sqrt(2 * nc)) sqrt(ades$bkni * nc)


  if (any(df < 1))
    stop("n too small. Degrees of freedom <1!")
  pow <- calc_power_theta(alpha,
                          ltheta1 = theta1,
                          ltheta2 = theta2,
                          ldiff = theta0,
                          sem = sd * se.fac,
                          df = df)
  return(pow)
}


calc_power_theta = function (alpha = 0.05,
                             ltheta1,
                             ltheta2,
                             diffm,
                             sem,
                             df)
{
  stopifnot(length(alpha) <= 2)
  tval <- qt(1 - alpha, df, lower.tail = TRUE)
  dl <- length(tval)
  delta1 <- (diffm - ltheta1)/sem
  delta2 <- (diffm - ltheta2)/sem
  delta1[is.nan(delta1)] <- 0
  delta2[is.nan(delta2)] <- 0
  R <- (delta1 - delta2) * sqrt(df)/(tval[1] + tval[dl])
  R[is.nan(R)] <- 0
  R[R <= 0] <- Inf
  if (min(df) > 10000) {
    tval <- qnorm(1 - alpha)
    p1 <- pnorm(tval[1] - delta1)
    p2 <- pnorm(-tval[dl] - delta2)
    pwr <- p2 - p1
    pwr[pwr < 0] <- 0
    return(pwr)
  }
  if (min(df) >= 5000 & min(df <= 10000)) {
    return(approx.power.TOST(alpha, ltheta1, ltheta2, diffm,
                              sem, df))
  }
  p1 <- vector("numeric", length(delta1))
  p2 <- vector("numeric", length(delta1))
  for (i in seq_along(delta1)) {
    p1[i] <- OwensQ(df[1], tval[1], delta1[i], 0, R[i])
    p2[i] <- OwensQ(df[1], -tval[dl], delta2[i], 0, R[i])
  }
  pwr <- p2 - p1
  pwr[pwr < 0] <- 0
  return(pwr)
}

approx.power.TOST = function (alpha = 0.05, ltheta1, ltheta2, diffm, sem, df)
{
  stopifnot(length(alpha) <= 2)
  tval <- qt(1 - alpha, df, lower.tail = TRUE, log.p = FALSE)
  delta1 <- (diffm - ltheta1)/sem
  delta2 <- (diffm - ltheta2)/sem
  delta1[is.nan(delta1)] <- 0
  delta2[is.nan(delta2)] <- 0
  pow <- suppressWarnings(pt(-tval[length(tval)], df, ncp = delta2) -
                            pt(tval[1], df, ncp = delta1))
  pow[pow < 0] <- 0
  return(pow)
}

OwensQ = function (nu, t, delta, a = 0, b)
{
  if (missing(b))
    stop("Upper integration limit missing.")
  if (a != 0)
    stop("Only a==0 implemented.")
  if (length(nu) > 1 | length(t) > 1 | length(delta) > 1 |
      length(a) > 1 | length(b) > 1)
    stop("Input must be scalars!")
  if (nu < 1)
    stop("nu must be >=1.")
  if (a == b)
    return(0)
  if (nu < 29 && abs(delta) > 37.62) {
    if (is.infinite(b)) {
      i_fun <- function(y) .Q.integrand(y/(1 - y), nu,
                                        t, delta) * 1/(1 - y)^2
      return(integrate(i_fun, lower = 0, upper = 1, subdivisions = 1000L,
                       rel.tol = 1e-08, stop.on.error = TRUE)[[1]])
    }
    else {
      return(OwensQOwen(nu, t, delta, 0, b))
    }
  }
  else {
    if (is.infinite(b)) {
      return(suppressWarnings(pt(t, df = nu, ncp = delta)))
    }
    else {
      i_fun <- function(y) {
        .Q.integrand(b + y/(1 - y), nu, t, delta)/(1 -
                                                     y)^2
      }
      Integral01 <- integrate(i_fun, lower = 0, upper = 1,
                              subdivisions = 1000L, rel.tol = 1e-08, stop.on.error = TRUE)
      return(suppressWarnings(pt(t, df = nu, ncp = delta)) -
               Integral01[[1]])
    }
  }
}

.Q.integrand = function (x, nu, t, delta)
{
  lnQconst <- -((nu/2) - 1) * log(2) - lgamma(nu/2)
  dens <- x
  x1 <- x[x != 0]
  dens[x != 0] <- sign(x1)^(nu - 1) * pnorm(t * x1/sqrt(nu) -
                                              delta, mean = 0, sd = 1, log.p = FALSE) * exp((nu -
                                                                                               1) * log(abs(x1)) - 0.5 * x1^2 + lnQconst)
  dens
}



