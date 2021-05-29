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


  if (type == "two.sample"){
    df = ifelse(length(n) > 1 , sum(n)-2, 2*n-2)
    #nc <- sum(1/n)
    nc = ifelse(length(n) > 1 , sum(1/n), 2*(1/(n)))

  } else {
    df = n-1
    nc <- sum(1/n)
  }
  se.fac <- sqrt(1 * nc) #sqrt(ades$bkni * nc) ** error for power.TOST


  if (any(df < 1))
    stop("n too small. Degrees of freedom <1!")
  pow <- calc_power_theta(alpha,
                          ltheta1 = theta1,
                          ltheta2 = theta2,
                          diffm = theta0,
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

OwensQOwen =  function (nu, t, delta, a = 0, b)
{
  if (nu < 1)
    stop("nu must be >=1!")
  if (a != 0)
    stop("Only a=0 implemented!")
  if (!is.finite(b))
    return(pt(t, df = nu, ncp = delta))
  if (!is.finite(delta))
    delta <- sign(delta) * 1e+20
  A <- t/sqrt(nu)
  B <- nu/(nu + t * t)
  upr <- nu - 2
  av <- vector(mode = "numeric", length = nu)
  for (k in seq_along(av)) {
    if (k == 1 | k == 2)
      av[k] <- 1
    else av[k] <- 1/((k - 2) * av[k - 1])
  }
  ll <- ifelse((upr - 1) > 0, upr - 1, 0)
  L <- vector(mode = "numeric", length = ll)
  if (is.finite(b)) {
    for (k in seq_along(L)) {
      if (k == 1)
        L[1] <- 0.5 * A * B * b * dnorm(b) * dnorm(A *
                                                     b - delta)
      else L[k] <- av[k + 3] * b * L[k - 1]
    }
  }
  ll <- ifelse((upr + 1) > 0, upr + 1, 0)
  H <- vector(mode = "numeric", length = ll)
  if (is.finite(b)) {
    for (k in seq_along(H)) {
      if (k == 1)
        H[1] <- -dnorm(b) * pnorm(A * b - delta)
      else H[k] <- av[k + 1] * b * H[k - 1]
    }
  }
  M <- vector(mode = "numeric", length = ll)
  sB <- sqrt(B)
  for (k in seq_along(M)) {
    if (k == 1)
      M[1] <- A * sB * dnorm(delta * sB) * (pnorm(delta *
                                                    A * sB) - pnorm((delta * A * B - b)/sB))
    if (k == 2)
      M[2] <- B * (delta * A * M[1] + A * dnorm(delta *
                                                  sB) * (dnorm(delta * A * sB) - dnorm((delta *
                                                                                          A * B - b)/sB)))
    if (k > 2)
      M[k] <- ((k - 2)/(k - 1)) * B * (av[k - 1] * delta *
                                         A * M[k - 1] + M[k - 2]) - L[k - 2]
  }
  sumt <- 0
  if (2 * (nu%/%2) != nu) {
    if (upr >= 1) {
      k <- seq(1, upr, by = 2)
      sumt <- sum(M[k + 1]) + sum(H[k + 1])
    }
    qv <- pnorm(b) - 2 * OwensT(b, A - delta/b) - 2 * OwensT(delta *
                                                               sB, (delta * A * B - b)/B/delta) + 2 * OwensT(delta *
                                                                                                               sB, A) - (delta >= 0) + 2 * sumt
  }
  else {
    if (upr >= 0) {
      k <- seq(0, upr, by = 2)
      sumt <- sum(M[k + 1]) + sum(H[k + 1])
    }
    qv <- pnorm(-delta) + sqrt(2 * pi) * sumt
  }
  return(qv)
}

OwensT = function (h, a)
{
  eps <- .Machine$double.eps
  if (abs(a) < eps | !is.finite(h) | abs(1 - abs(a)) < eps |
      abs(h) < eps | !is.finite(abs(a))) {
    if (abs(a) < eps)
      return(0)
    if (!is.finite(h))
      return(0)
    if (abs(1 - abs(a)) < eps)
      return(sign(a) * 0.5 * pnorm(h) * (1 - pnorm(h)))
    if (abs(h) < eps)
      return(atan(a)/2/pi)
    if (!is.finite(abs(a))) {
      if (h < 0)
        tha <- pnorm(h)/2
      else tha <- (1 - pnorm(h))/2
      return(sign(a) * tha)
    }
  }
  aa <- abs(a)
  if (aa <= 1) {
    tha <- tfn(h, a)
    return(tha)
  }
  else {
    ah <- aa * h
    gh <- pnorm(h)
    gah <- pnorm(ah)
    tha <- 0.5 * (gh + gah) - gh * gah - tfn(ah, 1/aa)
  }
  if (a < 0)
    tha <- -tha
  return(tha)
}

tfn = function (x, fx)
{
  ng <- 5
  r <- c(0.1477621, 0.1346334, 0.1095432, 0.0747257, 0.0333357)
  u <- c(0.0744372, 0.2166977, 0.3397048, 0.4325317, 0.4869533)
  tp <- 1/(2 * pi)
  tv1 <- .Machine$double.eps
  tv2 <- 15
  tv3 <- 15
  tv4 <- 1e-05
  if (tv2 < abs(x))
    return(0)
  xs <- -0.5 * x * x
  x2 <- fx
  fxs <- fx * fx
  if (tv3 <= log(1 + fxs) - xs * fxs) {
    x1 <- 0.5 * fx
    fxs <- 0.25 * fxs
    while (1) {
      rt <- fxs + 1
      x2 <- x1 + (xs * fxs + tv3 - log(rt))/(2 * x1 *
                                               (1/rt - xs))
      fxs <- x2 * x2
      if (abs(x2 - x1) < tv4)
        break
      x1 <- x2
    }
  }
  r1 <- 1 + fxs * (0.5 + u)^2
  r2 <- 1 + fxs * (0.5 - u)^2
  rt <- sum(r * (exp(xs * r1)/r1 + exp(xs * r2)/r2))
  return(rt * x2 * tp)
}


