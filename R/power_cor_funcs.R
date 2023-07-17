pow_corr = function (n = NULL, r = NULL, power = NULL, null = 0,
          alpha = NULL, alternative = c("two.sided", "less", "greater"))
{

  if (sum(sapply(list(n, r, power, alpha), is.null)) != 1)
    stop("exactly one of n, r, power, and alpha must be NULL")
  if (!is.null(alpha) && !is.numeric(alpha) || any(0 > alpha |
                                                   alpha > 1))
    stop(sQuote("alpha"), " must be numeric in [0, 1]")
  if (!is.null(power) && !is.numeric(power) || any(0 > power |
                                                   power > 1))
    stop(sQuote("power"), " must be numeric in [0, 1]")
  if (!is.null(n) && min(n) < 4)
    stop("number of observations must be at least 4")
  p=0
  alternative <- match.arg(alternative)
  tside <- switch(alternative, less = 1, two.sided = 2, greater = 3)
  if (tside == 2 && !is.null(r))
    r <- abs(r)
  if (tside == 3) {
    p.body <- quote({
      delta <- sqrt(n - 3 - p) * (log((1 + r)/(1 - r))/2 +
                                    r/(n - 1 - p)/2 * (1 + (5 + r^2)/(n - 1 - p)/4 +
                                                         (11 + 2 * r^2 + 3 * r^4)/(n - 1 - p)^2/8) -
                                    log((1 + null)/(1 - null))/2 - null/(n - 1 -
                                                                           p)/2)
      v <- (n - 3 - p)/(n - 1 - p) * (1 + (4 - r^2)/(n -
                                                       1 - p)/2 + (22 - 6 * r^2 - 3 * r^4)/(n - 1 -
                                                                                              p)^2/6)
      zalpha <- qnorm(1 - alpha)
      pnorm((delta - zalpha)/sqrt(v))
    })
  }
  if (tside == 1) {
    p.body <- quote({
      delta <- sqrt(n - 3 - p) * (log((1 + r)/(1 - r))/2 +
                                    r/(n - 1 - p)/2 * (1 + (5 + r^2)/(n - 1 - p)/4 +
                                                         (11 + 2 * r^2 + 3 * r^4)/(n - 1 - p)^2/8) -
                                    log((1 + null)/(1 - null))/2 - null/(n - 1 -
                                                                           p)/2)
      v <- (n - 3 - p)/(n - 1 - p) * (1 + (4 - r^2)/(n -
                                                       1 - p)/2 + (22 - 6 * r^2 - 3 * r^4)/(n - 1 -
                                                                                              p)^2/6)
      zalpha <- qnorm(1 - alpha)
      pnorm((-delta - zalpha)/sqrt(v))
    })
  }
  if (tside == 2) {
    p.body <- quote({
      delta <- sqrt(n - 3 - p) * (log((1 + r)/(1 - r))/2 +
                                    r/(n - 1 - p)/2 * (1 + (5 + r^2)/(n - 1 - p)/4 +
                                                         (11 + 2 * r^2 + 3 * r^4)/(n - 1 - p)^2/8) -
                                    log((1 + null)/(1 - null))/2 - null/(n - 1 -
                                                                           p)/2)
      v <- (n - 3 - p)/(n - 1 - p) * (1 + (4 - r^2)/(n -
                                                       1 - p)/2 + (22 - 6 * r^2 - 3 * r^4)/(n - 1 -
                                                                                              p)^2/6)
      zalpha <- qnorm(1 - alpha/2)
      pnorm((delta - zalpha)/sqrt(v)) + pnorm((-delta -
                                                 zalpha)/sqrt(v))
    })
  }
  if (is.null(power))
    power <- eval(p.body)
  else if (is.null(n))
    n <- uniroot(function(n) eval(p.body) - power, c(4 +
                                                       p + 1e-10, 1e+07))$root
  else if (is.null(r)) {
    if (tside == 2) {
      r <- uniroot(function(r) eval(p.body) - power, c(1e-10,
                                                       1 - 1e-10))$root
    }
    else {
      r <- uniroot(function(r) eval(p.body) - power, c(-1 +
                                                         1e-10, 1 - 1e-10))$root
    }
  }
  else if (is.null(alpha))
    alpha <- uniroot(function(alpha) eval(p.body) - power,
                     c(1e-10, 1 - 1e-10))$root
  else stop("internal error")
  METHOD <- "Power for Pearson Product-Moment Correlation"

  structure(list(n = n, rho = r,
                 alpha = alpha, beta = 1-power, power = power,
                 null = null, alternative = alternative,
                 method = METHOD),
            class = "power.htest")
}

pow_corr_tost = function (n = NULL, r = 0, power = NULL, null = NULL,
                       alpha = NULL)
{

  if (sum(sapply(list(n, r, power, alpha), is.null)) != 1)
    stop("exactly one of n, r, power, and alpha must be NULL")
  if(is.null(r)){
    stop("r cannot be set to NULL at this time.")
  }
  if (!is.null(alpha) && !is.numeric(alpha) || any(0 > alpha |
                                                   alpha > 1))
    stop(sQuote("alpha"), " must be numeric in [0, 1]")
  if (!is.null(power) && !is.numeric(power) || any(0 > power |
                                                   power > 1))
    stop(sQuote("power"), " must be numeric in [0, 1]")
  if (!is.null(n) && min(n) < 4)
    stop("number of observations must be at least 4")

  alternative <- "equivalence"

  p.body <- quote({
    se <- 1/sqrt(n-3)
    C1 = 0.5 * log((1+min(null))/(1-min(null)))
    C2 = 0.5 * log((1+max(null))/(1-max(null)))
    statistical_power1 <- 2 * (pnorm((abs(C1)/se) - qnorm(1-alpha)) + pnorm(-(abs(C1)/se) - qnorm(1-alpha))) - 1
    statistical_power2 <- 2 * (pnorm((abs(C2)/se) - qnorm(1-alpha)) + pnorm(-(abs(C2)/se) - qnorm(1-alpha))) - 1
    min(statistical_power1, statistical_power2)
  })

  if (is.null(power))
    power <- eval(p.body)
  else if (is.null(n))
    n <- uniroot(function(n) eval(p.body) - power, c(4 + 1e-10, 1e+07))$root
  else if (is.null(alpha))
    alpha <- uniroot(function(alpha) eval(p.body) - power,
                     c(1e-10, 1 - 1e-10))$root
  else stop("internal error")
  METHOD <- "Power for Pearson Product-Moment Correlation"

  structure(list(n = n, rho = r,
                 alpha = alpha, beta = 1-power, power = power,
                 null = null, alternative = alternative,
                 method = METHOD),
            class = "power.htest")
}


