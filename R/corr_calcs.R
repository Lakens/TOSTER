cor_to_ci <- function(cor, n, ci = 0.95,
                      method = "pearson",
                      correction = "fieller",
                      se = NULL, ...) {
  method <- match.arg(tolower(method),
                      c("pearson", "kendall", "spearman"),
                      several.ok = FALSE)

  if (method == "kendall") {
    out <- .cor_to_ci_kendall(cor, n,
                              ci = ci,
                              correction = correction,
                              se = se, ...)
  } else if (method == "spearman") {
    out <- .cor_to_ci_spearman(cor, n,
                               ci = ci,
                               correction = correction,
                               se = se, ...)
  } else {
    out <- .cor_to_ci_pearson(cor, n, ci = ci, se = se, ...)
  }

  out
}

# Kendall -----------------------------------------------------------------
#' @importFrom stats qnorm
.cor_to_ci_kendall <- function(cor, n,
                               ci = 0.95,
                               correction = "fieller",
                               se = NULL, ...) {
  # by @tsbaguley (https://rpubs.com/seriousstats/616206)

  if (!is.null(se)) {
    tau.se <- se
  } else if (correction == "fieller") {
    tau.se <- (0.437 / (n - 4))^0.5
  } else {
    tau.se <- 1 / (n - 3)^0.5
  }

  moe <- stats::qnorm(1 - (1 - ci) / 2) * tau.se
  zu <- atanh(cor) + moe
  zl <- atanh(cor) - moe

  # Convert back to r
  ci_low <- tanh(zl)
  ci_high <- tanh(zu)

  c(ci_low, ci_high)
}


# Spearman -----------------------------------------------------------------
.cor_to_ci_spearman <- function(cor, n,
                                ci = 0.95,
                                correction = "fieller",
                                se = NULL, ...) {
  # by @tsbaguley (https://rpubs.com/seriousstats/616206)

  if (!is.null(se)) {
    zrs.se <- se
  } else if (correction == "fieller") {
    zrs.se <- (1.06 / (n - 3))^0.5
  } else if (correction == "bw") {
    zrs.se <- ((1 + (cor^2) / 2) / (n - 3))^0.5
  } else {
    zrs.se <- 1 / (n - 3)^0.5
  }

  moe <- stats::qnorm(1 - (1 - ci) / 2) * zrs.se

  zu <- atanh(cor) + moe
  zl <- atanh(cor) - moe

  # Convert back to r
  ci_low <- tanh(zl)
  ci_high <- tanh(zu)

  c(ci_low, ci_high)
}


# Pearson -----------------------------------------------------------------
.cor_to_ci_pearson <- function(cor, n, ci = 0.95, se = NULL, ...) {
  z <- atanh(cor)
  if (is.null(se)) {
    se <- 1 / sqrt(n - 3) # Sample standard error
  }

  # CI
  alpha <- 1 - (1 - ci) / 2
  ci_low <- z - se * stats::qnorm(alpha)
  ci_high <- z + se * stats::qnorm(alpha)

  # Convert back to r
  ci_low <- tanh(ci_low)
  ci_high <- tanh(ci_high)

  c(ci_low, ci_high)
}
# TOSTER:::cor_to_ci(.5,18,"pearson",ci=.95)

corr_curv = function (r,
                      n,
                      type = "pearson",
                      steps = 5000) {
  intrvls <- (0:steps)/steps
  intrvls = subset(intrvls,intrvls>0 & intrvls<1)


  results <-
    suppressWarnings({
      lapply(
        intrvls,
        FUN = function(i)
          cor_to_ci(
            cor = r, n = n,
            method = type,
            correction = "fieller",
            ci = i
          )
      )
    })

  df <- data.frame(do.call(rbind, results))
  intrvl.limit <- c("lower.limit", "upper.limit")
  colnames(df) <- intrvl.limit
  df$intrvl.width <- (abs((df$upper.limit) - (df$lower.limit)))
  df$intrvl.level <- intrvls
  df$cdf <- (abs(df$intrvl.level/2)) + 0.5
  df$pvalue <- 1 - intrvls
  df$svalue <- -log2(df$pvalue)
  df <- head(df, -1)
  class(df) <- c("data.frame", "concurve")
  densdf <- data.frame(c(df$lower.limit, df$upper.limit))
  colnames(densdf) <- "x"
  densdf <- head(densdf, -1)
  class(densdf) <- c("data.frame", "concurve")

  return(list(df, densdf))

}

# plot_cor(.5,18)

.corboot <- function(isub, x, y, method){
  res <- cor(x[isub], y[isub], method = method)
  res
}

# Other corrs -----


# winsorized

wincor <- function(x, y, tr=0.2){

  df <- cbind(x,y)
  df <- df[complete.cases(df), ]
  n <- nrow(df)
  x <- df[,1]
  y <- df[,2]
  g <- floor(tr*n)
  xvec <- winval(x,tr)
  yvec <- winval(y,tr)
  corv <- cor(xvec,yvec)

  return(corv)
}

winval <- function(x, tr=0.2){
  y <- sort(x)
  n <- length(x)
  ibot <- floor(tr*n)+1
  itop <- length(x)-ibot+1
  xbot <- y[ibot]
  xtop <- y[itop]
  winval <- ifelse(x<=xbot,xbot,x)
  winval <- ifelse(winval>=xtop,xtop,winval)
  winval
}

# percentage bend measure of location
pbos <- function(x, beta=.2){
  n <- length(x)
  temp <- sort(abs(x-median(x)))
  omhatx <- temp[floor((1-beta)*n)]
  psi <- (x-median(x))/omhatx
  i1 <- length(psi[psi<(-1)])
  i2 <- length(psi[psi>1])
  sx <- ifelse(psi<(-1),0,x)
  sx <- ifelse(psi>1,0,sx)
  pbos <- (sum(sx) + omhatx*(i2-i1)) / (n-i1-i2)
  pbos
}

pbcor <- function(x, y, beta=.2){

  df <- cbind(x,y)
  df <- df[complete.cases(df), ]
  n <- nrow(df)
  x <- df[,1]
  y <- df[,2]
  temp <- sort(abs(x-median(x)))
  omhatx <- temp[floor((1-beta)*n)]
  temp <- sort(abs(y-median(y)))
  omhaty <- temp[floor((1-beta)*n)]
  a <- (x-pbos(x,beta))/omhatx
  b <- (y-pbos(y,beta))/omhaty
  a <- ifelse(a<=-1,-1,a)
  a <- ifelse(a>=1,1,a)
  b <- ifelse(b<=-1,-1,b)
  b <- ifelse(b>=1,1,b)
  corv <- sum(a*b)/sqrt(sum(a^2)*sum(b^2))

  return(corv)
}

.corboot_wincor <- function(isub, x, y, ...){
  res <- wincor(x[isub], y[isub], ...)
  res
}

.corboot_pbcor <- function(isub, x, y, ...){
  res <- pbcor(x[isub], y[isub], ...)
  res
}

# Studentized bootstrap CI for correlations -----

#' Studentized bootstrap CI on the correlation scale
#'
#' Computes a studentized (bootstrap-t) confidence interval by pivoting on the
#' Fisher z scale and then back-transforming.
#'
#' @param tvec Numeric vector of bootstrap pivots: (z_star - z_obs) / se_star
#' @param t0_z Observed Fisher z value: atanh(est)
#' @param se_obs Analytical SE of z_obs
#' @param alpha Two-tailed significance level (e.g., 0.05 for 95% CI)
#' @return Numeric vector of length 2: c(lower, upper) on the correlation scale
#' @keywords internal
stud_ci <- function(tvec, t0_z, se_obs, alpha) {
  qs <- quantile(tvec, probs = c(1 - alpha / 2, alpha / 2), names = FALSE)
  z_bounds <- t0_z - qs * se_obs
  tanh(z_bounds)
}

# Bootstrap p-value dispatch -----

#' Compute a bootstrap p-value consistent with the selected CI method
#'
#' Dispatches to the appropriate p-value calculation based on the CI method,
#' ensuring that p < alpha if and only if the corresponding CI excludes the null.
#'
#' @param bvec Numeric vector of bootstrap correlation estimates
#' @param est Observed correlation estimate
#' @param null Null hypothesis value (single numeric)
#' @param alternative One of "two.sided", "greater", "less"
#' @param boot_ci One of "perc", "basic", "bca", "stud"
#' @param tvec Bootstrap pivots (required for "stud")
#' @param se_obs Analytical SE on z scale (required for "stud")
#' @param z0 BCa bias correction (required for "bca")
#' @param acc BCa acceleration (required for "bca")
#' @param nboot Number of bootstrap replicates
#' @return A single p-value
#' @keywords internal
boot_pvalue <- function(bvec, est, null, alternative,
                        boot_ci, tvec = NULL, se_obs = NULL,
                        z0 = NULL, acc = NULL, nboot) {

  if (boot_ci == "perc") {
    sig <- .pval_perc(bvec, null, alternative, nboot)
  } else if (boot_ci == "basic") {
    sig <- .pval_basic(bvec, est, null, alternative, nboot)
  } else if (boot_ci == "bca") {
    sig <- .pval_bca(bvec, est, null, alternative, nboot,
                     z0 = z0, acc = acc)
  } else if (boot_ci == "stud") {
    sig <- .pval_stud(tvec, est, null, alternative, se_obs, nboot)
  }

  sig
}

# Percentile p-value (Wilcox method) -----
.pval_perc <- function(bvec, null, alternative, nboot) {
  if (alternative == "two.sided") {
    phat <- (sum(bvec < null) + 0.5 * sum(bvec == null)) / nboot
    sig <- 2 * min(phat, 1 - phat)
  } else if (alternative == "greater") {
    sig <- 1 - sum(bvec >= null) / nboot
  } else { # less
    sig <- 1 - sum(bvec <= null) / nboot
  }
  sig
}

# Basic (reflected) p-value -----
.pval_basic <- function(bvec, est, null, alternative, nboot) {
  reflected <- 2 * est - bvec
  if (alternative == "two.sided") {
    phat <- (sum(reflected < null) + 0.5 * sum(reflected == null)) / nboot
    sig <- 2 * min(phat, 1 - phat)
  } else if (alternative == "greater") {
    sig <- (sum(reflected < null) + 0.5 * sum(reflected == null)) / nboot
  } else { # less
    sig <- (sum(reflected > null) + 0.5 * sum(reflected == null)) / nboot
  }
  sig
}

# BCa inverted p-value -----
.pval_bca <- function(bvec, est, null, alternative, nboot,
                      z0, acc) {
  # Continuity-corrected proportion below null
  p0 <- (sum(bvec < null) + 0.5) / (nboot + 1)

  z_p0 <- qnorm(p0)

  # BCa-adjusted cumulative probability at the null
  u <- z_p0 - z0
  p_bca <- pnorm((u * (1 - acc * z0) - z0) / (1 + acc * u))

  if (alternative == "two.sided") {
    sig <- 2 * min(p_bca, 1 - p_bca)
  } else if (alternative == "greater") {
    sig <- p_bca
  } else { # less
    sig <- 1 - p_bca
  }
  sig
}

# Studentized (pivot) p-value -----
.pval_stud <- function(tvec, est, null, alternative, se_obs, nboot) {
  t_obs <- (atanh(est) - atanh(null)) / se_obs

  if (alternative == "two.sided") {
    sig <- 2 * min(mean(tvec >= t_obs), mean(tvec <= t_obs))
  } else if (alternative == "greater") {
    sig <- mean(tvec <= t_obs)
  } else { # less
    sig <- mean(tvec >= t_obs)
  }
  sig
}

# SE on Fisher z scale for studentized bootstrap -----

#' Analytical SE on the Fisher z scale for a given correlation method
#'
#' @param r_star Correlation estimate (can be a vector for bootstrap replicates)
#' @param n Sample size
#' @param method One of "pearson", "kendall", "spearman"
#' @return SE on the Fisher z scale
#' @keywords internal
.fisher_z_se <- function(r_star, n, method) {
  if (method == "pearson") {
    1 / sqrt(n - 3)
  } else if (method == "spearman") {
    sqrt((1 + r_star^2 / 2) / (n - 3))
  } else if (method == "kendall") {
    sqrt(0.437 / (n - 4))
  }
}
