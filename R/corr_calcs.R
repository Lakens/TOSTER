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

#' Studentized bootstrap CI for a correlation
#'
#' Computes a studentized (bootstrap-t) confidence interval on the working
#' scale (Fisher z or correlation) and back-transforms it to the correlation
#' scale.
#'
#' @param tvec Numeric vector of bootstrap pivots on the working scale:
#'   (t_star - t0) / se_star
#' @param t0 Observed estimate on the working scale
#' @param se_obs Observed SE on the working scale
#' @param alpha Two-tailed significance level (e.g., 0.05 for 95% CI)
#' @param back Function mapping the working scale back to the correlation
#'   scale (`tanh` for Fisher z, `identity` for the correlation scale)
#' @return Numeric vector of length 2: c(lower, upper) on the correlation
#'   scale, truncated to \[-1, 1\]
#' @keywords internal
stud_ci <- function(tvec, t0, se_obs, alpha, back = identity) {
  qs <- quantile(tvec, probs = c(1 - alpha / 2, alpha / 2),
                 names = FALSE, na.rm = TRUE)
  pmin(pmax(back(t0 - qs * se_obs), -1), 1)
}

# Bootstrap p-value dispatch -----

#' Compute a bootstrap p-value consistent with the selected CI method
#'
#' Dispatches to the appropriate p-value calculation based on the CI method,
#' ensuring that p < alpha if and only if the corresponding CI excludes the null.
#' Used by all bootstrap functions in the package (correlations, t-tests, SMD,
#' SES, TOST, and log-TOST).
#'
#' @param bvec Numeric vector of bootstrap estimates (on the working scale)
#' @param est Observed estimate (on the same scale as bvec)
#' @param null Null hypothesis value (single numeric, on the same scale as bvec)
#' @param alternative One of "two.sided", "greater", "less"
#' @param boot_ci One of "perc", "basic", "bca", "stud"
#' @param tvec Bootstrap pivots (required for "stud"):
#'   typically `(bvec - est) / bootstrap_se`
#' @param se_obs Observed standard error on the working scale (required for "stud")
#' @param z0 BCa bias correction from `bca_params()` (required for "bca")
#' @param acc BCa acceleration from `bca_params()` (required for "bca")
#' @param nboot Number of bootstrap replicates
#' @return A single p-value
#' @keywords internal
boot_pvalue <- function(bvec, est, null, alternative,
                        boot_ci, tvec = NULL, se_obs = NULL,
                        z0 = NULL, acc = NULL, nboot) {

  boot_ci <- match.arg(boot_ci, c("perc", "basic", "bca", "stud"))

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
# Note on tie handling: the two-sided case uses the standard 0.5 continuity
# correction for ties at the null (Wilcox, 2022). For one-sided tests, the
# strict inequality convention (>= / <=) is used without tie correction,
# following the standard percentile bootstrap p-value definition. This is a
# deliberate methodological choice; see .pval_basic() for an alternative
# that applies tie correction in the one-sided cases.
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
  t_obs <- (est - null) / se_obs
  tvec <- tvec[!is.na(tvec)]

  if (alternative == "two.sided") {
    sig <- 2 * min(mean(tvec >= t_obs), mean(tvec <= t_obs))
  } else if (alternative == "greater") {
    sig <- mean(tvec >= t_obs)
  } else { # less
    sig <- mean(tvec <= t_obs)
  }
  sig
}

# Influence-function SEs for studentized bootstrap -----

#' Influence-function standard error of a correlation coefficient
#'
#' Estimates the standard error of a Pearson, Spearman, or Kendall correlation
#' from the data, without assuming bivariate normality. The SE is
#' \eqn{\sqrt{\sum_i \psi_i^2} / n}, where \eqn{\psi_i} is the empirical
#' influence function of the coefficient evaluated at observation \eqn{i}
#' (for Pearson, inflated by an HC4-type leverage correction).
#' Because the SE is re-estimated from each bootstrap resample, it can be used
#' to studentize the bootstrap. See the "Studentized bootstrap" section of
#' [boot_cor_test()] for the formulas.
#'
#' @param x,y Numeric vectors of equal length (no missing values)
#' @param method One of "pearson", "spearman", "kendall"
#' @return SE on the correlation scale
#' @keywords internal
.cor_se <- function(x, y, method) {
  psi <- switch(method,
                pearson = .cor_if_pearson(x, y) * .hc4_weight(x, y),
                spearman = .cor_if_spearman(x, y),
                kendall = .cor_if_kendall(x, y))
  sqrt(sum(psi^2)) / length(x)
}

# HC4-type leverage inflation, (1 - h_i)^(-delta_i / 2), with leverage
# h_i = 1/n + MD_i^2 / (n - 1) of each point in (x, y), so sum(h) = p = 3
.hc4_weight <- function(x, y) {
  n <- length(x)
  X <- cbind(x, y)
  md2 <- tryCatch(stats::mahalanobis(X, colMeans(X), stats::cov(X)),
                  error = function(e) rep(NA_real_, n))
  # cap leverage so a lone extreme point in a tiny resample stays finite
  h <- pmin(1 / n + md2 / (n - 1), 0.99)
  delta <- pmin(4, n * h / 3)
  (1 - h)^(-delta / 2)
}

# Pearson: asymptotic distribution-free (fourth-moment) influence function
.cor_if_pearson <- function(x, y) {
  cx <- x - mean(x)
  cy <- y - mean(y)
  zx <- cx / sqrt(mean(cx^2))
  zy <- cy / sqrt(mean(cy^2))
  r <- mean(zx * zy)
  zx * zy - r * (zx^2 + zy^2) / 2
}

# Spearman: influence function of the Pearson correlation of the
# mid-distribution transforms U = F*(X), V = G*(Y) (i.e., Pearson on midranks).
# Each moment's influence has a direct term plus a term from estimating F*, G*.
# Without ties this reduces to the influence function of 12 E[F(X)G(Y)] - 3.
.cor_if_spearman <- function(x, y) {
  n <- length(x)
  u <- (rank(x) - 0.5) / n
  v <- (rank(y) - 0.5) / n
  mu <- mean(u)
  mv <- mean(v)
  su <- mean(u^2) - mu^2
  sv <- mean(v^2) - mv^2
  rho <- (mean(u * v) - mu * mv) / sqrt(su * sv)
  # influence of each moment: direct + transform-estimation terms
  d_u <- (u - mu) + (.upper_mid_sum(x, rep(1, n)) / n - mu)
  d_v <- (v - mv) + (.upper_mid_sum(y, rep(1, n)) / n - mv)
  d_uv <- (u * v - mean(u * v)) +
    (.upper_mid_sum(x, v) / n - mean(u * v)) +
    (.upper_mid_sum(y, u) / n - mean(u * v))
  d_uu <- (u^2 - mean(u^2)) + 2 * (.upper_mid_sum(x, u) / n - mean(u^2))
  d_vv <- (v^2 - mean(v^2)) + 2 * (.upper_mid_sum(y, v) / n - mean(v^2))
  # delta method for rho = C / sqrt(Su * Sv)
  d_c <- d_uv - mv * d_u - mu * d_v
  d_su <- d_uu - 2 * mu * d_u
  d_sv <- d_vv - 2 * mv * d_v
  d_c / sqrt(su * sv) - rho / 2 * (d_su / su + d_sv / sv)
}

# For each i: sum_j w_j * (1{x_j > x_i} + 1{x_j == x_i} / 2), in O(n log n)
.upper_mid_sum <- function(x, w) {
  g <- match(x, sort(unique(x)))
  gs <- rowsum(w, g, reorder = TRUE)[, 1]
  above <- rev(cumsum(rev(gs))) - gs
  (above + gs / 2)[g]
}

# Kendall: Hoeffding projection of tau-b = A / sqrt(Bx * By), where A, Bx, By
# are U-statistics (concordance and non-tied pairs), via the delta method
.cor_if_kendall <- function(x, y) {
  n <- length(x)
  sx <- sign(outer(x, x, "-"))
  sy <- sign(outer(y, y, "-"))
  a <- rowSums(sx * sy) / (n - 1)
  bx <- rowSums(sx != 0) / (n - 1)
  by <- rowSums(sy != 0) / (n - 1)
  A <- mean(a)
  Bx <- mean(bx)
  By <- mean(by)
  tau <- A / sqrt(Bx * By)
  2 * ((a - A) / sqrt(Bx * By) -
         tau / 2 * ((bx - Bx) / Bx + (by - By) / By))
}
