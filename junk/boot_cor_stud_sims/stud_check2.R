# Vectorized: X, Y are R x n matrices of resampled data
mom <- function(X, Y) {
  n <- ncol(X)
  cx <- X - rowMeans(X); cy <- Y - rowMeans(Y)
  zx <- cx / sqrt(rowMeans(cx^2)); zy <- cy / sqrt(rowMeans(cy^2))
  r <- rowMeans(zx * zy)
  v <- (r^2/4) * (rowMeans(zx^4) + rowMeans(zy^4) + 2*rowMeans(zx^2*zy^2)) -
    r * (rowMeans(zx^3*zy) + rowMeans(zx*zy^3)) + rowMeans(zx^2*zy^2)
  list(r = r, se_r = sqrt(pmax(v, 1e-12) / n))
}
jack_se_z <- function(X, Y) { # leave-one-out via sums, vectorized over rows
  n <- ncol(X)
  Sx <- rowSums(X); Sy <- rowSums(Y); Sxx <- rowSums(X^2); Syy <- rowSums(Y^2); Sxy <- rowSums(X*Y)
  zj <- matrix(0, nrow(X), n)
  for (j in 1:n) {
    m <- n - 1
    sx <- Sx - X[, j]; sy <- Sy - Y[, j]
    cxy <- (Sxy - X[, j]*Y[, j]) - sx*sy/m
    vx <- (Sxx - X[, j]^2) - sx^2/m; vy <- (Syy - Y[, j]^2) - sy^2/m
    zj[, j] <- atanh(pmin(pmax(cxy/sqrt(vx*vy), -0.9999999), 0.9999999))
  }
  sqrt((n-1)/n * rowSums((zj - rowMeans(zj))^2))
}
one <- function(x, y, R = 999, alpha = .05) {
  n <- length(x); idx <- matrix(sample(n, n*R, TRUE), R)
  X <- matrix(x[idx], R); Y <- matrix(y[idx], R)
  mb <- mom(X, Y); m0 <- mom(matrix(x, 1), matrix(y, 1))
  r0 <- m0$r; z0 <- atanh(r0); zb <- atanh(mb$r)
  out <- list()
  qf <- function(t) quantile(t, c(1 - alpha/2, alpha/2), na.rm = TRUE, names = FALSE)
  # old: constant SE -> basic on z
  out$old <- tanh(z0 - qf(zb - z0))
  # ADF on z scale
  sez_b <- mb$se_r / (1 - mb$r^2); sez_0 <- m0$se_r / (1 - r0^2)
  out$adf_z <- tanh(z0 - qf((zb - z0)/sez_b) * sez_0)
  # ADF on r scale
  out$adf_r <- pmin(pmax(r0 - qf((mb$r - r0)/mb$se_r) * m0$se_r, -1), 1)
  # jackknife SE on z
  sj_b <- jack_se_z(X, Y); sj_0 <- jack_se_z(matrix(x, 1), matrix(y, 1))
  out$jack_z <- tanh(z0 - qf((zb - z0)/sj_b) * sj_0)
  out$perc <- quantile(mb$r, c(alpha/2, 1 - alpha/2), names = FALSE)
  out
}
gens <- list(
  normal = function(n, rho) { x <- rnorm(n); y <- rho*x + sqrt(1-rho^2)*rnorm(n); list(x=x,y=y) },
  hetero = function(n, rho) { x <- rnorm(n); y <- rho*x + rnorm(n)*sqrt(1+x^2); list(x=x,y=y) },
  t5     = function(n, rho) { x <- rt(n,5); y <- rho*x + sqrt(1-rho^2)*rt(n,5); list(x=x,y=y) }
)
set.seed(2026)
nsim <- 1000
for (g in names(gens)) {
  big <- gens[[g]](2e6, .5); rt <- cor(big$x, big$y)
  for (n in c(30, 80, 200)) {
    hits <- replicate(nsim, { d <- gens[[g]](n, .5); ci <- one(d$x, d$y)
      sapply(ci, function(c) c[1] < rt & rt < c[2]) })
    cat(sprintf("%-7s n=%3d rho=%.3f  ", g, n, rt)); print(round(rowMeans(hits), 3))
  }
}
