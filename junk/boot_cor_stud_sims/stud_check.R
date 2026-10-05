# ADF (distribution-free) SE of Pearson r, on Fisher z scale
adf_se_z <- function(x, y) {
  n <- length(x)
  zx <- (x - mean(x)) / sqrt(mean((x - mean(x))^2))
  zy <- (y - mean(y)) / sqrt(mean((y - mean(y))^2))
  r <- mean(zx * zy)
  m40 <- mean(zx^4); m04 <- mean(zy^4); m22 <- mean(zx^2 * zy^2)
  m31 <- mean(zx^3 * zy); m13 <- mean(zx * zy^3)
  v <- (r^2 / 4) * (m40 + m04 + 2 * m22) - r * (m31 + m13) + m22
  sqrt(v / n) / (1 - r^2)
}
ci_one <- function(x, y, R = 999, alpha = .05, type) {
  n <- length(x); r <- cor(x, y); z0 <- atanh(r)
  idx <- matrix(sample(n, n * R, TRUE), R)
  zs <- se <- numeric(R)
  for (b in 1:R) { i <- idx[b, ]; zs[b] <- atanh(cor(x[i], y[i]))
    se[b] <- if (type == "old") 1 / sqrt(n - 3) else adf_se_z(x[i], y[i]) }
  se0 <- if (type == "old") 1 / sqrt(n - 3) else adf_se_z(x, y)
  q <- quantile((zs - z0) / se, c(1 - alpha / 2, alpha / 2), na.rm = TRUE)
  tanh(z0 - q * se0)
}
set.seed(1)
gen <- function(n, rho) { # heteroscedastic, non-normal: classic failure case
  x <- rnorm(n); e <- rnorm(n) * sqrt(1 + x^2)
  y <- rho * x + e; list(x = x, y = y) }
# true rho for this DGP
big <- gen(2e6, .5); rho_true <- cor(big$x, big$y)
cov <- sapply(c("old", "adf"), function(tp) mean(replicate(400, {
  d <- gen(40, .5); ci <- ci_one(d$x, d$y, type = tp); ci[1] < rho_true & rho_true < ci[2] })))
print(round(rho_true, 3)); print(cov)
