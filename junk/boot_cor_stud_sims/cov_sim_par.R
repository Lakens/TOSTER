library(parallel)
# Smoke-test runs (SMOKE=1) write to a temp file so they never touch results
out_file <- if (Sys.getenv("SMOKE") == "1") tempfile(fileext = ".out") else
  "C:/GitHub/TOSTER/junk/boot_cor_stud_sims/cov_sim_par.out"
cat("", file = out_file)
work <- function(cell) {
  se_fun <- TOSTER:::.cor_se
  # HC-type leverage-corrected Pearson SE (delta: "hc3" or "hc4")
  se_hc <- function(x, y, type) {
    n <- length(x)
    psi <- TOSTER:::.cor_if_pearson(x, y)
    X <- cbind(x, y); md2 <- tryCatch(mahalanobis(X, colMeans(X), cov(X)), error = function(e) rep(NA_real_, n))
    h <- pmin(1 / n + md2 / (n - 1), 0.99)
    d <- if (type == "hc3") 2 else pmin(4, n * h / 3)
    sqrt(sum((psi * (1 - h)^(-d / 2))^2)) / n
  }
  one <- function(x, y, m, R = 999, a = .05) {
    n <- length(x); idx <- matrix(sample(n, n * R, TRUE), R)
    est <- cor(x, y, method = m); z0 <- atanh(est)
    b <- apply(idx, 1, function(i) cor(x[i], y[i], method = m)); zb <- atanh(b)
    q <- function(t) quantile(t[is.finite(t)], c(1 - a/2, a/2), names = FALSE)
    st <- function(f) { s0 <- f(x, y) / (1 - est^2)
      sb <- apply(idx, 1, function(i) f(x[i], y[i])) / (1 - b^2)
      tanh(z0 - q((zb - z0) / sb) * s0) }
    out <- list(basic_z = tanh(z0 - q(zb - z0)),
                stud_IF = st(function(x, y) se_fun(x, y, m)))
    if (m == "pearson") {
      out$stud_HC3 <- st(function(x, y) se_hc(x, y, "hc3"))
      out$stud_HC4 <- st(function(x, y) se_hc(x, y, "hc4"))
    }
    if (m == "spearman") {
      out$stud_ADFrank <- st(function(x, y) se_fun(rank(x), rank(y), "pearson"))
      f <- function(r) sqrt((1 + r^2/2) / (n - 3))
      out$old_stud <- tanh(z0 - q((zb - z0) / f(b)) * f(est))
    }
    jk <- vapply(1:n, function(j) cor(x[-j], y[-j], method = m), 0)
    out$bca <- TOSTER:::bca_ci(b, est, jk, a)
    out
  }
  gens <- list(
    normal = function(n) { x <- rnorm(n); list(x = x, y = .5*x + sqrt(.75)*rnorm(n)) },
    hetero = function(n) { x <- rnorm(n); list(x = x, y = .5*x + rnorm(n)*sqrt(1+x^2)) },
    t5 = function(n) { x <- rt(n, 5); list(x = x, y = .6*x + .8*rt(n, 5)) },
    t3 = function(n) { x <- rt(n, 3); list(x = x, y = .6*x + .8*rt(n, 3)) },
    discrete = function(n) { x <- sample(1:5, n, TRUE); list(x = x, y = pmin(5, x + sample(0:2, n, TRUE))) })
  g <- cell$g; m <- cell$m; n <- cell$n
  set.seed(cell$seed)
  analytic <- c(normal = .5, hetero = 1/3, t5 = .6, t3 = .6)
  truth <- if (m == "pearson" && g %in% names(analytic)) analytic[[g]] else
    if (m == "kendall") mean(replicate(20, { d <- gens[[g]](3000); cor(d$x, d$y, method = m) })) else
    { d <- gens[[g]](3e5); cor(d$x, d$y, method = m) }
  hits <- replicate(as.integer(Sys.getenv("NSIM", "1000")), { d <- gens[[g]](n); ci <- one(d$x, d$y, m)
    sapply(ci, function(c) c[1] < truth & truth < c[2]) })
  line <- sprintf("%-8s %-8s n=%2d  %s\n", g, m, n,
                  paste(sprintf("%s=%.3f", rownames(hits), rowMeans(hits)), collapse = "  "))
  cat(line, file = out_file, append = TRUE)
  line
}
cells <- expand.grid(n = c(80, 30), m = c("kendall", "spearman", "pearson"),
                     g = c("normal", "hetero", "t5", "t3", "discrete"), stringsAsFactors = FALSE)
cells <- cells[order(cells$m != "kendall", -cells$n), ]
cells$seed <- seq_len(nrow(cells)) + 2026
# Smoke test: SMOKE=1 NSIM=3 runs the t3 Pearson cells serially
if (Sys.getenv("SMOKE") == "1") { devtools::load_all("C:/GitHub/TOSTER", quiet = TRUE); for (i in which(cells$m == "pearson" & cells$g == "t3")) print(work(cells[i, ])); quit() }
cl <- makeCluster(12)
clusterEvalQ(cl, devtools::load_all("C:/GitHub/TOSTER", quiet = TRUE))
clusterExport(cl, "out_file")
res <- parLapplyLB(cl, split(cells, seq_len(nrow(cells))), work)
stopCluster(cl)
cat("DONE\n", file = out_file, append = TRUE)
