# Simulation: range-preserving intervals for the Brunner-Munzel relative effect --------
#
# Question: is a probit option, or a permutation-calibrated logit interval, a
# useful addition to brunner_munzel()? Compares 95% CI coverage of the true
# relative effect p = P(X > Y) + 0.5 * P(X = Y) for:
#   - t          : brunner_munzel(test_method = "t")       (package)
#   - logit      : brunner_munzel(test_method = "logit")   (package)
#   - probit     : delta-method probit interval            (prototype below)
#   - perm       : brunner_munzel(test_method = "perm")    (package; clamped to [0, 1])
#   - perm_logit : studentized permutation test on the logit scale, inverted to
#                  a CI and back-transformed (prototype below; range-preserving)
#
# Coverage at p = 0.5 equals 1 - Type I error of the two-sided test of H0: p = 0.5.
# Scenarios with p = 0.8, 0.9, 0.95 probe behaviour near the boundary, which is
# where range-preserving intervals are supposed to help.
#
# Results are summarized in junk/bm_logit_probit_permlogit.md

# Settings --------

devtools::load_all(".", quiet = TRUE)

nsim <- as.integer(Sys.getenv("SIM_NSIM", 2000))
R_perm <- 999
alpha <- 0.05
n_cores <- max(1, parallel::detectCores() - 1)
out_file <- "junk/sim_brunner_munzel_logit_probit_permlogit_results.rds"

p_grid <- c(0.5, 0.8, 0.9, 0.95)

# Prototype helpers (mirror brunner_munzel() internals) --------

# Two-sample: estimate, SE, and Satterthwaite df exactly as in brunner_munzel()
bm2_stats <- function(x, y) {
  nx <- length(x); ny <- length(y); N <- nx + ny
  rxy <- rank(c(x, y)); rx <- rank(x); ry <- rank(y)
  pl2 <- (rxy[1:nx] - rx) / ny
  pl1 <- (rxy[(nx + 1):N] - ry) / nx
  pd <- mean(pl2)
  s1 <- var(pl2) / nx
  s2 <- var(pl1) / ny
  V <- N * (s1 + s2)
  if (V == 0) V <- N * 0.5 / (nx * ny)^2
  df <- (s1 + s2)^2 / (s1^2 / (nx - 1) + s2^2 / (ny - 1))
  if (is.nan(df)) df <- 1000
  list(pd = pd, se = sqrt(V / N), df = df, eps = 0.5 / (nx * ny))
}

# Paired: estimate and SE exactly as in brunner_munzel(paired = TRUE)
bmp_stats <- function(x, y) {
  n <- length(x)
  all_data <- c(y, x)
  rx <- rank(all_data)
  BM1 <- (rx[1:n] - rank(y)) / n
  BM2 <- (rx[(n + 1):(2 * n)] - rank(x)) / n
  BM3 <- BM1 - BM2
  pd <- mean(BM2)
  v <- (sum(BM3^2) - n * mean(BM3)^2) / (n - 1)
  if (v == 0) v <- 1 / n
  list(pd = pd, se = sqrt(v / n), df = n - 1, eps = 0.5 / n^2)
}

# Move an estimate of exactly 0 or 1 half a step inside (continuity clamp)
clamp01 <- function(p, eps) pmin(pmax(p, eps), 1 - eps)

# Delta-method interval on a transformed scale (link = "logit" or "probit")
link_ci <- function(st, link, alpha) {
  pl <- clamp01(st$pd, st$eps)
  if (link == "logit") {
    L <- qlogis(pl); seL <- st$se / (pl * (1 - pl)); inv <- plogis
  } else {
    L <- qnorm(pl); seL <- st$se / dnorm(L); inv <- pnorm
  }
  q <- qt(1 - alpha / 2, st$df)
  inv(L + c(-1, 1) * q * seL)
}

# Studentized permutation interval on the logit scale
# perms: list of (x, y) permuted datasets; stats_fun: bm2_stats or bmp_stats
perm_logit_ci <- function(st_obs, perm_stats, alpha, p_method) {
  lstat <- function(st) {
    pl <- clamp01(st$pd, st$eps)
    c(L = qlogis(pl), seL = st$se / (pl * (1 - pl)))
  }
  obs <- lstat(st_obs)
  # permuted statistics centred at logit(0.5) = 0 (exchangeability)
  Tperm <- vapply(perm_stats, function(st) {
    s <- lstat(st); s[["L"]] / s[["seL"]]
  }, numeric(1))
  crit <- TOSTER:::perm_crit(Tperm, alpha, p_method, "abs")
  plogis(obs[["L"]] + c(-1, 1) * obs[["seL"]] * crit)
}

# Permutation sets matching brunner_munzel(): two-sample samples R
# reassignments; paired enumerates 2^n swaps for n <= 13, else samples R
perm_stats_two <- function(x, y, R) {
  z <- c(x, y); nx <- length(x); N <- length(z)
  lapply(seq_len(R), function(i) {
    idx <- sample.int(N, nx)
    bm2_stats(z[idx], z[-idx])
  })
}
perm_stats_paired <- function(x, y, R) {
  n <- length(x)
  S <- if (n <= 13) {
    as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), n)))
  } else {
    matrix(sample(c(FALSE, TRUE), n * R, replace = TRUE), nrow = R)
  }
  lapply(seq_len(nrow(S)), function(i) {
    s <- S[i, ]
    bmp_stats(ifelse(s, y, x), ifelse(s, x, y))
  })
}

# Data generation with a known relative effect --------

# Two-sample normal: X ~ N(delta, sdx), Y ~ N(0, 1), p = pnorm(delta / sqrt(sdx^2 + 1))
gen_two <- function(nx, ny, p, sdx) {
  delta <- sqrt(sdx^2 + 1) * qnorm(p)
  list(x = rnorm(nx, delta, sdx), y = rnorm(ny, 0, 1))
}
# Paired bivariate normal, correlation rho, marginals N(delta, 1) and N(0, 1)
gen_paired <- function(n, p, rho) {
  delta <- sqrt(2) * qnorm(p)
  z1 <- rnorm(n); z2 <- rho * z1 + sqrt(1 - rho^2) * rnorm(n)
  list(x = z1 + delta, y = z2)
}

# One replication --------

one_rep <- function(cell, alpha, R_perm) {
  paired <- cell$paired
  d <- if (paired) gen_paired(cell$nx, cell$p, cell$rho) else
    gen_two(cell$nx, cell$ny, cell$p, cell$sdx)

  pkg <- function(method) {
    tryCatch(suppressWarnings(suppressMessages(
      brunner_munzel(d$x, d$y, paired = paired, test_method = method,
                     R = R_perm, alpha = alpha)))$conf.int,
      error = function(e) c(NA_real_, NA_real_))
  }
  st <- if (paired) bmp_stats(d$x, d$y) else bm2_stats(d$x, d$y)
  n_exact <- paired && cell$nx <= 13
  ps <- if (paired) perm_stats_paired(d$x, d$y, R_perm) else
    perm_stats_two(d$x, d$y, R_perm)

  cis <- list(
    t = pkg("t"),
    logit = pkg("logit"),
    probit = link_ci(st, "probit", alpha),
    perm = pkg("perm"),
    perm_logit = perm_logit_ci(st, ps, alpha, if (n_exact) "exact" else "plusone")
  )
  cover <- vapply(cis, function(ci) ci[1] <= cell$p && cell$p <= ci[2], logical(1))
  width <- vapply(cis, function(ci) ci[2] - ci[1], numeric(1))
  at_bound <- vapply(cis, function(ci) ci[1] <= 0 || ci[2] >= 1, logical(1))
  c(setNames(cover, paste0("cov_", names(cis))),
    setNames(width, paste0("wid_", names(cis))),
    setNames(at_bound, paste0("bnd_", names(cis))))
}

# Cells --------

two <- expand.grid(p = p_grid, design = c("7v7", "10v10", "15v15", "30v30", "10v20"),
                   sdx = c(1, 3), stringsAsFactors = FALSE)
two$nx <- as.integer(sub("v.*", "", two$design))
two$ny <- as.integer(sub(".*v", "", two$design))
two$paired <- FALSE
two$rho <- NA

pr <- expand.grid(p = p_grid, nx = c(7L, 10L, 15L, 30L), rho = 0.5)
pr$ny <- pr$nx
pr$design <- paste0(pr$nx, " pairs")
pr$sdx <- 1
pr$paired <- TRUE

cells <- rbind(two[, c("paired", "design", "nx", "ny", "sdx", "rho", "p")],
               pr[, c("paired", "design", "nx", "ny", "sdx", "rho", "p")])

# Sanity check: prototypes reproduce the package's logit interval --------

set.seed(1)
chk <- replicate(20, {
  d <- gen_two(8, 12, 0.8, 1)
  a <- link_ci(bm2_stats(d$x, d$y), "logit", alpha)
  b <- suppressWarnings(suppressMessages(
    brunner_munzel(d$x, d$y, test_method = "logit")))$conf.int
  dp <- gen_paired(9, 0.8, 0.5)
  ap <- link_ci(bmp_stats(dp$x, dp$y), "logit", alpha)
  bp <- suppressWarnings(suppressMessages(
    brunner_munzel(dp$x, dp$y, paired = TRUE, test_method = "logit")))$conf.int
  max(abs(a - b), abs(ap - bp))
})
cat("Max |prototype logit CI - package logit CI| over 20 datasets:",
    signif(max(chk), 3), "\n")
cat("(differences arise only when the estimate is exactly 0 or 1, where the",
    "prototype uses a half-step continuity clamp)\n")

# Run --------

cl <- parallel::makeCluster(n_cores)
invisible(parallel::clusterEvalQ(cl, devtools::load_all(".", quiet = TRUE)))
parallel::clusterExport(cl, c("one_rep", "bm2_stats", "bmp_stats", "clamp01",
                              "link_ci", "perm_logit_ci", "perm_stats_two",
                              "perm_stats_paired", "gen_two", "gen_paired"))
parallel::clusterSetRNGStream(cl, 20260928)

t0 <- Sys.time()
results <- list()
for (k in seq_len(nrow(cells))) {
  cell <- as.list(cells[k, ])
  reps <- parallel::parSapply(cl, seq_len(nsim), function(i, cell, alpha, R_perm) {
    one_rep(cell, alpha, R_perm)
  }, cell = cell, alpha = alpha, R_perm = R_perm)
  results[[k]] <- cbind(cells[k, ], t(rowMeans(reps, na.rm = TRUE)))
  message(k, "/", nrow(cells), " ", cell$design, " sdx=", cell$sdx, " p=", cell$p,
          " done (", format(round(Sys.time() - t0, 1)), ")")
  saveRDS(list(summary = do.call(rbind, results), nsim = nsim, R = R_perm,
               alpha = alpha, complete = FALSE), out_file)
}
parallel::stopCluster(cl)

summary_tab <- do.call(rbind, results)
rownames(summary_tab) <- NULL

# Summary --------

cat("\nnsim =", nsim, " R =", R_perm, " nominal coverage =", 1 - alpha,
    " (MC SE ~", round(sqrt(0.05 * 0.95 / nsim), 4), ")\n\n")
print(summary_tab[, c("paired", "design", "sdx", "p",
                      grep("^cov_", names(summary_tab), value = TRUE))],
      digits = 3, row.names = FALSE)
cat("\nMean width:\n")
print(summary_tab[, c("paired", "design", "sdx", "p",
                      grep("^wid_", names(summary_tab), value = TRUE))],
      digits = 3, row.names = FALSE)
cat("\nProportion of intervals touching 0 or 1:\n")
print(summary_tab[, c("paired", "design", "sdx", "p",
                      grep("^bnd_", names(summary_tab), value = TRUE))],
      digits = 3, row.names = FALSE)

saveRDS(list(summary = summary_tab, nsim = nsim, R = R_perm, alpha = alpha,
             complete = TRUE, date = Sys.time()), out_file)
