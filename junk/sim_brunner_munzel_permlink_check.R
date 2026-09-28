# Simulation: perm-logit vs perm-probit intervals near the boundary --------
#
# Follow-up to junk/sim_brunner_munzel_logit_probit_permlogit.R. The
# permutation-calibrated logit interval undercovered in some small-sample cells
# with p = 0.8-0.95. This checks whether the probit link does better. Both
# intervals are computed from the SAME permuted datasets, so the comparison is
# paired. Only boundary-relevant cells are run (p = 0.8, 0.9, 0.95; n <= 15 per
# group, plus 10 v 20).
#
# Results are added to junk/bm_logit_probit_permlogit.md

# Settings --------

devtools::load_all(".", quiet = TRUE)

nsim <- as.integer(Sys.getenv("SIM_NSIM", 2000))
R_perm <- 999
alpha <- 0.05
n_cores <- max(1, parallel::detectCores() - 1)
out_file <- "junk/sim_brunner_munzel_permlink_check_results.rds"

# Helpers (same as sim_brunner_munzel_logit_probit_permlogit.R) --------

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

clamp01 <- function(p, eps) pmin(pmax(p, eps), 1 - eps)

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

gen_two <- function(nx, ny, p, sdx) {
  delta <- sqrt(sdx^2 + 1) * qnorm(p)
  list(x = rnorm(nx, delta, sdx), y = rnorm(ny, 0, 1))
}
gen_paired <- function(n, p, rho) {
  delta <- sqrt(2) * qnorm(p)
  z1 <- rnorm(n); z2 <- rho * z1 + sqrt(1 - rho^2) * rnorm(n)
  list(x = z1 + delta, y = z2)
}

# Studentized permutation interval on a link scale ("logit" or "probit")
perm_link_ci <- function(st_obs, perm_stats, alpha, p_method, link) {
  fwd <- if (link == "logit") qlogis else qnorm
  inv <- if (link == "logit") plogis else pnorm
  dlink <- if (link == "logit") function(p) 1 / (p * (1 - p)) else
    function(p) 1 / dnorm(qnorm(p))
  lstat <- function(st) {
    pl <- clamp01(st$pd, st$eps)
    c(L = fwd(pl), seL = st$se * dlink(pl))
  }
  obs <- lstat(st_obs)
  # permuted statistics centred at link(0.5) = 0 (exchangeability)
  Tperm <- vapply(perm_stats, function(st) {
    s <- lstat(st); s[["L"]] / s[["seL"]]
  }, numeric(1))
  crit <- TOSTER:::perm_crit(Tperm, alpha, p_method, "abs")
  inv(obs[["L"]] + c(-1, 1) * obs[["seL"]] * crit)
}

# One replication --------

one_rep <- function(cell, alpha, R_perm) {
  paired <- cell$paired
  d <- if (paired) gen_paired(cell$nx, cell$p, cell$rho) else
    gen_two(cell$nx, cell$ny, cell$p, cell$sdx)
  st <- if (paired) bmp_stats(d$x, d$y) else bm2_stats(d$x, d$y)
  ps <- if (paired) perm_stats_paired(d$x, d$y, R_perm) else
    perm_stats_two(d$x, d$y, R_perm)
  pm <- if (paired && cell$nx <= 13) "exact" else "plusone"

  cis <- list(perm_logit = perm_link_ci(st, ps, alpha, pm, "logit"),
              perm_probit = perm_link_ci(st, ps, alpha, pm, "probit"))
  cover <- vapply(cis, function(ci) ci[1] <= cell$p && cell$p <= ci[2], logical(1))
  width <- vapply(cis, function(ci) ci[2] - ci[1], numeric(1))
  c(setNames(cover, paste0("cov_", names(cis))),
    setNames(width, paste0("wid_", names(cis))),
    both_differ = unname(cover[1] != cover[2]))
}

# Cells --------

p_grid <- c(0.8, 0.9, 0.95)
two <- expand.grid(p = p_grid, design = c("7v7", "10v10", "15v15", "10v20"),
                   sdx = c(1, 3), stringsAsFactors = FALSE)
two$nx <- as.integer(sub("v.*", "", two$design))
two$ny <- as.integer(sub(".*v", "", two$design))
two$paired <- FALSE
two$rho <- NA

pr <- expand.grid(p = p_grid, nx = c(7L, 10L, 15L), rho = 0.5)
pr$ny <- pr$nx
pr$design <- paste0(pr$nx, " pairs")
pr$sdx <- 1
pr$paired <- TRUE

cells <- rbind(two[, c("paired", "design", "nx", "ny", "sdx", "rho", "p")],
               pr[, c("paired", "design", "nx", "ny", "sdx", "rho", "p")])

# Run --------

cl <- parallel::makeCluster(n_cores)
invisible(parallel::clusterEvalQ(cl, devtools::load_all(".", quiet = TRUE)))
parallel::clusterExport(cl, c("one_rep", "bm2_stats", "bmp_stats", "clamp01",
                              "perm_link_ci", "perm_stats_two",
                              "perm_stats_paired", "gen_two", "gen_paired"))
parallel::clusterSetRNGStream(cl, 20260929)

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

cat("\nnsim =", nsim, " R =", R_perm, " nominal coverage =", 1 - alpha,
    " (MC SE ~", round(sqrt(0.05 * 0.95 / nsim), 4), ")\n\n")
print(summary_tab[, c("paired", "design", "sdx", "p", "cov_perm_logit",
                      "cov_perm_probit", "wid_perm_logit", "wid_perm_probit",
                      "both_differ")], digits = 3, row.names = FALSE)

saveRDS(list(summary = summary_tab, nsim = nsim, R = R_perm, alpha = alpha,
             complete = TRUE, date = Sys.time()), out_file)
