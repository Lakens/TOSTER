# Simulation: Brunner-Munzel equivalence (TOST) and minimal effect testing --------
#
# Question: which test_method gives valid and efficient equivalence and minimal
# effect tests for the relative effect p = P(X > Y) + 0.5 * P(X = Y)?
#
# For each set of bounds, data are generated with the true p at the lower bound,
# at the upper bound, and at 0.5 (the centre). At a bound:
#   - equivalence (TOST) rejection rate = Type I error (should be <= alpha)
#   - minimal effect rejection rate     = Type I error (should be <= alpha)
# At p = 0.5 the equivalence rejection rate is power.
#
# Methods:
#   - t, logit : brunner_munzel(test_method = "t" / "logit")  (package)
#   - perm     : studentized permutation TOST on the probability scale
#                (prototype of the package's perm method; validated below
#                against brunner_munzel(test_method = "perm") in fully
#                enumerated cases)
#   - perm_logit : the same TOST on the logit scale (prototype)
# perm and perm_logit are computed from the SAME permutation set per dataset.
#
# Results are summarized in junk/bm_logit_probit_permlogit.md

# Settings --------

devtools::load_all(".", quiet = TRUE)

nsim <- as.integer(Sys.getenv("SIM_NSIM", 1500))
R_perm <- 999
alpha <- 0.05
n_cores <- max(1, parallel::detectCores() - 1)
out_file <- "junk/sim_brunner_munzel_tost_results.rds"

bounds_list <- list(`0.4-0.6` = c(0.4, 0.6), `0.3-0.7` = c(0.3, 0.7),
                    `0.2-0.8` = c(0.2, 0.8))

# Prototype helpers (mirror brunner_munzel() internals) --------

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
  list(pd = pd, se = sqrt(V / N), eps = 0.5 / (nx * ny))
}

bmp_stats <- function(x, y) {
  n <- length(x)
  rx <- rank(c(y, x))
  BM1 <- (rx[1:n] - rank(y)) / n
  BM2 <- (rx[(n + 1):(2 * n)] - rank(x)) / n
  BM3 <- BM1 - BM2
  pd <- mean(BM2)
  v <- (sum(BM3^2) - n * mean(BM3)^2) / (n - 1)
  if (v == 0) v <- 1 / n
  list(pd = pd, se = sqrt(v / n), eps = 0.5 / n^2)
}

clamp01 <- function(p, eps) pmin(pmax(p, eps), 1 - eps)

# Permutation sets matching brunner_munzel(): two-sample enumerates when
# choose(N, nx) <= R, else samples R reassignments; paired enumerates 2^n swaps
# for n <= 13, else samples R swap patterns. Returns list(stats, exact).
perm_set_two <- function(x, y, R) {
  z <- c(x, y); nx <- length(x); N <- length(z)
  if (choose(N, nx) <= R) {
    idx <- utils::combn(N, nx)
    st <- lapply(seq_len(ncol(idx)), function(i) bm2_stats(z[idx[, i]], z[-idx[, i]]))
    return(list(stats = st, exact = TRUE))
  }
  st <- lapply(seq_len(R), function(i) {
    idx <- sample.int(N, nx)
    bm2_stats(z[idx], z[-idx])
  })
  list(stats = st, exact = FALSE)
}
perm_set_paired <- function(x, y, R) {
  n <- length(x)
  exact <- n <= 13
  S <- if (exact) {
    as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), n)))
  } else {
    matrix(sample(c(FALSE, TRUE), n * R, replace = TRUE), nrow = R)
  }
  st <- lapply(seq_len(nrow(S)), function(i) {
    s <- S[i, ]
    bmp_stats(ifelse(s, y, x), ifelse(s, x, y))
  })
  list(stats = st, exact = exact)
}

# Permutation TOST and minimal effect p-values on a given scale
# scale = "identity" (package perm method) or "logit"
perm_tost <- function(st_obs, pset, bounds, scale) {
  f <- function(st) {
    if (scale == "identity") {
      c(L = st$pd, seL = st$se, centre = 0.5)
    } else {
      pl <- clamp01(st$pd, st$eps)
      c(L = qlogis(pl), seL = st$se / (pl * (1 - pl)), centre = 0)
    }
  }
  g <- if (scale == "identity") identity else qlogis
  obs <- f(st_obs)
  Tperm <- vapply(pset$stats, function(st) {
    s <- f(st); (s[["L"]] - s[["centre"]]) / s[["seL"]]
  }, numeric(1))
  t_low <- (obs[["L"]] - g(bounds[1])) / obs[["seL"]]
  t_high <- (obs[["L"]] - g(bounds[2])) / obs[["seL"]]
  R <- length(Tperm)
  pm <- if (pset$exact) "exact" else "plusone"
  pv <- function(b) TOSTER:::compute_perm_pval(b, R, pm)
  cnt <- function(t, type) TOSTER:::perm_count(Tperm, t, type)
  c(equ = max(pv(cnt(t_low, "ge")), pv(cnt(t_high, "le"))),
    met = min(pv(cnt(t_low, "le")), pv(cnt(t_high, "ge"))))
}

# Data generation with a known relative effect --------

gen_two <- function(nx, ny, p, sdx) {
  delta <- sqrt(sdx^2 + 1) * qnorm(p)
  list(x = rnorm(nx, delta, sdx), y = rnorm(ny, 0, 1))
}
gen_paired <- function(n, p, rho) {
  delta <- sqrt(2) * qnorm(p)
  z1 <- rnorm(n); z2 <- rho * z1 + sqrt(1 - rho^2) * rnorm(n)
  list(x = z1 + delta, y = z2)
}

# Validation: prototype perm TOST equals the package in enumerated cases --------

set.seed(1)
val <- replicate(10, {
  d <- gen_two(5, 5, 0.7, 1)
  proto <- perm_tost(bm2_stats(d$x, d$y), perm_set_two(d$x, d$y, 999),
                     c(0.3, 0.7), "identity")
  pk_e <- suppressWarnings(suppressMessages(brunner_munzel(
    d$x, d$y, test_method = "perm", alternative = "equivalence",
    mu = c(0.3, 0.7), R = 999)))$p.value
  pk_m <- suppressWarnings(suppressMessages(brunner_munzel(
    d$x, d$y, test_method = "perm", alternative = "minimal.effect",
    mu = c(0.3, 0.7), R = 999)))$p.value
  dp <- gen_paired(8, 0.7, 0.5)
  proto_p <- perm_tost(bmp_stats(dp$x, dp$y), perm_set_paired(dp$x, dp$y, 999),
                       c(0.3, 0.7), "identity")
  pkp_e <- suppressWarnings(suppressMessages(brunner_munzel(
    dp$x, dp$y, paired = TRUE, test_method = "perm",
    alternative = "equivalence", mu = c(0.3, 0.7))))$p.value
  pkp_m <- suppressWarnings(suppressMessages(brunner_munzel(
    dp$x, dp$y, paired = TRUE, test_method = "perm",
    alternative = "minimal.effect", mu = c(0.3, 0.7))))$p.value
  max(abs(c(proto - c(pk_e, pk_m), proto_p - c(pkp_e, pkp_m))))
})
cat("Max |prototype perm p - package perm p| (enumerated 5v5 and 8 pairs):",
    signif(max(val), 3), "\n")

# One replication --------

one_rep <- function(cell, alpha, R_perm) {
  paired <- cell$paired
  bounds <- c(cell$low, cell$high)
  d <- if (paired) gen_paired(cell$nx, cell$p, cell$rho) else
    gen_two(cell$nx, cell$ny, cell$p, cell$sdx)

  pkg <- function(method, alt) {
    tryCatch(suppressWarnings(suppressMessages(
      brunner_munzel(d$x, d$y, paired = paired, test_method = method,
                     alternative = alt, mu = bounds, alpha = alpha)))$p.value,
      error = function(e) NA_real_)
  }
  st <- if (paired) bmp_stats(d$x, d$y) else bm2_stats(d$x, d$y)
  pset <- if (paired) perm_set_paired(d$x, d$y, R_perm) else
    perm_set_two(d$x, d$y, R_perm)
  pp <- perm_tost(st, pset, bounds, "identity")
  pl <- perm_tost(st, pset, bounds, "logit")

  equ <- c(t = pkg("t", "equivalence"), logit = pkg("logit", "equivalence"),
           perm = pp[["equ"]], perm_logit = pl[["equ"]])
  met <- c(t = pkg("t", "minimal.effect"), logit = pkg("logit", "minimal.effect"),
           perm = pp[["met"]], perm_logit = pl[["met"]])
  c(setNames(equ <= alpha, paste0("equ_", names(equ))),
    setNames(met <= alpha, paste0("met_", names(met))))
}

# Cells --------

two <- expand.grid(design = c("7v7", "10v10", "15v15", "30v30", "10v20"),
                   sdx = c(1, 3), stringsAsFactors = FALSE)
two$nx <- as.integer(sub("v.*", "", two$design))
two$ny <- as.integer(sub(".*v", "", two$design))
two$paired <- FALSE
two$rho <- NA

pr <- data.frame(design = paste0(c(7, 10, 15, 30), " pairs"), sdx = 1,
                 nx = c(7L, 10L, 15L, 30L), ny = c(7L, 10L, 15L, 30L),
                 paired = TRUE, rho = 0.5)

designs <- rbind(two[, names(pr)], pr)

cells <- do.call(rbind, lapply(names(bounds_list), function(bn) {
  b <- bounds_list[[bn]]
  do.call(rbind, lapply(c("low", "centre", "high"), function(where) {
    p <- switch(where, low = b[1], centre = 0.5, high = b[2])
    cbind(designs, bounds = bn, low = b[1], high = b[2], truth = where, p = p)
  }))
}))
rownames(cells) <- NULL

# Run --------

cl <- parallel::makeCluster(n_cores)
invisible(parallel::clusterEvalQ(cl, devtools::load_all(".", quiet = TRUE)))
parallel::clusterExport(cl, c("one_rep", "bm2_stats", "bmp_stats", "clamp01",
                              "perm_set_two", "perm_set_paired", "perm_tost",
                              "gen_two", "gen_paired"))
parallel::clusterSetRNGStream(cl, 20260929)

t0 <- Sys.time()
results <- list()
for (k in seq_len(nrow(cells))) {
  cell <- as.list(cells[k, ])
  reps <- parallel::parSapply(cl, seq_len(nsim), function(i, cell, alpha, R_perm) {
    one_rep(cell, alpha, R_perm)
  }, cell = cell, alpha = alpha, R_perm = R_perm)
  results[[k]] <- cbind(cells[k, ], t(rowMeans(reps, na.rm = TRUE)),
                        n_na = sum(is.na(reps)))
  message(k, "/", nrow(cells), " ", cell$design, " sdx=", cell$sdx, " bounds=",
          cell$bounds, " truth=", cell$truth, " done (",
          format(round(Sys.time() - t0, 1)), ")")
  saveRDS(list(summary = do.call(rbind, results), nsim = nsim, R = R_perm,
               alpha = alpha, validation = max(val), complete = FALSE), out_file)
}
parallel::stopCluster(cl)

summary_tab <- do.call(rbind, results)
rownames(summary_tab) <- NULL

cat("\nnsim =", nsim, " R =", R_perm, " alpha =", alpha,
    " (MC SE near alpha ~", round(sqrt(alpha * (1 - alpha) / nsim), 4), ")\n\n")
print(summary_tab[, c("design", "sdx", "bounds", "truth",
                      grep("^(equ|met)_", names(summary_tab), value = TRUE))],
      digits = 3, row.names = FALSE)

saveRDS(list(summary = summary_tab, nsim = nsim, R = R_perm, alpha = alpha,
             validation = max(val), complete = TRUE, date = Sys.time()), out_file)
