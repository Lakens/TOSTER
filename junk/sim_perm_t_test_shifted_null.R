# Simulation: shifted vs unshifted permutation null for perm_t_test TOST --------
#
# Question: would a perm_null = "shifted" option in perm_t_test() improve
# equivalence (TOST) and minimal effect testing?
#
# At a null value theta (e.g., an equivalence bound):
#   - unshifted (current perm_t_test): observed statistic (d - theta) / SE is
#     compared with the permutation distribution of the UNSHIFTED data
#     (two-sample: pooled x and y; paired/one-sample: sign flips of the
#     differences centred at their mean)
#   - shifted: the data are shifted so the null value becomes zero before
#     permuting (two-sample: x - theta and y; paired/one-sample: sign flips of
#     d - theta). Exact under a pure location shift of size theta (two-sample)
#     or symmetry of the differences about theta (paired).
#   - welch: stats::t.test reference.
# Both permutation versions use the SAME permutation indices per dataset.
#
# Truth is placed at a bound (theta = +B or -B), so:
#   - TOST Type I error at the bound = the one-sided test in the equivalence
#     direction (should be <= alpha)
#   - MET Type I error at the bound = the one-sided test in the opposite
#     direction (should be <= alpha)
#   - their sum = 1 - coverage of the 1 - 2*alpha CI at that bound
# Truth at 0 gives TOST power (both one-sided tests reject).
#
# Results are summarized in junk/perm_t_test_shifted_null.md

# Settings --------

devtools::load_all(".", quiet = TRUE)

nsim <- as.integer(Sys.getenv("SIM_NSIM", 4000))
R_perm <- 999
alpha <- 0.05
n_cores <- max(1, parallel::detectCores() - 1)
out_file <- "junk/sim_perm_t_test_shifted_null_results.rds"

# Shapes: mean 0, SD 1 --------

draw <- list(
  normal = function(n) rnorm(n),
  skewed = function(n) rexp(n) - 1,
  heavy  = function(n) rt(n, df = 3) / sqrt(3)
)

# Permutation machinery (vectorized) --------

# Two-sample: index sets for group x; enumerate when choose(N, nx) <= R
two_perm_mask <- function(nx, ny, R) {
  N <- nx + ny
  exact <- choose(N, nx) <= R
  idx <- if (exact) t(utils::combn(N, nx)) else
    t(replicate(R, sample.int(N, nx)))
  B <- nrow(idx)
  mask <- matrix(FALSE, B, N)
  mask[cbind(rep(seq_len(B), nx), as.vector(idx))] <- TRUE
  list(mask = mask, exact = exact)
}

# Welch-studentized statistic for each permutation of pooled data z
two_perm_T <- function(z, mask, nx) {
  N <- length(z); ny <- N - nx
  Z <- matrix(z, nrow(mask), N, byrow = TRUE)
  Zc <- Z * mask; Zo <- Z * !mask
  mx <- rowSums(Zc) / nx; my <- rowSums(Zo) / ny
  vx <- (rowSums(Zc^2) - nx * mx^2) / (nx - 1)
  vy <- (rowSums(Zo^2) - ny * my^2) / (ny - 1)
  se <- sqrt(pmax(vx / nx + vy / ny, 0))
  ifelse(se > 1e-12, (mx - my) / se, 0)
}

# Paired / one-sample: sign matrices; enumerate when 2^n <= R
sign_matrix <- function(n, R) {
  exact <- 2^n <= R
  S <- if (exact) as.matrix(expand.grid(rep(list(c(-1, 1)), n))) else
    matrix(sample(c(-1, 1), n * R, replace = TRUE), nrow = R)
  list(S = S, exact = exact)
}

# Studentized mean for each sign pattern applied to |v| (as in perm_t_test)
paired_perm_T <- function(v, S) {
  n <- length(v)
  P <- S * matrix(abs(v), nrow(S), n, byrow = TRUE)
  m <- rowMeans(P)
  s <- sqrt(pmax((rowSums(P^2) - n * m^2) / (n - 1), 0) / n)
  ifelse(s > 1e-12, m / s, 0)
}

# One-sided p-values (package counting rules, with tie tolerance)
one_sided <- function(Tperm, t_obs, exact) {
  R <- length(Tperm)
  pm <- if (exact) "exact" else "plusone"
  c(less = TOSTER:::compute_perm_pval(TOSTER:::perm_count(Tperm, t_obs, "le"), R, pm),
    greater = TOSTER:::compute_perm_pval(TOSTER:::perm_count(Tperm, t_obs, "ge"), R, pm))
}

# p-values for testing a null value theta: unshifted, shifted, welch
test_two <- function(x, y, theta, pm) {
  nx <- length(x)
  se <- sqrt(var(x) / nx + var(y) / length(y))
  t_obs <- (mean(x) - mean(y) - theta) / se
  un <- one_sided(two_perm_T(c(x, y), pm$mask, nx), t_obs, pm$exact)
  sh <- one_sided(two_perm_T(c(x - theta, y), pm$mask, nx), t_obs, pm$exact)
  we <- c(less = t.test(x, y, mu = theta, alternative = "less")$p.value,
          greater = t.test(x, y, mu = theta, alternative = "greater")$p.value)
  list(unshifted = un, shifted = sh, welch = we)
}

test_paired <- function(d, theta, sm) {
  n <- length(d)
  t_obs <- (mean(d) - theta) / (sd(d) / sqrt(n))
  un <- one_sided(paired_perm_T(d - mean(d), sm$S), t_obs, sm$exact)
  sh <- one_sided(paired_perm_T(d - theta, sm$S), t_obs, sm$exact)
  we <- c(less = t.test(d, mu = theta, alternative = "less")$p.value,
          greater = t.test(d, mu = theta, alternative = "greater")$p.value)
  list(unshifted = un, shifted = sh, welch = we)
}

# Validation: prototype "unshifted" equals the package --------

set.seed(1)
val <- replicate(10, {
  x <- rexp(6) - 1 + 0.4; y <- rexp(6) - 1
  pm <- two_perm_mask(6, 6, 999)        # choose(12, 6) = 924: enumerated
  proto <- test_two(x, y, 0.5, pm)$unshifted
  pk <- sapply(c("less", "greater"), function(a) suppressMessages(
    perm_t_test(x, y, alternative = a, mu = 0.5, R = 999))$p.value)
  d <- rexp(8) - 1 + 0.3
  sm <- sign_matrix(8, 999)             # 2^8 = 256: enumerated
  proto_p <- test_paired(d, 0.5, sm)$unshifted
  pkp <- sapply(c("less", "greater"), function(a) suppressMessages(
    perm_t_test(d, alternative = a, mu = 0.5, R = 999))$p.value)
  max(abs(c(proto - pk, proto_p - pkp)))
})
cat("Max |prototype unshifted p - package p| (enumerated two-sample and one-sample):",
    signif(max(val), 3), "\n")

# One replication --------

one_rep <- function(cell, alpha, R_perm) {
  B <- cell$B
  f <- draw[[cell$shape]]
  if (cell$design == "two") {
    x <- f(cell$nx) * cell$sdx + cell$truth
    y <- f(cell$ny)
    pm <- two_perm_mask(cell$nx, cell$ny, R_perm)
    at <- function(theta) test_two(x, y, theta, pm)
  } else {
    d <- f(cell$nx) + cell$truth
    sm <- sign_matrix(cell$nx, R_perm)
    at <- function(theta) test_paired(d, theta, sm)
  }

  methods <- c("unshifted", "shifted", "welch")
  if (cell$truth > 0) {
    r <- at(B)
    tost <- sapply(methods, function(m) r[[m]][["less"]] <= alpha)
    met <- sapply(methods, function(m) r[[m]][["greater"]] <= alpha)
  } else if (cell$truth < 0) {
    r <- at(-B)
    tost <- sapply(methods, function(m) r[[m]][["greater"]] <= alpha)
    met <- sapply(methods, function(m) r[[m]][["less"]] <= alpha)
  } else {
    lo <- at(-B); hi <- at(B)
    tost <- sapply(methods, function(m)
      max(lo[[m]][["greater"]], hi[[m]][["less"]]) <= alpha)
    met <- setNames(rep(NA, 3), methods)
  }
  c(setNames(tost, paste0("tost_", methods)), setNames(met, paste0("met_", methods)))
}

# Cells --------

two <- expand.grid(shape = c("normal", "skewed", "heavy"),
                   n = c(6, 10, 20), B = c(0.5, 1, 2),
                   stringsAsFactors = FALSE)
two$nx <- two$n; two$ny <- two$n; two$sdx <- 1; two$scenario <- "equal shape"

het <- expand.grid(shape = "normal", n = c(8, 16), B = c(0.5, 1, 2),
                   stringsAsFactors = FALSE)
het$nx <- het$n; het$ny <- 2 * het$n; het$sdx <- 3
het$scenario <- "heteroskedastic (small group SD 3)"

two <- rbind(two, het)
two$design <- "two"

pr <- expand.grid(shape = c("normal", "skewed"), n = c(8, 15, 30),
                  B = c(0.5, 1), stringsAsFactors = FALSE)
pr$nx <- pr$n; pr$ny <- pr$n; pr$sdx <- 1; pr$scenario <- "paired differences"
pr$design <- "paired"

base <- rbind(two, pr)
cells <- do.call(rbind, lapply(c(-1, 0, 1), function(s) {
  b <- base; b$truth <- s * b$B; b$where <- c("lower", "centre", "upper")[s + 2]; b
}))
rownames(cells) <- NULL

# Run --------

cl <- parallel::makeCluster(n_cores)
invisible(parallel::clusterEvalQ(cl, devtools::load_all(".", quiet = TRUE)))
parallel::clusterExport(cl, c("one_rep", "draw", "two_perm_mask", "two_perm_T",
                              "sign_matrix", "paired_perm_T", "one_sided",
                              "test_two", "test_paired"))
parallel::clusterSetRNGStream(cl, 20260929)

t0 <- Sys.time()
results <- list()
for (k in seq_len(nrow(cells))) {
  cell <- as.list(cells[k, ])
  reps <- parallel::parSapply(cl, seq_len(nsim), function(i, cell, alpha, R_perm) {
    one_rep(cell, alpha, R_perm)
  }, cell = cell, alpha = alpha, R_perm = R_perm)
  results[[k]] <- cbind(cells[k, ], t(rowMeans(reps)))
  message(k, "/", nrow(cells), " done (", format(round(Sys.time() - t0, 1)), ")")
  saveRDS(list(summary = do.call(rbind, results), nsim = nsim, R = R_perm,
               alpha = alpha, validation = max(val), complete = FALSE), out_file)
}
parallel::stopCluster(cl)

summary_tab <- do.call(rbind, results)
rownames(summary_tab) <- NULL

cat("\nnsim =", nsim, " R =", R_perm, " alpha =", alpha,
    " (MC SE near alpha ~", round(sqrt(alpha * (1 - alpha) / nsim), 4), ")\n\n")
print(summary_tab[, c("design", "scenario", "shape", "nx", "ny", "B", "where",
                      grep("^(tost|met)_", names(summary_tab), value = TRUE))],
      digits = 3, row.names = FALSE)

saveRDS(list(summary = summary_tab, nsim = nsim, R = R_perm, alpha = alpha,
             validation = max(val), complete = TRUE, date = Sys.time()), out_file)
