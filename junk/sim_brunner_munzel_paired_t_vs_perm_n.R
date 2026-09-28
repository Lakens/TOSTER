# Simulation: paired Brunner-Munzel, t vs logit vs perm by sample size --------
#
# Question: for the paired Brunner-Munzel test, how do test_method = "t",
# "logit", and "perm" compare in small samples, and from what n does "t" give
# essentially the same inference as "perm"?
#
# The relative effect p = P(X > Y) + 0.5 * P(X = Y) is defined on the marginal
# distributions (X and Y drawn independently from their marginals), so it does
# not depend on the within-pair correlation. Pairs are generated with a
# Gaussian copula: correlated normals (rho) are mapped to the target marginals
# through their quantile functions. In every scenario p = 0.5, so every
# rejection of H0: p = 0.5 is a Type I error. All methods are applied to the
# same datasets.
#
# The paired permutation test swaps x and y within pairs. It enumerates all 2^n
# swaps when n <= 13 (p_method "exact") and samples R swap patterns otherwise
# (p_method "plusone").
#
# Metrics per cell:
#   - t, logit, perm: two-sided Type I error rate
#   - diff: t - perm Type I error (paired, same datasets)
#   - disagree: proportion of datasets where t and perm reach different decisions
#   - mad_p: mean absolute difference between the t and perm p-values
#
# Results are summarized in junk/bm_paired_t_vs_perm_sample_size.md

# Settings --------

nsim <- as.integer(Sys.getenv("SIM_NSIM", 4000))  # replications per cell
R_perm <- 999
alpha <- 0.05
n_cores <- max(1, parallel::detectCores() - 1)
out_file <- Sys.getenv("SIM_OUT", "junk/sim_brunner_munzel_paired_t_vs_perm_n_results.rds")

n_grid <- c(5, 7, 10, 13, 15, 20, 30, 50)
rho_grid <- c(0, 0.5, 0.8)

# Shift giving p = 0.5 for the skewed pair --------

# p(c) = P(X - c > Y) = integral of F_Y(x - c) f_X(x) dx (marginals only)
solve_shift <- function(fX, FY, lower = 0, upper = Inf) {
  pfun <- function(c) {
    stats::integrate(function(x) FY(x - c) * fX(x), lower, upper)$value - 0.5
  }
  stats::uniroot(pfun, c(-20, 20), tol = 1e-10)$root
}
c_skew_scale <- solve_shift(function(x) dexp(x, 1 / 3), function(z) pexp(z, 1))

# Marginal quantile functions (x listed first) --------
# Each entry maps a pair of uniforms (u1, u2) to (x, y).

margins <- list(
  identical_normal = local(function(u1, u2) {
    list(x = qnorm(u1), y = qnorm(u2))
  }),
  identical_ordinal = local({
    # 5-point ordinal with many ties; identical marginals
    cp <- cumsum(c(0.1, 0.2, 0.4, 0.2, 0.1))
    q <- function(u) findInterval(u, cp[-length(cp)]) + 1
    function(u1, u2) list(x = q(u1), y = q(u2))
  }),
  normal_sd3_vs_sd1 = local(function(u1, u2) {
    list(x = qnorm(u1, 0, 3), y = qnorm(u2, 0, 1))
  }),
  skewed_scale3_vs_1 = local({
    shift <- c_skew_scale
    function(u1, u2) list(x = 3 * qexp(u1) - shift, y = qexp(u2))
  })
)

# Gaussian copula pairs with correlation rho
gen_pairs <- function(margin, n, rho) {
  z1 <- rnorm(n)
  z2 <- rho * z1 + sqrt(1 - rho^2) * rnorm(n)
  margin(pnorm(z1), pnorm(z2))
}

# Monte Carlo check that p = 0.5 --------

set.seed(1)
p_check <- sapply(margins, function(m) {
  d <- m(runif(2e5), runif(2e5))  # independent draws from the marginals
  mean(d$x > d$y) + 0.5 * mean(d$x == d$y)
})
cat("Monte Carlo check of true relative effect (should be ~0.5):\n")
print(round(p_check, 4))

# One replication --------

one_rep <- function(margin, n, rho, R_perm) {
  d <- gen_pairs(margin, n, rho)
  pv <- function(method) {
    tryCatch(
      suppressWarnings(suppressMessages(
        brunner_munzel(d$x, d$y, paired = TRUE, test_method = method, R = R_perm)
      ))$p.value,
      error = function(e) NA_real_)
  }
  c(p_t = pv("t"), p_logit = pv("logit"), p_perm = pv("perm"))
}

# Run --------

cells <- expand.grid(n = n_grid, rho = rho_grid, dist = names(margins),
                     stringsAsFactors = FALSE)
# optional subset for timing, e.g. SIM_CELLS=1:3
if (nzchar(Sys.getenv("SIM_CELLS"))) {
  cells <- cells[eval(parse(text = Sys.getenv("SIM_CELLS"))), ]
}

cl <- parallel::makeCluster(n_cores)
invisible(parallel::clusterEvalQ(cl, devtools::load_all(".", quiet = TRUE)))
parallel::clusterExport(cl, c("one_rep", "gen_pairs"))
parallel::clusterSetRNGStream(cl, 20260928)

t0 <- Sys.time()
results <- list()
for (k in seq_len(nrow(cells))) {
  margin <- margins[[cells$dist[k]]]
  n <- cells$n[k]
  rho <- cells$rho[k]
  reps <- parallel::parSapply(cl, seq_len(nsim), function(i, margin, n, rho, R_perm) {
    one_rep(margin, n, rho, R_perm)
  }, margin = margin, n = n, rho = rho, R_perm = R_perm)

  ok <- stats::complete.cases(t(reps[c("p_t", "p_perm"), , drop = FALSE]))
  rej_t <- reps["p_t", ok] <= alpha
  rej_p <- reps["p_perm", ok] <= alpha
  results[[k]] <- data.frame(
    dist = cells$dist[k], rho = rho, n = n,
    t = mean(rej_t),
    logit = mean(reps["p_logit", ] <= alpha, na.rm = TRUE),
    perm = mean(rej_p),
    diff = mean(rej_t) - mean(rej_p),
    disagree = mean(rej_t != rej_p),
    mad_p = mean(abs(reps["p_t", ok] - reps["p_perm", ok])),
    n_na = sum(!ok))

  message(cells$dist[k], " n=", n, " rho=", rho, " done (",
          format(round(Sys.time() - t0, 1)), ")")
  # save partial results as we go
  saveRDS(list(summary = do.call(rbind, results), p_check = p_check,
               nsim = nsim, R = R_perm, alpha = alpha, complete = FALSE),
          out_file)
}
parallel::stopCluster(cl)

summary_tab <- do.call(rbind, results)
summary_tab$se_diff <- sqrt(summary_tab$disagree / nsim)
summary_tab$z_diff <- ifelse(summary_tab$se_diff > 0,
                             summary_tab$diff / summary_tab$se_diff, 0)

# Smallest n from which t is practically equivalent to perm --------
# |t - perm| <= tol and disagreement <= 1% at this and every larger n

equiv_n <- function(tol) {
  do.call(rbind, lapply(split(summary_tab, list(summary_tab$dist, summary_tab$rho),
                              drop = TRUE), function(d) {
    d <- d[order(d$n), ]
    ok <- abs(d$diff) <= tol & d$disagree <= 0.01
    idx <- which(rev(cumprod(rev(ok))) == 1)
    data.frame(dist = d$dist[1], rho = d$rho[1],
               n = if (length(idx)) d$n[min(idx)] else NA)
  }))
}
eq <- merge(equiv_n(0.005), equiv_n(0.01), by = c("dist", "rho"),
            suffixes = c("_tol005", "_tol010"))

# Summary --------

mc_se <- sqrt(alpha * (1 - alpha) / nsim)
cat("\nnsim =", nsim, " R =", R_perm, " alpha =", alpha,
    " (MC SE near alpha ~", round(mc_se, 4), ")\n\n")
print(summary_tab[order(summary_tab$dist, summary_tab$rho, summary_tab$n), ],
      digits = 3, row.names = FALSE)
cat("\nSmallest n from which t is practically equivalent to perm:\n")
print(eq[order(eq$dist, eq$rho), ], row.names = FALSE)

saveRDS(list(summary = summary_tab, equiv_n = eq, p_check = p_check,
             nsim = nsim, R = R_perm, alpha = alpha, complete = TRUE,
             date = Sys.time()), out_file)
