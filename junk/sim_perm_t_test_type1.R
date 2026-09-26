# Simulation: perm_t_test type I error, exact vs randomized vs Welch --------
#
# Two-sided type I error at the true (mean) difference for:
#   - welch: stats::t.test (Welch)
#   - perm_random: perm_t_test with R = 999 sampled permutations (plusone)
#   - perm_exact: perm_t_test with full enumeration (b/R)
# All three methods are applied to the same simulated datasets. Sample sizes
# are kept small enough that full enumeration is feasible (~3000 permutations).
#
# Under exchangeability (identical distributions) both permutation versions are
# exact; with unequal variances/shapes they are only asymptotically valid.

# Settings --------

nsim <- as.integer(Sys.getenv("SIM_NSIM", 1000))  # replications per scenario
R_perm <- 999
alpha <- 0.05
n_cores <- max(1, parallel::detectCores() - 1)
out_file <- "junk/sim_perm_t_test_type1_results.rds"

# Scenarios --------

scenarios <- list(
  equal_normal_7v7 = list(
    gen = function() list(x = rnorm(7, 0, 1), y = rnorm(7, 0, 1)),
    truth = 0),
  bf_small_n_large_var_5v10 = list(
    gen = function() list(x = rnorm(5, 0, 4), y = rnorm(10, 0, 1)),
    truth = 0),
  bf_small_n_small_var_5v10 = list(
    gen = function() list(x = rnorm(5, 0, 1), y = rnorm(10, 0, 4)),
    truth = 0),
  equal_n_unequal_var_7v7 = list(
    gen = function() list(x = rnorm(7, 0, 1), y = rnorm(7, 0, 4)),
    truth = 0),
  skewed_unequal_5v10 = list(
    # lognormal(0, 1) mean = exp(0.5); exponential(1) mean = 1
    gen = function() list(x = rlnorm(5, 0, 1), y = rexp(10, 1)),
    truth = exp(0.5) - 1)
)

# One replication --------

one_rep <- function(sc, alpha, R_perm) {
  d <- sc$gen()
  q <- function(expr) suppressMessages(expr)
  c(welch = stats::t.test(d$x, d$y, mu = sc$truth)$p.value <= alpha,
    perm_random = q(perm_t_test(d$x, d$y, mu = sc$truth, R = R_perm,
                                keep_perm = FALSE))$p.value <= alpha,
    perm_exact = q(perm_t_test(d$x, d$y, mu = sc$truth, R = NULL,
                               keep_perm = FALSE))$p.value <= alpha)
}

# Run --------

cl <- parallel::makeCluster(n_cores)
invisible(parallel::clusterEvalQ(cl, devtools::load_all(".", quiet = TRUE)))
parallel::clusterExport(cl, "one_rep")
parallel::clusterSetRNGStream(cl, 20260926)

t0 <- Sys.time()
results <- lapply(names(scenarios), function(nm) {
  sc <- scenarios[[nm]]
  reps <- parallel::parSapply(cl, seq_len(nsim), function(i, sc, alpha, R_perm) {
    one_rep(sc, alpha, R_perm)
  }, sc = sc, alpha = alpha, R_perm = R_perm)
  message(nm, " done (", format(round(Sys.time() - t0, 1)), ")")
  data.frame(scenario = nm, t(rowMeans(reps)))
})
parallel::stopCluster(cl)

summary_tab <- do.call(rbind, results)

# Summary --------

mc_se <- sqrt(alpha * (1 - alpha) / nsim)
cat("\nnsim =", nsim, " R =", R_perm, " alpha =", alpha,
    " (MC SE near alpha ~", round(mc_se, 4), ")\n\n")
print(summary_tab, digits = 3, row.names = FALSE)

saveRDS(list(summary = summary_tab, nsim = nsim, R = R_perm, alpha = alpha,
             date = Sys.time()), out_file)
