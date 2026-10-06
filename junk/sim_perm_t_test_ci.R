# Simulation: perm_t_test CI vs p-value alignment (issue #120) --------
#
# Compares the new test-inversion CI in perm_t_test() against the old
# percentile CI of the raw permuted differences (rebuilt from perm.eff).
#
# Metrics per scenario:
#   - cov_new / cov_old: CI coverage of the true (trimmed) mean difference
#   - type1: two-sided p-value rejection rate at the true value
#   - dis_new / dis_old: rate at which "true value outside CI" disagrees with
#     "p <= alpha" (should be 0 for the new CI)
#   - tost_rej: equivalence rejection rate when the true value sits on the
#     lower bound (should be <= alpha)
#   - tost_dis_new / tost_dis_old: disagreement between the TOST p-value decision
#     and "90% CI within bounds"
#
# Settings --------

nsim <- as.integer(Sys.getenv("SIM_NSIM", 1000))  # replications per scenario
R_perm <- 999      # permutations per test
alpha <- 0.05
n_cores <- max(1, parallel::detectCores() - 1)
out_file <- "junk/sim_perm_t_test_ci_results.rds"

# Scenarios --------

# Each scenario returns a list(x, y, paired) and records the true value
scenarios <- list(
  equal_normal = list(
    gen = function() list(x = rnorm(15, 0, 1), y = rnorm(15, 0, 1)),
    truth = 0, tr = 0, paired = FALSE),
  bf_small_n_large_var = list(
    gen = function() list(x = rnorm(8, 0, 4), y = rnorm(30, 0, 1)),
    truth = 0, tr = 0, paired = FALSE),
  bf_small_n_small_var = list(
    gen = function() list(x = rnorm(8, 0, 1), y = rnorm(30, 0, 4)),
    truth = 0, tr = 0, paired = FALSE),
  equal_n_unequal_var = list(
    gen = function() list(x = rnorm(12, 0, 1), y = rnorm(12, 0, 4)),
    truth = 0, tr = 0, paired = FALSE),
  skewed_unequal = list(
    # lognormal(0, 1) mean = exp(0.5); exponential(1) mean = 1
    gen = function() list(x = rlnorm(10, 0, 1), y = rexp(30, 1)),
    truth = exp(0.5) - 1, tr = 0, paired = FALSE),
  bf_trimmed = list(
    # symmetric distributions, so trimmed mean difference = mean difference
    gen = function() list(x = rnorm(10, 0, 4), y = rnorm(30, 0, 1)),
    truth = 0, tr = 0.2, paired = FALSE),
  one_sample_normal = list(
    gen = function() list(x = rnorm(15, 0, 2), y = NULL),
    truth = 0, tr = 0, paired = FALSE),
  paired_normal = list(
    gen = function() {
      x <- rnorm(15, 0, 1)
      list(x = x, y = x + rnorm(15, 0, 2))
    },
    truth = 0, tr = 0, paired = TRUE)
)

# One replication --------

one_rep <- function(sc, alpha, R_perm) {
  d <- sc$gen()
  truth <- sc$truth

  run <- function(alternative, mu) {
    suppressMessages(perm_t_test(d$x, d$y, paired = sc$paired, tr = sc$tr,
                                 alternative = alternative, mu = mu,
                                 alpha = alpha, R = R_perm))
  }

  # Two-sided test at the true value
  two <- run("two.sided", truth)
  ci_new <- two$conf.int
  ci_old <- quantile(two$perm.eff, c(alpha / 2, 1 - alpha / 2), names = FALSE)
  rej <- two$p.value <= alpha
  out_new <- truth < ci_new[1] || truth > ci_new[2]
  out_old <- truth < ci_old[1] || truth > ci_old[2]

  # Equivalence test with the true value on the lower bound
  sd_ref <- if (is.null(d$y)) sd(d$x) else sd(c(d$x, d$y))
  bounds <- c(truth, truth + 2 * sd_ref)
  eq <- run("equivalence", bounds)
  eq_rej <- eq$p.value <= alpha
  eq_ci_old <- quantile(eq$perm.eff, c(alpha, 1 - alpha), names = FALSE)
  in_new <- eq$conf.int[1] > bounds[1] && eq$conf.int[2] < bounds[2]
  in_old <- eq_ci_old[1] > bounds[1] && eq_ci_old[2] < bounds[2]

  c(cov_new = !out_new,
    cov_old = !out_old,
    type1 = rej,
    dis_new = out_new != rej,
    dis_old = out_old != rej,
    tost_rej = eq_rej,
    tost_dis_new = in_new != eq_rej,
    tost_dis_old = in_old != eq_rej)
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
    " (MC SE for a rate near alpha ~", round(mc_se, 4), ")\n\n")
print(summary_tab, digits = 3, row.names = FALSE)

saveRDS(list(summary = summary_tab, nsim = nsim, R = R_perm, alpha = alpha,
             date = Sys.time()), out_file)
