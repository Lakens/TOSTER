# Simulation: boot_cor_test type I error by boot_ci vs z_cor_test --------
#
# Pearson only. Data: y = b*x + e with x, e standardized (mean 0, var 1), so
# rho = b / sqrt(b^2 + 1) is known analytically.
#
# Hypotheses:
#   - two.sided at rho = 0 (null = 0)
#   - equivalence at the bound: rho = 0.3, bounds +/- 0.3
#     ("one_sided" = CI upper < 0.3, i.e. the near-bound one-sided test, which
#     should reject at ~alpha; TOST itself is conservative at small n)
#
# Each dataset is analysed with every boot_ci type using the same seed, so all
# CI types see identical bootstrap resamples (RNG state restored per call).
#
# Note: for Pearson, boot_ci = "stud" uses the constant SE 1/sqrt(n - 3) on the
# Fisher-z scale, so it is the basic bootstrap on the z scale.

# Settings --------

nsim <- as.integer(Sys.getenv("SIM_NSIM", 1000))
R_boot <- 999
alpha <- 0.05
n_cores <- max(1, parallel::detectCores() - 1)
out_file <- "junk/sim_boot_cor_test_type1_results.rds"
ns <- c(10, 20, 50)
boot_types <- c("bca", "stud", "basic", "perc")

# Data generation --------

gen_xy <- function(n, dist, b) {
  if (dist == "normal") {
    x <- rnorm(n)
    e <- rnorm(n)
  } else if (dist == "skewed") {
    x <- rexp(n) - 1
    e <- rt(n, df = 4) / sqrt(2)
  } else if (dist == "hetero") {
    x <- rnorm(n)
    e <- rnorm(n) * sqrt((1 + x^2) / 2)
  }
  list(x = x, y = b * x + e)
}

hyps <- list(
  two_sided = list(rho = 0, alternative = "two.sided", null = 0),
  equivalence = list(rho = 0.3, alternative = "equivalence", null = 0.3)
)

cells <- expand.grid(dist = c("normal", "skewed", "hetero"),
                     hyp = names(hyps), n = ns, stringsAsFactors = FALSE)

# One replication --------

summarise_fit <- function(res, rho, hyp, alpha) {
  if (inherits(res, "error")) return(c(rej = NA, lo_miss = NA, hi_miss = NA,
                                       agree = NA, est = NA))
  ci <- as.numeric(res$conf.int)
  rej <- res$p.value <= alpha
  ci_rej <- if (hyp$alternative == "two.sided") {
    ci[1] > hyp$null || ci[2] < hyp$null
  } else {
    ci[1] > -hyp$null && ci[2] < hyp$null
  }
  c(rej = rej, lo_miss = ci[1] > rho, hi_miss = ci[2] < rho,
    agree = rej == ci_rej, est = unname(res$estimate))
}

one_rep <- function(cell, hyps, alpha, R_boot, boot_types) {
  hyp <- hyps[[cell$hyp]]
  b <- hyp$rho / sqrt(1 - hyp$rho^2)
  d <- gen_xy(cell$n, cell$dist, b)
  rng_state <- get(".Random.seed", envir = globalenv())

  out <- lapply(boot_types, function(bt) {
    assign(".Random.seed", rng_state, envir = globalenv())
    res <- tryCatch(suppressWarnings(
      boot_cor_test(d$x, d$y, alternative = hyp$alternative,
                    method = "pearson", null = hyp$null, alpha = alpha,
                    boot_ci = bt, R = R_boot)), error = function(e) e)
    summarise_fit(res, hyp$rho, hyp, alpha)
  })
  names(out) <- paste0("boot_", boot_types)

  res_z <- tryCatch(suppressWarnings(
    z_cor_test(d$x, d$y, alternative = hyp$alternative, method = "pearson",
               null = hyp$null, alpha = alpha)), error = function(e) e)
  out$fisher_z <- summarise_fit(res_z, hyp$rho, hyp, alpha)

  unlist(out)
}

# Run --------

cl <- parallel::makeCluster(n_cores)
invisible(parallel::clusterEvalQ(cl, devtools::load_all(".", quiet = TRUE)))
parallel::clusterExport(cl, c("one_rep", "gen_xy", "summarise_fit"))
parallel::clusterSetRNGStream(cl, 20261005)

t0 <- Sys.time()
results <- lapply(seq_len(nrow(cells)), function(i) {
  cell <- cells[i, ]
  reps <- parallel::parSapply(cl, seq_len(nsim), function(j, cell, hyps, alpha,
                                                          R_boot, boot_types) {
    one_rep(cell, hyps, alpha, R_boot, boot_types)
  }, cell = cell, hyps = hyps, alpha = alpha, R_boot = R_boot,
  boot_types = boot_types)
  m <- rowMeans(reps, na.rm = TRUE)
  n_fail <- rowSums(is.na(reps))
  methods <- unique(sub("\\..*$", "", names(m)))
  tab <- do.call(rbind, lapply(methods, function(mt) {
    data.frame(dist = cell$dist, hyp = cell$hyp, n = cell$n, method = mt,
               rej = m[[paste0(mt, ".rej")]],
               lo_miss = m[[paste0(mt, ".lo_miss")]],
               hi_miss = m[[paste0(mt, ".hi_miss")]],
               agree = m[[paste0(mt, ".agree")]],
               mean_est = m[[paste0(mt, ".est")]],
               n_fail = n_fail[[paste0(mt, ".rej")]])
  }))
  message(cell$dist, " / ", cell$hyp, " / n = ", cell$n, " done (",
          format(round(Sys.time() - t0, 1)), ")")
  tab
})
parallel::stopCluster(cl)

summary_tab <- do.call(rbind, results)

# Summary --------

mc_se <- sqrt(alpha * (1 - alpha) / nsim)
cat("\nnsim =", nsim, " R =", R_boot, " alpha =", alpha,
    " (MC SE near alpha ~", round(mc_se, 4), ")\n\n")
print(summary_tab, digits = 3, row.names = FALSE)

saveRDS(list(summary = summary_tab, nsim = nsim, R = R_boot, alpha = alpha,
             date = Sys.time()), out_file)
