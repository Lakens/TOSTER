# Simulation: boot_smd_calc type I error by boot_ci vs smd_calc --------
#
# Two independent groups, equal n per group, function defaults otherwise
# (var.equal = FALSE -> d_av, Hedges' bias correction). Data are standardized
# draws (mean 0, var 1) scaled by group SDs, so the population d_av is exact:
#   d_av = (mu_x - mu_y) / sqrt((sd_x^2 + sd_y^2) / 2)
#
# Hypotheses:
#   - two.sided at d = 0 (null.value = 0)
#   - equivalence at the bound: d = 0.5, bounds +/- 0.5
#     (hi_miss = CI upper < 0.5, i.e. the near-bound one-sided test)
#
# Each dataset is analysed with every boot_ci type from the same RNG state, so
# all CI types see identical bootstrap resamples.
#
# smd_calc comparators: smd_ci/test_method = z/z, t/t, and nct (default CI;
# its p-value uses test_method = "z", so p and CI need not agree).

# Settings --------

nsim <- as.integer(Sys.getenv("SIM_NSIM", 1000))
R_boot <- 999
alpha <- 0.05
n_cores <- max(1, parallel::detectCores() - 1)
out_file <- Sys.getenv("SIM_OUT", "junk/sim_boot_smd_calc_type1_results.rds")
dist_filter <- Sys.getenv("SIM_DIST", "")  # e.g. "skewed" to run one DGP
ns <- c(10, 20, 50)
boot_types <- c("stud", "basic", "perc", "bca")

# Data generation --------

std_draw <- function(n, dist) {
  if (dist == "skewed") {
    s <- 0.8
    m <- exp(s^2 / 2)
    v <- (exp(s^2) - 1) * exp(s^2)
    (rlnorm(n, 0, s) - m) / sqrt(v)
  } else {
    rnorm(n)
  }
}

gen_xy <- function(n, dist, d) {
  sds <- if (dist == "hetero") c(1, 3) else c(1, 1)
  shift <- d * sqrt(sum(sds^2) / 2)
  list(x = shift + sds[1] * std_draw(n, dist),
       y = sds[2] * std_draw(n, dist))
}

hyps <- list(
  two_sided = list(d = 0, alternative = "two.sided", null = 0),
  equivalence = list(d = 0.5, alternative = "equivalence", null = c(-0.5, 0.5))
)

cells <- expand.grid(dist = c("normal", "skewed", "hetero"),
                     hyp = names(hyps), n = ns, stringsAsFactors = FALSE)
if (nzchar(dist_filter)) cells <- cells[cells$dist %in% dist_filter, ]

# One replication --------

summarise_fit <- function(res, hyp, alpha) {
  if (inherits(res, "error")) return(c(rej = NA, ci_rej = NA, lo_miss = NA,
                                       hi_miss = NA, agree = NA, est = NA))
  ci <- as.numeric(res$conf.int)
  rej <- res$p.value <= alpha
  ci_rej <- if (hyp$alternative == "two.sided") {
    ci[1] > hyp$null || ci[2] < hyp$null
  } else {
    ci[1] > min(hyp$null) && ci[2] < max(hyp$null)
  }
  c(rej = rej, ci_rej = ci_rej, lo_miss = ci[1] > hyp$d,
    hi_miss = ci[2] < hyp$d, agree = rej == ci_rej,
    est = unname(res$estimate))
}

one_rep <- function(cell, hyps, alpha, R_boot, boot_types) {
  hyp <- hyps[[cell$hyp]]
  d <- gen_xy(cell$n, cell$dist, hyp$d)
  q <- function(expr) suppressMessages(suppressWarnings(expr))
  rng_state <- get(".Random.seed", envir = globalenv())

  out <- lapply(boot_types, function(bt) {
    assign(".Random.seed", rng_state, envir = globalenv())
    res <- tryCatch(q(
      boot_smd_calc(d$x, d$y, alternative = hyp$alternative,
                    null.value = hyp$null, alpha = alpha,
                    boot_ci = bt, R = R_boot)), error = function(e) e)
    summarise_fit(res, hyp, alpha)
  })
  names(out) <- paste0("boot_", boot_types)

  analytic <- list(smd_z = c("z", "z"), smd_t = c("t", "t"),
                   smd_nct = c("nct", "z"))
  for (nm in names(analytic)) {
    res <- tryCatch(q(
      smd_calc(d$x, d$y, alternative = hyp$alternative,
               null.value = hyp$null, alpha = alpha,
               smd_ci = analytic[[nm]][1], test_method = analytic[[nm]][2])),
      error = function(e) e)
    out[[nm]] <- summarise_fit(res, hyp, alpha)
  }

  unlist(out)
}

# Run --------

cl <- parallel::makeCluster(n_cores)
invisible(parallel::clusterEvalQ(cl, devtools::load_all(".", quiet = TRUE)))
parallel::clusterExport(cl, c("one_rep", "gen_xy", "std_draw",
                              "summarise_fit"))
parallel::clusterSetRNGStream(cl, 20261005)

# Per-cell checkpoint so an interrupted run can resume
ckpt_file <- sub("_results.rds$", "_partial.rds", out_file)
ckpt <- if (file.exists(ckpt_file)) readRDS(ckpt_file) else list()
if (!is.null(ckpt$nsim) && ckpt$nsim != nsim) ckpt <- list()
ckpt$nsim <- nsim

t0 <- Sys.time()
results <- lapply(seq_len(nrow(cells)), function(i) {
  cell <- cells[i, ]
  key <- paste(cell$dist, cell$hyp, cell$n, sep = "_")
  if (!is.null(ckpt[[key]])) return(ckpt[[key]])
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
               ci_rej = m[[paste0(mt, ".ci_rej")]],
               lo_miss = m[[paste0(mt, ".lo_miss")]],
               hi_miss = m[[paste0(mt, ".hi_miss")]],
               agree = m[[paste0(mt, ".agree")]],
               mean_est = m[[paste0(mt, ".est")]],
               n_fail = n_fail[[paste0(mt, ".rej")]])
  }))
  message(cell$dist, " / ", cell$hyp, " / n = ", cell$n, " done (",
          format(round(Sys.time() - t0, 1)), ")")
  ckpt[[key]] <<- tab
  saveRDS(ckpt, ckpt_file)
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
