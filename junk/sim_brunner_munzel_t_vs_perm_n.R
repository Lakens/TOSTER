# Simulation: when does the Brunner-Munzel t method match the permutation test? --------
#
# Question: at what sample size does brunner_munzel(test_method = "t") give
# essentially the same inference as test_method = "perm"?
#
# In every scenario the true relative effect p = P(X > Y) + 0.5 * P(X = Y)
# equals 0.5, so every rejection of H0: p = 0.5 is a Type I error. Both methods
# are applied to the same datasets.
#
# Metrics per cell:
#   - t, perm: two-sided Type I error rate
#   - diff: t - perm Type I error (paired, same datasets)
#   - disagree: proportion of datasets where the two methods reach different
#     decisions at alpha
#   - mad_p: mean absolute difference between the two p-values
#
# Results are summarized in junk/bm_t_vs_perm_sample_size.md

# Settings --------

nsim <- as.integer(Sys.getenv("SIM_NSIM", 4000))  # replications per cell
R_perm <- 999
alpha <- 0.05
n_cores <- max(1, parallel::detectCores() - 1)
out_file <- Sys.getenv("SIM_OUT", "junk/sim_brunner_munzel_t_vs_perm_n_results.rds")

# Smaller group size; designs are balanced (n, n) and unbalanced (n, 2n), (2n, n)
n_grid <- c(5, 7, 10, 15, 20, 25, 30, 40, 50)
# optional subset, e.g. SIM_NMAX=10 to rerun only the small-n cells
if (nzchar(Sys.getenv("SIM_NMAX"))) n_grid <- n_grid[n_grid <= as.numeric(Sys.getenv("SIM_NMAX"))]
ratios <- list(`1:1` = c(1, 1), `1:2` = c(1, 2), `2:1` = c(2, 1))

# Shift giving p = 0.5 for the skewed pair --------

# p(c) = P(X - c > Y) = integral of F_Y(x - c) f_X(x) dx
solve_shift <- function(fX, FY, lower = 0, upper = Inf) {
  pfun <- function(c) {
    stats::integrate(function(x) FY(x - c) * fX(x), lower, upper)$value - 0.5
  }
  stats::uniroot(pfun, c(-20, 20), tol = 1e-10)$root
}
c_skew_scale <- solve_shift(function(x) dexp(x, 1 / 3), function(z) pexp(z, 1))

# Distribution pairs (x listed first) --------

dists <- list(
  identical_normal = local(function(nx, ny) {
    list(x = rnorm(nx), y = rnorm(ny))
  }),
  identical_ordinal = local({
    pr <- c(0.1, 0.2, 0.4, 0.2, 0.1)
    function(nx, ny) list(x = sample(1:5, nx, TRUE, pr), y = sample(1:5, ny, TRUE, pr))
  }),
  normal_sd3_vs_sd1 = local(function(nx, ny) {
    list(x = rnorm(nx, 0, 3), y = rnorm(ny, 0, 1))
  }),
  skewed_scale3_vs_1 = local({
    shift <- c_skew_scale
    function(nx, ny) list(x = 3 * rexp(nx, 1) - shift, y = rexp(ny, 1))
  })
)

# One replication --------

one_rep <- function(gen, nx, ny, alpha, R_perm) {
  d <- gen(nx, ny)
  pv <- function(method) {
    tryCatch(
      suppressWarnings(suppressMessages(
        brunner_munzel(d$x, d$y, test_method = method, R = R_perm)
      ))$p.value,
      error = function(e) NA_real_)
  }
  c(p_t = pv("t"), p_perm = pv("perm"))
}

# Run --------

cells <- expand.grid(n = n_grid, ratio = names(ratios), dist = names(dists),
                     stringsAsFactors = FALSE)
# optional subset of distributions, e.g. SIM_DISTS=identical_ordinal
if (nzchar(Sys.getenv("SIM_DISTS"))) {
  cells <- cells[cells$dist %in% strsplit(Sys.getenv("SIM_DISTS"), ",")[[1]], ]
}

cl <- parallel::makeCluster(n_cores)
invisible(parallel::clusterEvalQ(cl, devtools::load_all(".", quiet = TRUE)))
parallel::clusterExport(cl, "one_rep")
parallel::clusterSetRNGStream(cl, 20260926)

t0 <- Sys.time()
results <- list()
for (k in seq_len(nrow(cells))) {
  gen <- dists[[cells$dist[k]]]
  sz <- cells$n[k] * ratios[[cells$ratio[k]]]
  reps <- parallel::parSapply(cl, seq_len(nsim), function(i, gen, sz, alpha, R_perm) {
    one_rep(gen, sz[1], sz[2], alpha, R_perm)
  }, gen = gen, sz = sz, alpha = alpha, R_perm = R_perm)

  ok <- stats::complete.cases(t(reps))
  rej_t <- reps["p_t", ok] <= alpha
  rej_p <- reps["p_perm", ok] <= alpha
  results[[k]] <- data.frame(
    dist = cells$dist[k], ratio = cells$ratio[k], n = cells$n[k],
    nx = sz[1], ny = sz[2],
    t = mean(rej_t), perm = mean(rej_p), diff = mean(rej_t) - mean(rej_p),
    disagree = mean(rej_t != rej_p),
    mad_p = mean(abs(reps["p_t", ok] - reps["p_perm", ok])),
    n_na = sum(!ok))

  message(cells$dist[k], " ", sz[1], "v", sz[2], " done (",
          format(round(Sys.time() - t0, 1)), ")")
  # save partial results as we go
  saveRDS(list(summary = do.call(rbind, results), nsim = nsim, R = R_perm,
               alpha = alpha, complete = FALSE), out_file)
}
parallel::stopCluster(cl)

summary_tab <- do.call(rbind, results)

# Smallest n from which t stays in Bradley's stringent band --------
# Bradley (1978) stringent criterion: 0.9 * alpha to 1.1 * alpha
band <- c(0.9, 1.1) * alpha
stable_n <- do.call(rbind, lapply(split(summary_tab, list(summary_tab$dist, summary_tab$ratio), drop = TRUE),
  function(d) {
    d <- d[order(d$n), ]
    inside <- d$t >= band[1] & d$t <= band[2]
    # first n such that all n' >= n are inside the band
    idx <- which(rev(cumprod(rev(inside))) == 1)
    data.frame(dist = d$dist[1], ratio = d$ratio[1],
               t_stable_n = if (length(idx)) d$n[min(idx)] else NA)
  }))
rownames(stable_n) <- NULL

# Summary --------

mc_se <- sqrt(alpha * (1 - alpha) / nsim)
cat("\nnsim =", nsim, " R =", R_perm, " alpha =", alpha,
    " (MC SE near alpha ~", round(mc_se, 4), ")\n\n")
print(summary_tab, digits = 3, row.names = FALSE)
cat("\nSmallest n (smaller group) from which the t method stays within",
    "Bradley's stringent band [", band[1], ",", band[2], "]:\n")
print(stable_n, row.names = FALSE)

saveRDS(list(summary = summary_tab, stable_n = stable_n, nsim = nsim,
             R = R_perm, alpha = alpha, band = band, complete = TRUE,
             date = Sys.time()), out_file)
