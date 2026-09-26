# Simulation: Brunner-Munzel type I error in small samples --------
#
# Compares brunner_munzel() test_method = "t", "logit", and "perm" on the same
# simulated datasets. In every scenario the true relative effect
# p = P(X > Y) + 0.5 * P(X = Y) equals 0.5, so every rejection of H0: p = 0.5
# is a Type I error.
#
# Scenarios cross distribution pairs (identical, unequal variance, unequal
# shape, skewed) with small sample sizes, including unbalanced designs in both
# directions (the first group listed is x).

# Settings --------

nsim <- as.integer(Sys.getenv("SIM_NSIM", 2000))  # replications per cell
R_perm <- 999
alpha <- 0.05
n_cores <- max(1, parallel::detectCores() - 1)
out_file <- "junk/sim_brunner_munzel_small_results.rds"

sizes <- list(c(7, 7), c(5, 10), c(10, 5), c(12, 12))

# Shifts giving p = 0.5 for asymmetric pairs --------

# p(c) = P(X - c > Y) = integral of F_Y(x - c) f_X(x) dx
solve_shift <- function(fX, FY, lower = 0, upper = Inf) {
  pfun <- function(c) {
    stats::integrate(function(x) FY(x - c) * fX(x), lower, upper)$value - 0.5
  }
  stats::uniroot(pfun, c(-20, 20), tol = 1e-10)$root
}

# X ~ Exp(1) - c vs Y ~ N(0, 1)
c_skew_norm <- solve_shift(function(x) dexp(x, 1), function(z) pnorm(z))
# X ~ 3 * Exp(1) - c vs Y ~ Exp(1)
c_skew_scale <- solve_shift(function(x) dexp(x, 1 / 3), function(z) pexp(z, 1))

# Distribution pairs --------
# Each returns a generator function(nx, ny) -> list(x, y). Closures are built
# with local() so the shift constants travel with them to the workers.

dists <- list(
  identical_normal = local(function(nx, ny) {
    list(x = rnorm(nx), y = rnorm(ny))
  }),
  identical_ordinal = local({
    # 5-point ordinal with many ties; exchangeable
    pr <- c(0.1, 0.2, 0.4, 0.2, 0.1)
    function(nx, ny) list(x = sample(1:5, nx, TRUE, pr), y = sample(1:5, ny, TRUE, pr))
  }),
  normal_sd3_vs_sd1 = local(function(nx, ny) {
    list(x = rnorm(nx, 0, 3), y = rnorm(ny, 0, 1))
  }),
  t3_vs_uniform = local(function(nx, ny) {
    # both symmetric about 0 -> p = 0.5; heavy tails vs bounded
    list(x = rt(nx, df = 3), y = runif(ny, -2, 2))
  }),
  skewed_vs_normal = local({
    shift <- c_skew_norm
    function(nx, ny) list(x = rexp(nx, 1) - shift, y = rnorm(ny))
  }),
  skewed_scale3_vs_1 = local({
    shift <- c_skew_scale
    function(nx, ny) list(x = 3 * rexp(nx, 1) - shift, y = rexp(ny, 1))
  })
)

# Monte Carlo check that p = 0.5 --------

set.seed(1)
p_check <- sapply(dists, function(g) {
  d <- g(2e5, 2e5)
  mean(d$x > d$y) + 0.5 * mean(d$x == d$y)
})
cat("Monte Carlo check of true relative effect (should be ~0.5):\n")
print(round(p_check, 4))

# One replication --------

one_rep <- function(gen, nx, ny, alpha, R_perm) {
  d <- gen(nx, ny)
  run <- function(method) {
    p <- tryCatch(
      suppressWarnings(suppressMessages(
        brunner_munzel(d$x, d$y, test_method = method, R = R_perm)
      ))$p.value,
      error = function(e) NA_real_)
    p <= alpha
  }
  c(t = run("t"), logit = run("logit"), perm = run("perm"))
}

# Run --------

cells <- expand.grid(dist = names(dists), size = seq_along(sizes),
                     stringsAsFactors = FALSE)

cl <- parallel::makeCluster(n_cores)
invisible(parallel::clusterEvalQ(cl, devtools::load_all(".", quiet = TRUE)))
parallel::clusterExport(cl, "one_rep")
parallel::clusterSetRNGStream(cl, 20260926)

t0 <- Sys.time()
results <- lapply(seq_len(nrow(cells)), function(k) {
  gen <- dists[[cells$dist[k]]]
  n <- sizes[[cells$size[k]]]
  reps <- parallel::parSapply(cl, seq_len(nsim), function(i, gen, n, alpha, R_perm) {
    one_rep(gen, n[1], n[2], alpha, R_perm)
  }, gen = gen, n = n, alpha = alpha, R_perm = R_perm)
  message(cells$dist[k], " ", n[1], "v", n[2], " done (",
          format(round(Sys.time() - t0, 1)), ")")
  data.frame(dist = cells$dist[k], nx = n[1], ny = n[2],
             t(rowMeans(reps, na.rm = TRUE)),
             n_na = sum(is.na(reps)))
})
parallel::stopCluster(cl)

summary_tab <- do.call(rbind, results)
summary_tab <- summary_tab[order(summary_tab$dist, summary_tab$nx + summary_tab$ny,
                                 summary_tab$nx), ]

# Summary --------

mc_se <- sqrt(alpha * (1 - alpha) / nsim)
cat("\nnsim =", nsim, " R =", R_perm, " alpha =", alpha,
    " (MC SE near alpha ~", round(mc_se, 4), ")\n\n")
print(summary_tab, digits = 3, row.names = FALSE)

saveRDS(list(summary = summary_tab, p_check = p_check, nsim = nsim,
             R = R_perm, alpha = alpha, sizes = sizes, date = Sys.time()),
        out_file)
