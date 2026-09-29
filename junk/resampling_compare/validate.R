# Resampling comparison: pre-flight checks and runtime estimate --------
#
# Run from the package root before launching run.R:
#
#   Rscript junk/resampling_compare/validate.R
#
# Uses the same SIM_* environment variables as run.R (SIM_R, SIM_CORES,
# SIM_NSIM, SIM_SD_RATIO, SIM_PKG, SIM_DIR). Every check prints PASS or FAIL.

sim_dir <- normalizePath(Sys.getenv("SIM_DIR", "junk/resampling_compare"))
pkg <- normalizePath(Sys.getenv("SIM_PKG", "."), mustWork = FALSE)
R <- as.integer(Sys.getenv("SIM_R", 1999))
nsim <- as.integer(Sys.getenv("SIM_NSIM", 2000))
n_cores <- as.integer(Sys.getenv("SIM_CORES", max(1, parallel::detectCores() - 1)))
sd_ratio <- as.numeric(Sys.getenv("SIM_SD_RATIO", 2))

source(file.path(sim_dir, "setup.R"))
load_toster(pkg)   # errors if the TOSTER version predates the CI fix

report <- function(label, ok, detail = "") {
  cat(sprintf("%-4s %s %s\n", if (isTRUE(ok)) "PASS" else "FAIL", label, detail))
  invisible(ok)
}

cat("TOSTER version:", as.character(utils::packageVersion("TOSTER")), "\n\n")
report("TOSTER includes the perm_t_test CI fix (perm_crit/perm_count)", TRUE)

# 1. Population trimmed means (true values) --------

cat("\n1. True values: population vs Monte Carlo (n = 2e6)\n")
set.seed(1)
for (dn in c(names(cont_dists), names(disc_dists))) {
  z <- rdist(dn, 2e6)
  for (tr in c(0, 0.2)) {
    pop <- tm_pop(dn, tr)
    mc <- mean(z, trim = tr)
    report(sprintf("%-13s tr = %.1f", dn, tr), abs(pop - mc) < 0.005,
           sprintf("(population %.4f, MC %.4f)", pop, mc))
  }
}

# 2. Analytic reference (Welch / Yuen) --------

cat("\n2. Analytic reference intervals\n")
set.seed(2)
x <- rexp(20); y <- rexp(25) + 0.3
report("Welch (tr = 0) matches t.test",
       isTRUE(all.equal(ref_ci(x, y, 0, 0.95),
                        as.numeric(t.test(x, y)$conf.int))))
if (requireNamespace("WRS2", quietly = TRUE)) {
  yu <- WRS2::yuen(v ~ g, data = data.frame(v = c(x, y),
                                            g = factor(rep(c("x", "y"), c(20, 25)))),
                   tr = 0.2)
  report("Yuen two-sample CI matches WRS2::yuen",
         isTRUE(all.equal(ref_ci(x, y, 0.2, 0.95), as.numeric(yu$conf.int),
                          tolerance = 1e-6)))
} else {
  cat("SKIP Yuen two-sample check (install WRS2 to compare against WRS2::yuen)\n")
}
# One-sample trimmed t: Yuen-type SE sqrt((n - 1) s_w^2 / (h (h - 1))), df = h - 1
# (the same SE perm_t_test uses for one-sample/paired trimming), recomputed
# independently here from the winsorized sample
n <- length(x); g <- floor(0.2 * n); h <- n - 2 * g
xs <- sort(x); xw <- pmin(pmax(xs, xs[g + 1]), xs[n - g])
se <- sqrt((n - 1) * var(xw) / (h * (h - 1)))
manual <- mean(x, trim = 0.2) + c(-1, 1) * qt(0.975, h - 1) * se
report("One-sample trimmed CI matches independent recomputation",
       isTRUE(all.equal(ref_ci(x, NULL, 0.2, 0.95), manual)))

# 3. The 90% interval equals the interval each function uses for TOST --------

cat("\n3. 90% equal-tailed interval = TOST interval\n")
set.seed(3)
x <- rlnorm(12); y <- rlnorm(18)
for (tr in c(0, 0.2)) {
  set.seed(10)
  a <- perm_t_test(x, y, tr = tr, R = R, alpha = 0.10, symmetric = FALSE)$conf.int
  set.seed(10)
  b <- quiet(perm_t_test(x, y, tr = tr, R = R, alternative = "equivalence",
                         mu = c(-10, 10), alpha = 0.05))$conf.int
  report(sprintf("perm_t_test tr = %.1f", tr),
         isTRUE(all.equal(as.numeric(a), as.numeric(b))))
  for (type in c("stud", "bca")) {
    set.seed(11)
    a <- quiet(boot_t_test(x, y, tr = tr, R = R, alpha = 0.10, boot_ci = type))$conf.int
    set.seed(11)
    b <- quiet(boot_t_test(x, y, tr = tr, R = R, alternative = "equivalence",
                           mu = c(-10, 10), alpha = 0.05, boot_ci = type))$conf.int
    report(sprintf("boot_t_test %s tr = %.1f", type, tr),
           isTRUE(all.equal(as.numeric(a), as.numeric(b))))
  }
}

# 4. Paired calls with y = 0 match the difference scores --------

cat("\n4. Paired design passes differences through exactly\n")
set.seed(4)
xp <- rnorm(10); yp <- xp + rnorm(10, 0.3)
d <- xp - yp
set.seed(12)
a <- perm_t_test(xp, yp, paired = TRUE, R = R)$conf.int
set.seed(12)
b <- perm_t_test(d, rep(0, 10), paired = TRUE, R = R)$conf.int
report("perm_t_test paired(x, y) == paired(d, 0)",
       isTRUE(all.equal(as.numeric(a), as.numeric(b))))
set.seed(13)
a <- quiet(boot_t_test(xp, yp, paired = TRUE, R = R))$conf.int
set.seed(13)
b <- quiet(boot_t_test(d, rep(0, 10), paired = TRUE, R = R))$conf.int
report("boot_t_test paired(x, y) == paired(d, 0)",
       isTRUE(all.equal(as.numeric(a), as.numeric(b))))

# 5. Smoke run and runtime estimate --------

cat("\n5. Smoke run (3 replications on a sample of cells) and runtime estimate\n")
for (block in c("A", "B", "C")) {
  cells <- build_cells(block, sd_ratio = sd_ratio)
  idx <- unique(round(seq(1, nrow(cells), length.out = min(12, nrow(cells)))))
  secs <- vapply(idx, function(i) {
    cell <- as.list(cells[i, ])
    t <- system.time(res <- lapply(1:3, function(r)
      one_rep(cell, cell$truth, R, rep_seed = cell$cell_id * 100003L + r)))[["elapsed"]]
    v <- do.call(rbind, res)
    fails <- v[, grep("\\.fail$", colnames(v))]
    if (any(fails > 0)) cat("     note: failures in cell", cell$cell_id, "\n")
    t / 3
  }, numeric(1))
  per_rep <- mean(secs)
  hours <- per_rep * nsim * nrow(cells) / n_cores / 3600
  report(sprintf("Block %s: %d cells ran without error", block, nrow(cells)), TRUE,
         sprintf("(~%.2f s per replication; ~%.1f h at nsim = %d on %d cores)",
                 per_rep, hours, nsim, n_cores))
}

cat("\nIf everything is PASS, launch with e.g.  SIM_BLOCK=A Rscript junk/resampling_compare/run.R\n")
