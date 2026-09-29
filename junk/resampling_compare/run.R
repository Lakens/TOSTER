# Resampling comparison: run one block --------
#
# Run from the package root (the folder containing DESCRIPTION), e.g.
#
#   SIM_BLOCK=A Rscript junk/resampling_compare/run.R
#
# Settings (environment variables):
#   SIM_BLOCK     "A" (core, tr = 0), "B" (core, tr = 0.2), "C" (different
#                 shapes + power). Required.
#   SIM_NSIM      replications per cell (default 2000; MC SE ~ .005 near .05)
#   SIM_R         permutations / bootstrap resamples (default 1999, both methods)
#   SIM_CORES     worker processes (default: all cores minus one)
#   SIM_SD_RATIO  SD ratio for the "unequal SD" cells (default 2)
#   SIM_PKG       path to the TOSTER source (default "."); if it has no
#                 DESCRIPTION, the installed TOSTER is used
#   SIM_DIR       folder containing setup.R (default "junk/resampling_compare")
#   SIM_OUT_DIR   where results are written (default "<SIM_DIR>/results")
#
# Results are saved after every cell to <SIM_OUT_DIR>/results_block_<X>.rds.
# Re-running the same block resumes: completed cells are skipped. Results are
# reproducible regardless of the number of cores (per-replication seeds).

# Settings --------

block <- toupper(Sys.getenv("SIM_BLOCK", ""))
if (!block %in% c("A", "B", "C")) stop("Set SIM_BLOCK to A, B, or C")
nsim <- as.integer(Sys.getenv("SIM_NSIM", 2000))
R <- as.integer(Sys.getenv("SIM_R", 1999))
n_cores <- as.integer(Sys.getenv("SIM_CORES", max(1, parallel::detectCores() - 1)))
sd_ratio <- as.numeric(Sys.getenv("SIM_SD_RATIO", 2))
sim_dir <- normalizePath(Sys.getenv("SIM_DIR", "junk/resampling_compare"))
out_dir <- Sys.getenv("SIM_OUT_DIR", file.path(sim_dir, "results"))
pkg <- normalizePath(Sys.getenv("SIM_PKG", "."), mustWork = FALSE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_file <- file.path(out_dir, paste0("results_block_", block, ".rds"))

source(file.path(sim_dir, "setup.R"))
Sys.setenv(SIM_PKG = pkg)
load_toster(pkg)

cells <- build_cells(block, sd_ratio = sd_ratio)

# Resume --------

done <- list()
if (file.exists(out_file)) {
  prev <- readRDS(out_file)
  if (!identical(prev$settings[c("nsim", "R", "sd_ratio")],
                 list(nsim = nsim, R = R, sd_ratio = sd_ratio))) {
    stop("Existing results in ", out_file, " used different settings. ",
         "Move or delete the file, or rerun with the same settings.")
  }
  done <- split(prev$summary, prev$summary$cell_id)
}
todo <- which(!as.character(cells$cell_id) %in% names(done))

cat("Block", block, ":", nrow(cells), "cells,", length(todo), "to run;",
    "nsim =", nsim, " R =", R, " cores =", n_cores, "\n")
if (length(todo) == 0) quit(save = "no")

# Run --------

cl <- parallel::makeCluster(n_cores)
parallel::clusterExport(cl, c("sim_dir", "pkg"))
invisible(parallel::clusterEvalQ(cl, {
  source(file.path(sim_dir, "setup.R"))
  load_toster(pkg)
}))

settings <- list(nsim = nsim, R = R, sd_ratio = sd_ratio)
t0 <- Sys.time()
for (k in seq_along(todo)) {
  i <- todo[k]
  cell <- as.list(cells[i, ])
  reps <- parallel::parSapply(cl, seq_len(nsim), function(r, cell, R) {
    one_rep(cell, cell$truth, R, rep_seed = cell$cell_id * 100003L + r)
  }, cell = cell, R = R)

  # reps: one row per metric, one column per replication
  row <- cbind(cells[i, ], t(rowMeans(reps, na.rm = TRUE)),
               n_valid = min(rowSums(!is.na(reps))))
  done[[as.character(cell$cell_id)]] <- row
  saveRDS(list(summary = do.call(rbind, done), cells = cells, settings = settings,
               complete = length(done) == nrow(cells),
               date = Sys.time()), out_file)

  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  cat(sprintf("[%s] %d/%d  %s %s %s n=%dv%d tr=%.1f  (%.1f min elapsed, ~%.0f min left)\n",
              block, k, length(todo), cell$design, cell$structure, cell$xdist,
              cell$nx, cell$ny, cell$tr, elapsed,
              elapsed / k * (length(todo) - k)))
}

parallel::stopCluster(cl)
cat("Done. Results in", out_file, "\n")
