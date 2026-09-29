# Resampling comparison: summarize results --------
#
# Run from the package root after one or more blocks have finished:
#
#   Rscript junk/resampling_compare/summarize.R
#
# Reads <SIM_OUT_DIR>/results_block_*.rds and writes:
#   <SIM_OUT_DIR>/summary_long.csv   one row per cell x method, all metrics
#   <SIM_OUT_DIR>/summary_tables.md  markdown tables for the report
#
# Metrics (see setup.R):
#   type1      two-sided Type I error of the standard test at the true value
#              (1 - 95% coverage); nominal .05
#   one_sided  larger of err_below / err_above: the one-sided alpha = .05 error,
#              i.e. the TOST / minimal effect Type I error if the true value
#              were an equivalence bound (worst direction); nominal .05
#   width95    mean width of the 95% CI
#   power0     rejection of H0: 0 by the 95% CI (power when truth != 0)
#   tost05, tost1  TOST rejection with bounds +/-0.5, +/-1 (power when truth = 0)
#   fail       proportion of data sets where the method failed

sim_dir <- normalizePath(Sys.getenv("SIM_DIR", "junk/resampling_compare"))
out_dir <- Sys.getenv("SIM_OUT_DIR", file.path(sim_dir, "results"))
files <- list.files(out_dir, pattern = "^results_block_[ABC]\\.rds$", full.names = TRUE)
if (length(files) == 0) stop("No results_block_*.rds files in ", out_dir)

res <- lapply(files, readRDS)
wide <- do.call(rbind, lapply(res, `[[`, "summary"))
nsim <- res[[1]]$settings$nsim
methods <- c("perm", "boot_stud", "boot_bca", "ref")
labels <- c(perm = "perm", boot_stud = "boot stud", boot_bca = "boot BCa",
            ref = "Welch/Yuen")

cat("Blocks found:", paste(sub(".*_([ABC])\\.rds$", "\\1", files), collapse = ", "), "\n")
for (r in res) {
  cat(sprintf("  block %s: %d of %d cells%s\n", r$summary$block[1], nrow(r$summary),
              nrow(r$cells), if (isTRUE(r$complete)) " (complete)" else ""))
}

# Long format --------

id_cols <- c("block", "cell_id", "design", "structure", "layout", "xdist", "ydist",
             "sdx", "shift", "tr", "nx", "ny", "truth", "n_valid")
long <- do.call(rbind, lapply(methods, function(m) {
  g <- function(x) wide[[paste0(m, ".", x)]]
  data.frame(wide[, id_cols], method = m,
             fail = g("fail"),
             type1 = 1 - g("cov95"),
             err_below = g("err_below"), err_above = g("err_above"),
             one_sided = pmax(g("err_below"), g("err_above")),
             width95 = g("width95"), power0 = g("rej0"),
             tost05 = g("tost05"), tost1 = g("tost1"),
             stringsAsFactors = FALSE)
}))
long <- long[order(long$block, long$cell_id, match(long$method, methods)), ]
write.csv(long, file.path(out_dir, "summary_long.csv"), row.names = FALSE)

# Markdown tables --------

mc_se <- sqrt(0.05 * 0.95 / nsim)
flag <- function(v) {
  # mark rates more than 2 MC SE above .05 (liberal) or below (conservative)
  s <- formatC(v, format = "f", digits = 3)
  ifelse(is.na(v), "NA",
         ifelse(v > 0.05 + 2 * mc_se, paste0("**", s, "**"),
                ifelse(v < 0.05 - 2 * mc_se, paste0("_", s, "_"), s)))
}

wide_by_method <- function(df, metric, fmt = flag) {
  keys <- unique(df[, c("cell_id", "block", "design", "structure", "layout",
                        "xdist", "ydist", "sdx", "tr", "nx", "ny", "truth")])
  for (m in methods) {
    v <- df[df$method == m, c("cell_id", metric)]
    keys[[labels[[m]]]] <- fmt(v[[metric]][match(keys$cell_id, v$cell_id)])
  }
  keys
}

md_table <- function(d, cols, header, con) {
  cat("| ", paste(header, collapse = " | "), " |\n", sep = "", file = con)
  cat("|", paste(rep("---", length(header)), collapse = "|"), "|\n", sep = "", file = con)
  for (i in seq_len(nrow(d))) {
    cat("| ", paste(vapply(cols, function(cl) as.character(d[[cl]][i]), ""),
                    collapse = " | "), " |\n", sep = "", file = con)
  }
  cat("\n", file = con)
}

scen_label <- function(d) {
  dist <- ifelse(d$design == "paired" | d$xdist == d$ydist | is.na(d$ydist),
                 d$xdist, paste(d$xdist, "vs", d$ydist))
  paste0(dist, ifelse(d$sdx != 1, paste0(" (SD x", d$sdx, ")"), ""))
}

con <- file(file.path(out_dir, "summary_tables.md"), "w")
cat("# Resampling comparison: summary tables\n\n", file = con)
cat(sprintf("nsim = %d per cell; MC SE near .05 = %.4f. Rates more than 2 MC SE above .05 are **bold** (liberal); more than 2 MC SE below are _italic_ (conservative).\n\n",
            nsim, mc_se), file = con)

err_cells <- long[!grepl("^power", long$structure), ]
for (blk in sort(unique(err_cells$block))) {
  for (des in c("two", "paired")) {
    sub <- err_cells[err_cells$block == blk & err_cells$design == des, ]
    if (nrow(sub) == 0) next
    for (metric in c("type1", "one_sided")) {
      w <- wide_by_method(sub, metric)
      w$scenario <- scen_label(w)
      w$n <- if (des == "paired") w$nx else paste0(w$nx, " v ", w$ny)
      w <- w[order(w$structure, w$scenario, w$layout, w$nx, w$ny), ]
      ttl <- if (metric == "type1") "two-sided Type I error at the true value" else
        "one-sided error (TOST / minimal effect Type I error at a bound, worst direction)"
      cat(sprintf("## Block %s, %s: %s\n\n", blk,
                  if (des == "two") "two-sample" else "paired", ttl), file = con)
      md_table(w, c("structure", "scenario", "layout", "n", "tr", labels),
               c("structure", "scenario", "layout", "n", "tr", unname(labels)), con)
    }
  }
}

pow <- long[grepl("^power", long$structure), ]
if (nrow(pow) > 0) {
  f3 <- function(v) formatC(v, format = "f", digits = 3)
  w <- wide_by_method(pow, "power0", fmt = f3)
  w$scenario <- scen_label(w)
  w$n <- ifelse(w$design == "paired", w$nx, paste0(w$nx, " v ", w$ny))
  cat("## Block C: power of the standard test (true shift 0.5)\n\n", file = con)
  md_table(w[order(w$design, w$scenario, w$tr, w$nx), ],
           c("design", "scenario", "n", "tr", labels),
           c("design", "scenario", "n", "tr", unname(labels)), con)
}

# TOST power from identical-distribution null cells (truth = 0)
tp <- long[long$structure %in% c("identical", "paired differences") &
             abs(long$truth) < 1e-12, ]
if (nrow(tp) > 0) {
  f3 <- function(v) formatC(v, format = "f", digits = 3)
  for (b in c("tost05", "tost1")) {
    w <- wide_by_method(tp, b, fmt = f3)
    w$scenario <- scen_label(w)
    w$n <- ifelse(w$design == "paired", w$nx, paste0(w$nx, " v ", w$ny))
    cat(sprintf("## TOST power at a true difference of 0, bounds %s\n\n",
                if (b == "tost05") "+/-0.5" else "+/-1"), file = con)
    md_table(w[order(w$block, w$design, w$scenario, w$nx, w$ny), ],
             c("block", "design", "scenario", "n", "tr", labels),
             c("block", "design", "scenario", "n", "tr", unname(labels)), con)
  }
}

# Failures
fl <- long[long$fail > 0, c("block", "cell_id", "design", "structure", "xdist",
                            "nx", "ny", "tr", "method", "fail")]
cat("## Method failures (proportion of data sets)\n\n", file = con)
if (nrow(fl) == 0) {
  cat("None.\n\n", file = con)
} else {
  fl$fail <- formatC(fl$fail, format = "f", digits = 4)
  md_table(fl, names(fl), names(fl), con)
}
close(con)

cat("Wrote", file.path(out_dir, "summary_long.csv"), "and",
    file.path(out_dir, "summary_tables.md"), "\n")
