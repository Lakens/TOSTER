# Merge Brunner-Munzel t vs perm sample-size results --------
#
# Combines the original run of sim_brunner_munzel_t_vs_perm_n.R with the reruns
# made after the floating point tie fix (tolerant permutation counts). Rerun
# cells replace the corresponding original cells:
#   - all cells with n <= 10 (SIM_NMAX=10)
#   - all identical_ordinal cells (tied data are affected at every n)
# Continuous-data cells with n > 10 are kept from the original run; with
# continuous data and sampled permutations, exact ties essentially never occur,
# so the fix does not change them.
#
# Prints markdown tables used in junk/bm_t_vs_perm_sample_size.md

orig <- readRDS("junk/sim_brunner_munzel_t_vs_perm_n_results.rds")
small <- readRDS("junk/sim_bm_t_vs_perm_n_rerun_small.rds")
ord <- readRDS("junk/sim_bm_t_vs_perm_n_rerun_ordinal.rds")
nsim <- orig$nsim

key <- function(d) paste(d$dist, d$ratio, d$n)
rerun <- rbind(small$summary, ord$summary)
rerun <- rerun[!duplicated(key(rerun)), ]
merged <- rbind(orig$summary[!key(orig$summary) %in% key(rerun), ], rerun)
merged$source <- ifelse(key(merged) %in% key(rerun), "rerun", "original")
merged <- merged[order(merged$dist, merged$ratio, merged$n), ]

# Paired precision of the t - perm difference (McNemar-type SE)
merged$se_diff <- sqrt(merged$disagree / nsim)
merged$z_diff <- ifelse(merged$se_diff > 0, merged$diff / merged$se_diff, 0)

saveRDS(list(summary = merged, nsim = nsim, R = orig$R, alpha = orig$alpha),
        "junk/sim_brunner_munzel_t_vs_perm_n_merged.rds")

# Markdown tables --------

fmt <- function(x, d = 3) formatC(x, format = "f", digits = d)

md_table <- function(df, cols, header) {
  cat("| ", paste(header, collapse = " | "), " |\n", sep = "")
  cat("|", paste(rep("---", length(header)), collapse = "|"), "|\n", sep = "")
  for (i in seq_len(nrow(df))) {
    cat("| ", paste(sapply(cols, function(cl) df[[cl]][i]), collapse = " | "), " |\n", sep = "")
  }
  cat("\n")
}

for (dd in unique(merged$dist)) {
  cat("###", dd, "\n\n")
  d <- merged[merged$dist == dd, ]
  d$design <- paste0(d$nx, " v ", d$ny)
  d$t_f <- fmt(d$t); d$perm_f <- fmt(d$perm)
  d$diff_f <- paste0(fmt(d$diff), " (", fmt(d$z_diff, 1), ")")
  d$dis_f <- fmt(d$disagree); d$mad_f <- fmt(d$mad_p)
  md_table(d, c("ratio", "design", "t_f", "perm_f", "diff_f", "dis_f", "mad_f"),
           c("ratio", "n (x v y)", "t", "perm", "t - perm (z)", "disagree", "mean |p_t - p_perm|"))
}

# Wide summary of z for the paired difference
wide <- reshape(merged[, c("dist", "ratio", "n", "z_diff")],
                idvar = c("dist", "ratio"), timevar = "n", direction = "wide")
print(wide, digits = 2, row.names = FALSE)
