# test_paired_score.R
# Manual verification of paired/one-sample score method against wilcox.test

library(TOSTER)

cat("=============================================================\n")
cat("  Paired Score Method: Manual Verification Against wilcox.test\n")
cat("=============================================================\n\n")

# Helper to compare ses_calc score vs wilcox.test
compare_paired <- function(x, y = NULL, mu = 0, label = "") {
  cat("---", label, "---\n")

  if (is.null(y)) {
    cat("  Design: one-sample, mu =", mu, "\n")
    cat("  x:", head(x, 10), if (length(x) > 10) "...", "\n")
    cat("  n =", length(x), "\n")
  } else {
    cat("  Design: paired\n")
    cat("  x:", head(x, 10), if (length(x) > 10) "...", "\n")
    cat("  y:", head(y, 10), if (length(y) > 10) "...", "\n")
    cat("  n =", length(x), "\n")
    d <- x - y
    cat("  # zero diffs:", sum(d == 0), "\n")
  }
  cat("\n")

  # --- Two-sided, no cc ---
  if (is.null(y)) {
    wt <- wilcox.test(x, mu = mu, exact = FALSE, correct = FALSE)
  } else {
    wt <- wilcox.test(x, y, paired = TRUE, exact = FALSE, correct = FALSE)
  }
  res <- ses_calc(x = x, y = y, paired = !is.null(y), ses = "cstat", mu = mu,
                  se_method = "score", correct = FALSE,
                  alternative = "two.sided", null.value = 0.5)

  cat("  Two-sided (no cc):\n")
  cat("    wilcox.test  V =", wt$statistic, " p =", format(wt$p.value, digits = 8), "\n")
  cat("    ses_calc     z =", format(res$statistic, digits = 5),
      " p =", format(res$p.value, digits = 8), "\n")
  cat("    p-value match:", all.equal(res$p.value, wt$p.value, tolerance = 1e-6), "\n")

  # --- Two-sided, with cc ---
  if (is.null(y)) {
    wt_cc <- wilcox.test(x, mu = mu, exact = FALSE, correct = TRUE)
  } else {
    wt_cc <- wilcox.test(x, y, paired = TRUE, exact = FALSE, correct = TRUE)
  }
  res_cc <- ses_calc(x = x, y = y, paired = !is.null(y), ses = "cstat", mu = mu,
                     se_method = "score", correct = TRUE,
                     alternative = "two.sided", null.value = 0.5)

  cat("  Two-sided (with cc):\n")
  cat("    wilcox.test  V =", wt_cc$statistic, " p =", format(wt_cc$p.value, digits = 8), "\n")
  cat("    ses_calc     z =", format(res_cc$statistic, digits = 5),
      " p =", format(res_cc$p.value, digits = 8), "\n")
  cat("    p-value match:", all.equal(res_cc$p.value, wt_cc$p.value, tolerance = 1e-6), "\n")

  # --- CI ---
  res_ci <- ses_calc(x = x, y = y, paired = !is.null(y), ses = "cstat", mu = mu,
                     se_method = "score", correct = FALSE)
  cat("  Score 95% CI (cstat): [", format(res_ci$conf.int[1], digits = 5), ",",
      format(res_ci$conf.int[2], digits = 5), "]\n")
  cat("  Estimate (cstat):", format(res_ci$estimate, digits = 5), "\n")

  # --- All scales ---
  res_rb <- ses_calc(x = x, y = y, paired = !is.null(y), ses = "rb", mu = mu,
                     se_method = "score", correct = FALSE)
  res_odds <- ses_calc(x = x, y = y, paired = !is.null(y), ses = "odds", mu = mu,
                       se_method = "score", correct = FALSE)
  res_lo <- ses_calc(x = x, y = y, paired = !is.null(y), ses = "logodds", mu = mu,
                     se_method = "score", correct = FALSE)

  cat("  rb:      est =", format(res_rb$estimate, digits = 5),
      " CI = [", format(res_rb$conf.int[1], digits = 5), ",",
      format(res_rb$conf.int[2], digits = 5), "]\n")
  cat("  odds:    est =", format(res_odds$estimate, digits = 5),
      " CI = [", format(res_odds$conf.int[1], digits = 5), ",",
      format(res_odds$conf.int[2], digits = 5), "]\n")
  cat("  logodds: est =", format(res_lo$estimate, digits = 5),
      " CI = [", format(res_lo$conf.int[1], digits = 5), ",",
      format(res_lo$conf.int[2], digits = 5), "]\n")

  # --- Verify scale consistency ---
  ci_c <- res_ci$conf.int
  expect_rb <- 2 * ci_c - 1
  expect_odds <- ci_c / (1 - ci_c)
  expect_lo <- log(ci_c / (1 - ci_c))
  cat("  Scale consistency (rb):", all.equal(as.numeric(res_rb$conf.int), as.numeric(expect_rb)), "\n")
  cat("  Scale consistency (odds):", all.equal(as.numeric(res_odds$conf.int), as.numeric(expect_odds)), "\n")
  cat("  Scale consistency (logodds):", all.equal(as.numeric(res_lo$conf.int), as.numeric(expect_lo)), "\n")

  # --- Compare with agresti ---
  res_agr <- suppressMessages(
    ses_calc(x = x, y = y, paired = !is.null(y), ses = "cstat", mu = mu,
             se_method = "agresti")
  )
  cat("  Agresti 95% CI (cstat): [", format(res_agr$conf.int[1], digits = 5), ",",
      format(res_agr$conf.int[2], digits = 5), "]\n")

  # --- Note string ---
  cat("  Note:", res_ci$note, "\n")

  cat("\n")
  invisible(res_ci)
}

# Test 1: sleep data --------
cat("=== Test 1: Sleep data (classic paired, no ties) ===\n\n")
data(sleep)
compare_paired(
  x = sleep$extra[sleep$group == 2],
  y = sleep$extra[sleep$group == 1],
  label = "sleep data"
)

# Test 2: Data with ties --------
cat("=== Test 2: Data with ties ===\n\n")
compare_paired(
  x = c(1, 2, 2, 3, 4, 5, 5, 6),
  y = c(2, 2, 3, 3, 3, 4, 5, 5),
  label = "tied data"
)

# Test 3: One-sample --------
cat("=== Test 3: One-sample ===\n\n")
compare_paired(
  x = c(1.5, 2.3, 3.1, 0.8, 4.2, 2.7, 3.9, 1.1, 5.0, 2.5),
  mu = 0,
  label = "one-sample (all positive)"
)

# Test 4: Heavy ties (ordinal) --------
cat("=== Test 4: Heavy ties (ordinal-like) ===\n\n")
set.seed(123)
compare_paired(
  x = sample(1:5, 20, replace = TRUE),
  y = sample(1:5, 20, replace = TRUE),
  label = "ordinal 1-5 (seed=123)"
)

# Test 5: Complete separation --------
cat("=== Test 5: Complete separation ===\n\n")
compare_paired(
  x = c(10, 11, 12, 13, 14),
  y = c(1, 2, 3, 4, 5),
  label = "all x > y"
)

compare_paired(
  x = c(1, 2, 3, 4, 5),
  y = c(10, 11, 12, 13, 14),
  label = "all y > x"
)

# Test 6: Data with zero differences --------
cat("=== Test 6: Zero differences (should be dropped) ===\n\n")
compare_paired(
  x = c(1, 2, 3, 4, 5, 6, 7),
  y = c(2, 2, 4, 4, 3, 7, 5),
  label = "2 zero diffs out of 7"
)

# Test 7: CI/p-value coherence --------
cat("=== Test 7: CI/p-value coherence for equivalence ===\n\n")
set.seed(123)
x7 <- sample(1:5, 20, replace = TRUE)
y7 <- sample(1:5, 20, replace = TRUE)

res_equiv <- ses_calc(x = x7, y = y7, paired = TRUE, ses = "cstat",
                      se_method = "score", alpha = 0.05,
                      alternative = "equivalence",
                      null.value = c(0.3, 0.7))
ci90 <- res_equiv$conf.int

cat("  90% CI (from equivalence, alpha=0.05):", format(ci90[1], digits = 6), ",",
    format(ci90[2], digits = 6), "\n")
cat("  Equivalence p-value:", format(res_equiv$p.value, digits = 6), "\n\n")

# Test at each CI bound
res_lo <- ses_calc(x = x7, y = y7, paired = TRUE, ses = "cstat",
                   se_method = "score",
                   alternative = "greater", null.value = ci90[1])
res_hi <- ses_calc(x = x7, y = y7, paired = TRUE, ses = "cstat",
                   se_method = "score",
                   alternative = "less", null.value = ci90[2])
cat("  p(cstat > lower CI bound):", format(res_lo$p.value, digits = 6),
    "(expect ~0.05)\n")
cat("  p(cstat < upper CI bound):", format(res_hi$p.value, digits = 6),
    "(expect ~0.05)\n")
cat("  Coherence check:", all.equal(res_lo$p.value, 0.05, tolerance = 1e-4),
    "/", all.equal(res_hi$p.value, 0.05, tolerance = 1e-4), "\n\n")

# Test 8: Equivalence on rb scale --------
cat("=== Test 8: Equivalence test on rb scale ===\n\n")
data(sleep)
x8 <- sleep$extra[sleep$group == 2]
y8 <- sleep$extra[sleep$group == 1]

res_eq_rb <- ses_calc(x = x8, y = y8, paired = TRUE, ses = "rb",
                      se_method = "score",
                      alternative = "equivalence",
                      null.value = c(-0.6, 0.6))
cat("  rb equivalence bounds: [-0.6, 0.6]\n")
cat("  Estimate:", format(res_eq_rb$estimate, digits = 5), "\n")
cat("  90% CI:", format(res_eq_rb$conf.int[1], digits = 5), ",",
    format(res_eq_rb$conf.int[2], digits = 5), "\n")
cat("  p-value:", format(res_eq_rb$p.value, digits = 6), "\n\n")

# Test 9: Larger sample --------
cat("=== Test 9: Larger sample (n=50) ===\n\n")
set.seed(999)
compare_paired(
  x = rnorm(50, mean = 0.3),
  y = rnorm(50, mean = 0),
  label = "n=50, small shift"
)

# Summary --------
cat("=============================================================\n")
cat("  Summary\n")
cat("=============================================================\n")
cat("  - Two-sided p-values match wilcox.test at pi0 = 0.5\n")
cat("  - Continuity correction matches wilcox.test correct = TRUE\n")
cat("  - Wilson-score CI has closed-form (no root-finding)\n")
cat("  - CI handles boundaries naturally (no Haldane correction)\n")
cat("  - CI/p-value coherence guaranteed by test inversion\n")
cat("  - Scale transformations (cstat/rb/odds/logodds) consistent\n")
cat("  - Zero differences dropped (matches wilcox.test behavior)\n")
cat("=============================================================\n")
