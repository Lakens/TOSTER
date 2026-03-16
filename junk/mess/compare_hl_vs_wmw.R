# Comparison: hodges_lehmann() vs simple_htest(test = "wilcox.test")
# for equivalence testing
#
# Purpose: Both functions perform nonparametric location tests related to
# the Wilcoxon-Mann-Whitney framework. This script tests whether they are
# assessing the same underlying hypothesis but with different methods.
#
# Key relationships:
#   - wilcox.test estimates the Hodges-Lehmann estimator as its location shift
#   - hodges_lehmann uses permutation/asymptotic inference on the HL estimator
#   - simple_htest wraps wilcox.test with TOST extensions
#   - Under H0, both should target the same location parameter
#
# Questions to answer:
#   1. Do the point estimates (HL estimator) agree?
#   2. Do the confidence intervals agree?
#   3. Do the p-values agree for standard alternatives?
#   4. Do the TOST (equivalence) results agree?
#   5. Where and why do they diverge?

library(TOSTER)

cat("=======================================================================\n")
cat("COMPARISON: hodges_lehmann() vs simple_htest(wilcox.test)\n")
cat("=======================================================================\n\n")

# =============================================================================
# SECTION 1: Two-sample independent test
# =============================================================================
cat("=== SECTION 1: TWO-SAMPLE INDEPENDENT TEST ===\n\n")

set.seed(2024)
x <- rnorm(20, mean = 0, sd = 1)
y <- rnorm(20, mean = 0.5, sd = 1.2)

cat("Group x: n =", length(x), ", mean =", round(mean(x), 4),
    ", median =", round(median(x), 4), "\n")
cat("Group y: n =", length(y), ", mean =", round(mean(y), 4),
    ", median =", round(median(y), 4), "\n\n")

# --- Two-sided test ---
cat("--- Two-sided test (H0: shift = 0) ---\n\n")

hl_2s <- hodges_lehmann(x, y, alternative = "two.sided", mu = 0)
sh_2s <- simple_htest(x, y, test = "wilcox.test", alternative = "two.sided", mu = 0)
wt_2s <- wilcox.test(x, y, alternative = "two.sided", mu = 0, conf.int = TRUE)

cat("Point estimates (Hodges-Lehmann):\n")
cat("  hodges_lehmann:  ", round(hl_2s$estimate, 6), "\n")
cat("  simple_htest:    ", round(sh_2s$estimate, 6), "\n")
cat("  wilcox.test:     ", round(wt_2s$estimate, 6), "\n\n")

cat("Confidence intervals:\n")
cat("  hodges_lehmann:  [", round(hl_2s$conf.int[1], 6), ",",
    round(hl_2s$conf.int[2], 6), "] (asymptotic)\n")
cat("  simple_htest:    [", round(sh_2s$conf.int[1], 6), ",",
    round(sh_2s$conf.int[2], 6), "]\n")
cat("  wilcox.test:     [", round(wt_2s$conf.int[1], 6), ",",
    round(wt_2s$conf.int[2], 6), "]\n\n")

cat("P-values:\n")
cat("  hodges_lehmann:  ", round(hl_2s$p.value, 6), "(asymptotic)\n")
cat("  simple_htest:    ", round(sh_2s$p.value, 6), "\n")
cat("  wilcox.test:     ", round(wt_2s$p.value, 6), "\n\n")

# --- With permutation ---
cat("--- Two-sided test (with permutation, R = 5000) ---\n\n")

set.seed(42)
hl_2s_perm <- hodges_lehmann(x, y, alternative = "two.sided", mu = 0, R = 5000)

cat("P-values:\n")
cat("  hodges_lehmann (perm): ", round(hl_2s_perm$p.value, 6), "\n")
cat("  hodges_lehmann (asym): ", round(hl_2s$p.value, 6), "\n")
cat("  simple_htest (wilcox): ", round(sh_2s$p.value, 6), "\n\n")

cat("Confidence intervals:\n")
cat("  hodges_lehmann (perm): [", round(hl_2s_perm$conf.int[1], 6), ",",
    round(hl_2s_perm$conf.int[2], 6), "]\n")
cat("  hodges_lehmann (asym): [", round(hl_2s$conf.int[1], 6), ",",
    round(hl_2s$conf.int[2], 6), "]\n")
cat("  simple_htest (wilcox): [", round(sh_2s$conf.int[1], 6), ",",
    round(sh_2s$conf.int[2], 6), "]\n\n")

# =============================================================================
# SECTION 2: Equivalence test (the main comparison)
# =============================================================================
cat("=== SECTION 2: EQUIVALENCE TESTING (TOST) ===\n\n")

eqbound <- 1  # equivalence bound of +/- 1

cat("Equivalence bounds: [", -eqbound, ",", eqbound, "]\n\n")

# --- Equivalence ---
hl_eq <- hodges_lehmann(x, y, alternative = "equivalence", mu = eqbound)
sh_eq <- simple_htest(x, y, test = "wilcox.test",
                      alternative = "equivalence", mu = eqbound)

cat("--- Equivalence test results ---\n\n")

cat("Point estimates:\n")
cat("  hodges_lehmann:  ", round(hl_eq$estimate, 6), "\n")
cat("  simple_htest:    ", round(sh_eq$estimate, 6), "\n\n")

cat("Null values:\n")
cat("  hodges_lehmann:  ", hl_eq$null.value, "\n")
cat("  simple_htest:    ", sh_eq$null.value, "\n\n")

cat("TOST p-values:\n")
cat("  hodges_lehmann:  ", round(hl_eq$p.value, 6), "\n")
cat("  simple_htest:    ", round(sh_eq$p.value, 6), "\n\n")

cat("Confidence intervals (1-2*alpha = 90%):\n")
cat("  hodges_lehmann:  [", round(hl_eq$conf.int[1], 6), ",",
    round(hl_eq$conf.int[2], 6), "]\n")
cat("  simple_htest:    [", round(sh_eq$conf.int[1], 6), ",",
    round(sh_eq$conf.int[2], 6), "]\n\n")

# --- With permutation ---
cat("--- Equivalence test with permutation (R = 5000) ---\n\n")

set.seed(42)
hl_eq_perm <- hodges_lehmann(x, y, alternative = "equivalence",
                              mu = eqbound, R = 5000)

cat("TOST p-values:\n")
cat("  hodges_lehmann (perm): ", round(hl_eq_perm$p.value, 6), "\n")
cat("  hodges_lehmann (asym): ", round(hl_eq$p.value, 6), "\n")
cat("  simple_htest (wilcox): ", round(sh_eq$p.value, 6), "\n\n")

cat("Confidence intervals (90%):\n")
cat("  hodges_lehmann (perm): [", round(hl_eq_perm$conf.int[1], 6), ",",
    round(hl_eq_perm$conf.int[2], 6), "]\n")
cat("  hodges_lehmann (asym): [", round(hl_eq$conf.int[1], 6), ",",
    round(hl_eq$conf.int[2], 6), "]\n")
cat("  simple_htest (wilcox): [", round(sh_eq$conf.int[1], 6), ",",
    round(sh_eq$conf.int[2], 6), "]\n\n")

# =============================================================================
# SECTION 3: Minimal effect test
# =============================================================================
cat("=== SECTION 3: MINIMAL EFFECT TESTING ===\n\n")

hl_met <- hodges_lehmann(x, y, alternative = "minimal.effect", mu = eqbound)
sh_met <- simple_htest(x, y, test = "wilcox.test",
                       alternative = "minimal.effect", mu = eqbound)

cat("MET p-values:\n")
cat("  hodges_lehmann:  ", round(hl_met$p.value, 6), "\n")
cat("  simple_htest:    ", round(sh_met$p.value, 6), "\n\n")

cat("Confidence intervals (90%):\n")
cat("  hodges_lehmann:  [", round(hl_met$conf.int[1], 6), ",",
    round(hl_met$conf.int[2], 6), "]\n")
cat("  simple_htest:    [", round(sh_met$conf.int[1], 6), ",",
    round(sh_met$conf.int[2], 6), "]\n\n")

# =============================================================================
# SECTION 4: One-sided tests
# =============================================================================
cat("=== SECTION 4: ONE-SIDED TESTS ===\n\n")

for (alt in c("less", "greater")) {
  hl_os <- hodges_lehmann(x, y, alternative = alt, mu = 0)
  sh_os <- simple_htest(x, y, test = "wilcox.test", alternative = alt, mu = 0)

  cat("Alternative:", alt, "\n")
  cat("  hodges_lehmann p-value: ", round(hl_os$p.value, 6), "\n")
  cat("  simple_htest p-value:   ", round(sh_os$p.value, 6), "\n")
  cat("  hodges_lehmann CI: [", round(hl_os$conf.int[1], 6), ",",
      round(hl_os$conf.int[2], 6), "]\n")
  cat("  simple_htest CI:   [", round(sh_os$conf.int[1], 6), ",",
      round(sh_os$conf.int[2], 6), "]\n\n")
}

# =============================================================================
# SECTION 5: Paired samples
# =============================================================================
cat("=== SECTION 5: PAIRED SAMPLES ===\n\n")

set.seed(2024)
before <- rnorm(15, mean = 50, sd = 10)
after <- before + rnorm(15, mean = 2, sd = 5)

cat("Paired data: n =", length(before), "pairs\n")
cat("Mean difference:", round(mean(after - before), 4), "\n\n")

# --- Two-sided ---
hl_p <- hodges_lehmann(before, after, paired = TRUE,
                       alternative = "two.sided", mu = 0)
sh_p <- simple_htest(before, after, paired = TRUE,
                     test = "wilcox.test", alternative = "two.sided", mu = 0)
wt_p <- wilcox.test(before, after, paired = TRUE,
                    alternative = "two.sided", mu = 0, conf.int = TRUE)

cat("--- Paired two-sided ---\n")
cat("Point estimates:\n")
cat("  hodges_lehmann:  ", round(hl_p$estimate, 6), "\n")
cat("  simple_htest:    ", round(sh_p$estimate, 6), "\n")
cat("  wilcox.test:     ", round(wt_p$estimate, 6), "\n\n")

cat("P-values:\n")
cat("  hodges_lehmann:  ", round(hl_p$p.value, 6), "\n")
cat("  simple_htest:    ", round(sh_p$p.value, 6), "\n")
cat("  wilcox.test:     ", round(wt_p$p.value, 6), "\n\n")

cat("Confidence intervals:\n")
cat("  hodges_lehmann:  [", round(hl_p$conf.int[1], 6), ",",
    round(hl_p$conf.int[2], 6), "] (asymptotic)\n")
cat("  simple_htest:    [", round(sh_p$conf.int[1], 6), ",",
    round(sh_p$conf.int[2], 6), "]\n")
cat("  wilcox.test:     [", round(wt_p$conf.int[1], 6), ",",
    round(wt_p$conf.int[2], 6), "]\n\n")

# --- Paired equivalence ---
cat("--- Paired equivalence (bounds = +/- 5) ---\n")
eqbound_paired <- 5

hl_peq <- hodges_lehmann(before, after, paired = TRUE,
                          alternative = "equivalence", mu = eqbound_paired)
sh_peq <- simple_htest(before, after, paired = TRUE,
                       test = "wilcox.test",
                       alternative = "equivalence", mu = eqbound_paired)

cat("TOST p-values:\n")
cat("  hodges_lehmann:  ", round(hl_peq$p.value, 6), "\n")
cat("  simple_htest:    ", round(sh_peq$p.value, 6), "\n\n")

cat("Confidence intervals (90%):\n")
cat("  hodges_lehmann:  [", round(hl_peq$conf.int[1], 6), ",",
    round(hl_peq$conf.int[2], 6), "]\n")
cat("  simple_htest:    [", round(sh_peq$conf.int[1], 6), ",",
    round(sh_peq$conf.int[2], 6), "]\n\n")

# =============================================================================
# SECTION 6: One-sample test
# =============================================================================
cat("=== SECTION 6: ONE-SAMPLE TEST ===\n\n")

set.seed(2024)
z <- rnorm(25, mean = 0.3, sd = 1)

cat("Data: n =", length(z), ", mean =", round(mean(z), 4),
    ", median =", round(median(z), 4), "\n\n")

# --- Two-sided ---
hl_1s <- hodges_lehmann(z, alternative = "two.sided", mu = 0)
sh_1s <- simple_htest(z, test = "wilcox.test", alternative = "two.sided", mu = 0)
wt_1s <- wilcox.test(z, alternative = "two.sided", mu = 0, conf.int = TRUE)

cat("--- One-sample two-sided ---\n")
cat("Point estimates (pseudomedian):\n")
cat("  hodges_lehmann:  ", round(hl_1s$estimate, 6), "\n")
cat("  simple_htest:    ", round(sh_1s$estimate, 6), "\n")
cat("  wilcox.test:     ", round(wt_1s$estimate, 6), "\n\n")

cat("P-values:\n")
cat("  hodges_lehmann:  ", round(hl_1s$p.value, 6), "\n")
cat("  simple_htest:    ", round(sh_1s$p.value, 6), "\n")
cat("  wilcox.test:     ", round(wt_1s$p.value, 6), "\n\n")

cat("Confidence intervals:\n")
cat("  hodges_lehmann:  [", round(hl_1s$conf.int[1], 6), ",",
    round(hl_1s$conf.int[2], 6), "] (asymptotic)\n")
cat("  simple_htest:    [", round(sh_1s$conf.int[1], 6), ",",
    round(sh_1s$conf.int[2], 6), "]\n")
cat("  wilcox.test:     [", round(wt_1s$conf.int[1], 6), ",",
    round(wt_1s$conf.int[2], 6), "]\n\n")

# --- One-sample equivalence ---
cat("--- One-sample equivalence (bounds = +/- 0.5) ---\n")
eqbound_1s <- 0.5

hl_1eq <- hodges_lehmann(z, alternative = "equivalence", mu = eqbound_1s)
sh_1eq <- simple_htest(z, test = "wilcox.test",
                       alternative = "equivalence", mu = eqbound_1s)

cat("TOST p-values:\n")
cat("  hodges_lehmann:  ", round(hl_1eq$p.value, 6), "\n")
cat("  simple_htest:    ", round(sh_1eq$p.value, 6), "\n\n")

cat("Confidence intervals (90%):\n")
cat("  hodges_lehmann:  [", round(hl_1eq$conf.int[1], 6), ",",
    round(hl_1eq$conf.int[2], 6), "]\n")
cat("  simple_htest:    [", round(sh_1eq$conf.int[1], 6), ",",
    round(sh_1eq$conf.int[2], 6), "]\n\n")

# =============================================================================
# SECTION 7: Systematic comparison across sample sizes
# =============================================================================
cat("=== SECTION 7: SYSTEMATIC COMPARISON (varying n) ===\n\n")

cat("Comparing equivalence test p-values and CIs across sample sizes\n")
cat("Equivalence bounds: +/- 0.8\n\n")

ns <- c(10, 20, 30, 50, 100)
eqb_sys <- 0.8

results <- data.frame(
  n = integer(),
  hl_asym_p = numeric(),
  sh_wilcox_p = numeric(),
  p_diff = numeric(),
  hl_estimate = numeric(),
  sh_estimate = numeric(),
  est_diff = numeric(),
  stringsAsFactors = FALSE
)

set.seed(2024)
for (n in ns) {
  x_sim <- rnorm(n, mean = 0, sd = 1)
  y_sim <- rnorm(n, mean = 0.3, sd = 1)

  hl_sim <- hodges_lehmann(x_sim, y_sim,
                            alternative = "equivalence", mu = eqb_sys)
  sh_sim <- simple_htest(x_sim, y_sim, test = "wilcox.test",
                         alternative = "equivalence", mu = eqb_sys)

  results <- rbind(results, data.frame(
    n = n,
    hl_asym_p = hl_sim$p.value,
    sh_wilcox_p = sh_sim$p.value,
    p_diff = abs(hl_sim$p.value - sh_sim$p.value),
    hl_estimate = hl_sim$estimate,
    sh_estimate = sh_sim$estimate,
    est_diff = abs(hl_sim$estimate - sh_sim$estimate)
  ))
}

cat("  n   | HL(asym) p | WMW p      | |p diff|  | HL est  | WMW est | |est diff|\n")
cat("------+------------+------------+-----------+---------+---------+----------\n")
for (i in seq_len(nrow(results))) {
  cat(sprintf("  %-4d| %-10.6f | %-10.6f | %-9.6f | %-7.4f | %-7.4f | %-8.6f\n",
              results$n[i],
              results$hl_asym_p[i],
              results$sh_wilcox_p[i],
              results$p_diff[i],
              results$hl_estimate[i],
              results$sh_estimate[i],
              results$est_diff[i]))
}

cat("\n")

# =============================================================================
# SECTION 8: Test under null (both samples from same distribution)
# =============================================================================
cat("=== SECTION 8: UNDER THE NULL (same distribution) ===\n\n")

set.seed(999)
x_null <- rnorm(30)
y_null <- rnorm(30)

eqb_null <- 0.5

hl_null_eq <- hodges_lehmann(x_null, y_null,
                              alternative = "equivalence", mu = eqb_null)
sh_null_eq <- simple_htest(x_null, y_null, test = "wilcox.test",
                           alternative = "equivalence", mu = eqb_null)

cat("Both groups: n = 30, drawn from N(0,1)\n")
cat("Equivalence bounds: +/-", eqb_null, "\n\n")

cat("Equivalence test:\n")
cat("  hodges_lehmann p: ", round(hl_null_eq$p.value, 6), "\n")
cat("  simple_htest p:   ", round(sh_null_eq$p.value, 6), "\n\n")

cat("Two-sided test:\n")
hl_null_2s <- hodges_lehmann(x_null, y_null, alternative = "two.sided")
sh_null_2s <- simple_htest(x_null, y_null, test = "wilcox.test",
                           alternative = "two.sided")
cat("  hodges_lehmann p: ", round(hl_null_2s$p.value, 6), "\n")
cat("  simple_htest p:   ", round(sh_null_2s$p.value, 6), "\n\n")

# =============================================================================
# SECTION 9: Non-normal data (where robustness matters)
# =============================================================================
cat("=== SECTION 9: NON-NORMAL DATA (contaminated normal) ===\n\n")

set.seed(2024)
# Contaminated normal: 80% N(0,1) + 20% N(5,1) -- outlier group
x_contam <- c(rnorm(24, 0, 1), rnorm(6, 5, 1))
y_contam <- rnorm(30, 0, 1)

cat("Group x: 80% N(0,1) + 20% N(5,1) contamination\n")
cat("Group y: N(0,1)\n")
cat("x mean =", round(mean(x_contam), 4),
    ", x median =", round(median(x_contam), 4), "\n")
cat("y mean =", round(mean(y_contam), 4),
    ", y median =", round(median(y_contam), 4), "\n\n")

eqb_contam <- 1.5

cat("Two-sided test:\n")
hl_c2s <- hodges_lehmann(x_contam, y_contam, alternative = "two.sided")
sh_c2s <- simple_htest(x_contam, y_contam, test = "wilcox.test",
                       alternative = "two.sided")
cat("  hodges_lehmann:  p =", round(hl_c2s$p.value, 6),
    ", est =", round(hl_c2s$estimate, 4), "\n")
cat("  simple_htest:    p =", round(sh_c2s$p.value, 6),
    ", est =", round(sh_c2s$estimate, 4), "\n\n")

cat("Equivalence test (bounds = +/-", eqb_contam, "):\n")
hl_ceq <- hodges_lehmann(x_contam, y_contam,
                          alternative = "equivalence", mu = eqb_contam)
sh_ceq <- simple_htest(x_contam, y_contam, test = "wilcox.test",
                       alternative = "equivalence", mu = eqb_contam)
cat("  hodges_lehmann:  p =", round(hl_ceq$p.value, 6), "\n")
cat("  simple_htest:    p =", round(sh_ceq$p.value, 6), "\n\n")

cat("CI comparison (90%):\n")
cat("  hodges_lehmann:  [", round(hl_ceq$conf.int[1], 4), ",",
    round(hl_ceq$conf.int[2], 4), "]\n")
cat("  simple_htest:    [", round(sh_ceq$conf.int[1], 4), ",",
    round(sh_ceq$conf.int[2], 4), "]\n\n")

# =============================================================================
# SECTION 10: Heavy-tailed data (Cauchy)
# =============================================================================
cat("=== SECTION 10: HEAVY-TAILED DATA (Cauchy) ===\n\n")

set.seed(2024)
x_cauchy <- rcauchy(30, location = 0, scale = 1)
y_cauchy <- rcauchy(30, location = 0.5, scale = 1)

cat("Cauchy(0,1) vs Cauchy(0.5,1), n = 30 each\n\n")

hl_cauchy <- hodges_lehmann(x_cauchy, y_cauchy, alternative = "two.sided")
sh_cauchy <- simple_htest(x_cauchy, y_cauchy, test = "wilcox.test",
                          alternative = "two.sided")

cat("Two-sided test:\n")
cat("  hodges_lehmann:  p =", round(hl_cauchy$p.value, 6),
    ", est =", round(hl_cauchy$estimate, 4), "\n")
cat("  simple_htest:    p =", round(sh_cauchy$p.value, 6),
    ", est =", round(sh_cauchy$estimate, 4), "\n\n")

eqb_cauchy <- 2
hl_cauchy_eq <- hodges_lehmann(x_cauchy, y_cauchy,
                                alternative = "equivalence", mu = eqb_cauchy)
sh_cauchy_eq <- simple_htest(x_cauchy, y_cauchy, test = "wilcox.test",
                             alternative = "equivalence", mu = eqb_cauchy)

cat("Equivalence test (bounds = +/-", eqb_cauchy, "):\n")
cat("  hodges_lehmann:  p =", round(hl_cauchy_eq$p.value, 6), "\n")
cat("  simple_htest:    p =", round(sh_cauchy_eq$p.value, 6), "\n\n")

# =============================================================================
# SUMMARY
# =============================================================================
cat("=======================================================================\n")
cat("SUMMARY\n")
cat("=======================================================================\n\n")

cat("WHAT'S THE SAME:\n")
cat("  - Both target the Hodges-Lehmann estimator (location shift)\n")
cat("  - Point estimates should be identical (both compute the HL estimator)\n")
cat("  - Both implement TOST via two one-sided tests\n")
cat("  - Both use 1-2*alpha CIs for equivalence/MET\n")
cat("  - Under the same distributional assumptions, they test the same null\n\n")

cat("WHAT'S DIFFERENT:\n")
cat("  - Inference method:\n")
cat("      hodges_lehmann: asymptotic (kernel density) or permutation-based\n")
cat("      simple_htest:   wraps wilcox.test (rank-based exact/asymptotic)\n")
cat("  - Confidence intervals:\n")
cat("      hodges_lehmann (asymptotic): based on kernel density variance est.\n")
cat("      wilcox.test: based on inverting the rank test\n")
cat("      These will generally DIFFER even though the point estimates match\n")
cat("  - P-values:\n")
cat("      hodges_lehmann: Z-test from HL estimator / permutation\n")
cat("      wilcox.test: rank-sum (exact or normal approx with continuity)\n")
cat("      These test the same hypothesis but via different statistics\n")
cat("  - Scale estimation (permutation mode):\n")
cat("      hodges_lehmann: robust scale (S1 or S2, Fried & Dehling 2011)\n")
cat("      wilcox.test: no explicit scale estimator (rank-based)\n\n")

cat("CONCLUSION:\n")
cat("  Both functions test the SAME hypothesis about location shift using\n")
cat("  the Hodges-Lehmann estimator, but they use DIFFERENT inference\n")
cat("  procedures. This means:\n")
cat("    - Point estimates will agree\n")
cat("    - P-values and CIs may differ (sometimes substantially)\n")
cat("    - They are complementary rather than redundant\n")
cat("    - hodges_lehmann offers more flexibility (permutation, scale choice)\n")
cat("    - simple_htest(wilcox) leverages well-established rank test theory\n")
