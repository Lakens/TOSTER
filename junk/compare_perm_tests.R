# Comparison of permutation t-test implementations
# TOSTER::perm_t_test vs MKinfer::perm.t.test vs coin::oneway_test
#
# Purpose: Demonstrate differences between studentized (TOSTER, coin)
# and non-studentized (MKinfer) permutation t-tests

library(TOSTER)
library(MKinfer)
library(coin)

cat("=======================================================================\n")
cat("PERMUTATION T-TEST COMPARISON: TOSTER vs MKinfer vs coin\n")
cat("=======================================================================\n\n")

# -----------------------------------------------------------------------------
# ONE-SAMPLE TEST
# -----------------------------------------------------------------------------
cat("=== ONE-SAMPLE TEST ===\n\n")

set.seed(12345)
x <- rnorm(8, mean = 0.5, sd = 1)

cat("Data: n =", length(x), "\n")
cat("Mean:", round(mean(x), 4), "\n")
cat("SD:", round(sd(x), 4), "\n\n")

# TOSTER (studentized)
result_toster1 <- perm_t_test(x, mu = 0, alternative = "two.sided", R = 9999)
result_toster2 <- perm_t_test(x, mu = 0, alternative = "two.sided", R = 9999,
                              perm_se = TRUE,
                              p_method = "o")
result_toster3 <- perm_t_test(x, mu = 0, alternative = "two.sided", R = 9999,
                              perm_se = FALSE,
                              p_method = "o")
# MKinfer (non-studentized)
result_mkinfer <- perm.t.test(x, mu = 0, alternative = "two.sided", R = 9999)

# coin doesn't have a direct one-sample permutation test,
# but we can use symmetry_test which is similar
# For one-sample, we test if the distribution is symmetric around mu
#result_coin <- wilcoxsign_test(x ~ 1, distribution = "exact")

cat("One-sample permutation t-test (H0: mu = 0):\n")
cat("-----------------------------------------\n")
cat("TOSTER p-value:  ", round(result_toster1$p.value, 4), " (studentized)\n")
cat("TOSTER p-value:  ", round(result_toster2$p.value, 4), " (studentized, original perm)\n")
cat("TOSTER p-value:  ", round(result_toster3$p.value, 4), " (non-studentized, original perm)\n")
cat("MKinfer p-value: ", round(result_mkinfer$perm.p.value, 4), " (non-studentized)\n")
#cat("coin p-value:    ", round(pvalue(result_coin), 4), " (Wilcoxon signed-rank, exact)\n")
cat("Parametric t p:  ", round(t.test(x, mu = 0)$p.value, 4), "\n\n")

# -----------------------------------------------------------------------------
# TWO-SAMPLE INDEPENDENT TEST
# -----------------------------------------------------------------------------
cat("=== TWO-SAMPLE INDEPENDENT TEST ===\n\n")

set.seed(42)
x <- rnorm(12, mean = 0, sd = 1)
y <- rnorm(10, mean = 0.8, sd = 1.2)

cat("Group 1: n =", length(x), ", mean =", round(mean(x), 4), ", sd =", round(sd(x), 4), "\n")
cat("Group 2: n =", length(y), ", mean =", round(mean(y), 4), ", sd =", round(sd(y), 4), "\n\n")

# TOSTER (studentized)
result_toster1 <- perm_t_test(x, y, alternative = "two.sided", var.equal = FALSE, R = 9999)
result_toster2 <- perm_t_test(x, y, alternative = "two.sided", var.equal = FALSE, R = 9999,
                              perm_se = TRUE,
                              p_method = "o")
result_toster3 <- perm_t_test(x, y, alternative = "two.sided", var.equal = FALSE, R = 9999,
                              perm_se = FALSE,
                              p_method = "o")

# MKinfer (their version)
result_mkinfer <- perm.t.test(x, y, alternative = "two.sided", var.equal = FALSE, R = 9999)

# coin (studentized permutation test)
df <- data.frame(
  value = c(x, y),
  group = factor(c(rep("A", length(x)), rep("B", length(y))))
)
result_coin <- oneway_test(value ~ group, data = df,
                           distribution = approximate(nresample = 9999))

cat("Two-sample permutation t-test (H0: mu_x = mu_y):\n")
cat("------------------------------------------------\n")
cat("TOSTER p-value:  ", round(result_toster1$p.value, 4), " (studentized)\n")
cat("TOSTER p-value:  ", round(result_toster2$p.value, 4), " (studentized, original perm)\n")
cat("TOSTER p-value:  ", round(result_toster3$p.value, 4), " (non-studentized, original perm)\n")
cat("MKinfer p-value: ", round(result_mkinfer$perm.p.value, 4), "\n")
cat("coin p-value:    ", round(pvalue(result_coin), 4), " (asymptotic studentized)\n")
cat("Parametric t p:  ", round(t.test(x, y, var.equal = FALSE)$p.value, 4), "\n\n")

# -----------------------------------------------------------------------------
# TWO-SAMPLE WITH EQUAL VARIANCES
# -----------------------------------------------------------------------------
cat("=== TWO-SAMPLE WITH EQUAL VARIANCES ===\n\n")

# TOSTER (studentized, pooled variance)
result_toster_eq1 <- perm_t_test(x, y, alternative = "two.sided", var.equal = TRUE, R = 9999)
result_toster_eq2 <- perm_t_test(x, y, alternative = "two.sided", var.equal = TRUE, R = 9999,
                                perm_se = TRUE,
                                p_method = "o")
result_toster_eq3 <- perm_t_test(x, y, alternative = "two.sided", var.equal = TRUE, R = 9999,
                                perm_se = FALSE,
                                p_method = "o")
# MKinfer (pooled variance)
result_mkinfer_eq <- perm.t.test(x, y, alternative = "two.sided", var.equal = TRUE, R = 9999)

# coin with pooled variance assumption
result_coin_eq <- oneway_test(value ~ group, data = df,
                              distribution = approximate(nresample = 9999))

cat("Two-sample permutation t-test (pooled variance):\n")
cat("------------------------------------------------\n")
cat("TOSTER p-value:  ", round(result_toster_eq1$p.value, 4), " (studentized)\n")
cat("TOSTER p-value:  ", round(result_toster_eq2$p.value, 4), " (studentized, original perm)\n")
cat("TOSTER p-value:  ", round(result_toster_eq3$p.value, 4), " (non-studentized, original perm)\n")
cat("MKinfer p-value: ", round(result_mkinfer_eq$perm.p.value, 4), "\n")
cat("coin p-value:    ", round(pvalue(result_coin_eq), 4), "\n")
cat("Parametric t p:  ", round(t.test(x, y, var.equal = TRUE)$p.value, 4), "\n\n")

# -----------------------------------------------------------------------------
# PAIRED SAMPLES TEST
# -----------------------------------------------------------------------------
cat("=== PAIRED SAMPLES TEST ===\n\n")

set.seed(123)
before <- rnorm(10, mean = 100, sd = 15)
after <- before + rnorm(10, mean = 5, sd = 8)  # Small improvement

cat("Paired data: n =", length(before), "pairs\n")
cat("Mean difference:", round(mean(after - before), 4), "\n")
cat("SD of differences:", round(sd(after - before), 4), "\n\n")

# TOSTER (studentized)
result_toster_paired <- perm_t_test(before, after, paired = TRUE,
                                     alternative = "two.sided", R = 9999)

# MKinfer
result_mkinfer_paired <- perm.t.test(before, after, paired = TRUE,
                                      alternative = "two.sided", R = 9999)

# coin (using difference scores)
#diff_scores <- after - before
#result_coin_paired <- wilcoxsign_test(after ~ before,
#                                       distribution = approximate(nresample = 9999))

cat("Paired samples permutation t-test:\n")
cat("----------------------------------\n")
cat("TOSTER p-value:  ", round(result_toster_paired$p.value, 4), " (studentized)\n")
cat("MKinfer p-value: ", round(result_mkinfer_paired$perm.p.value, 4), "\n")
#cat("coin p-value:    ", round(pvalue(result_coin_paired), 4), " (Wilcoxon signed-rank)\n")
cat("Parametric t p:  ", round(t.test(before, after, paired = TRUE)$p.value, 4), "\n\n")

# -----------------------------------------------------------------------------
# EXACT PERMUTATION (SMALL SAMPLE)
# -----------------------------------------------------------------------------
cat("=== EXACT PERMUTATION TEST (SMALL SAMPLE) ===\n\n")

set.seed(999)
x_small <- rnorm(6, mean = 0, sd = 1)
y_small <- rnorm(5, mean = 1, sd = 1)

cat("Group 1: n =", length(x_small), "\n")
cat("Group 2: n =", length(y_small), "\n")
cat("Total permutations possible:", choose(11, 6), "\n\n")

# TOSTER (exact)
result_toster_exact <- perm_t_test(x_small, y_small, alternative = "two.sided",
                                    var.equal = TRUE, R = 9999)

# MKinfer (exact)
result_mkinfer_exact <- perm.t.test(x_small, y_small, alternative = "two.sided",
                                     var.equal = TRUE, R = 9999)

# coin (exact)
df_small <- data.frame(
  value = c(x_small, y_small),
  group = factor(c(rep("A", length(x_small)), rep("B", length(y_small))))
)
result_coin_exact <- oneway_test(value ~ group, data = df_small,
                                  distribution = "exact")

cat("Exact permutation test:\n")
cat("-----------------------\n")
cat("TOSTER p-value:  ", round(result_toster_exact$p.value, 4), "\n")
cat("MKinfer p-value: ", round(result_mkinfer_exact$perm.p.value, 4), "\n")
cat("coin p-value:    ", round(pvalue(result_coin_exact), 4), " (exact)\n")
cat("Parametric t p:  ", round(t.test(x_small, y_small, var.equal = TRUE)$p.value, 4), "\n\n")

# -----------------------------------------------------------------------------
# SUMMARY
# -----------------------------------------------------------------------------
cat("=======================================================================\n")
cat("SUMMARY\n")
cat("=======================================================================\n\n")

cat("Key differences between implementations:\n\n")

cat("1. TOSTER::perm_t_test\n")
cat("   - Uses STUDENTIZED permutation test (Janssen 1997, Chung & Romano 2013)\n")
cat("   - Recomputes both mean AND standard error for each permutation\n")
cat("   - More robust to heteroscedasticity\n")
cat("   - Supports trimmed means\n\n")

cat("2. MKinfer::perm.t.test\n")
cat("   - Uses NON-STUDENTIZED approach for one-sample test\n")
cat("   - Uses original standard error, only permutes the mean\n")
cat("   - Follows Efron & Tibshirani (1993) Chapter 15 approach\n\n")

cat("3. coin::oneway_test\n")
cat("   - Uses studentized permutation test\n")
cat("   - Highly flexible framework\n")
cat("   - Results should be similar to TOSTER for two-sample tests\n\n")

cat("For two-sample tests, TOSTER and coin should give similar results.\n")
cat("For one-sample tests, TOSTER (studentized) and MKinfer (non-studentized)\n")
cat("will differ systematically due to their different approaches.\n")
