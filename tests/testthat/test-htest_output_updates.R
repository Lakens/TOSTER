# Tests for htest output updates:
# - Estimate labeling and mean difference appending
# - Sample size in output
# - Formula interface group name substitution
# - quote_if_numeric() quoting of numeric factor levels

hush = function(code) {
  sink(nullfile())
  on.exit(sink())
  suppressMessages(code)
}

# ===========================================================================
# simple_htest tests
# ===========================================================================

test_that("simple_htest t-test: two-sample estimate has 3 elements with correct structure", {
  res <- hush(simple_htest(1:10, y = c(7:20), mu = 0))

  # Should have 3 estimates: mean of group x, mean of group y, mean difference
  expect_equal(length(res$estimate), 3)
  expect_equal(names(res$estimate)[1], "mean of group x")
  expect_equal(names(res$estimate)[2], "mean of group y")
  expect_true(grepl("mean difference", names(res$estimate)[3]))
  expect_true(grepl("x - y", names(res$estimate)[3]))

  # Third element should equal first minus second
  expect_equal(unname(res$estimate[3]),
               unname(res$estimate[1] - res$estimate[2]))

  # Backwards compatibility: positional indexing still works
  expect_equal(unname(res$estimate[1]), mean(1:10))
  expect_equal(unname(res$estimate[2]), mean(7:20))

  # Named indexing with new labels
  expect_equal(unname(res$estimate["mean of group x"]), mean(1:10))
  expect_equal(unname(res$estimate["mean of group y"]), mean(7:20))
})

test_that("simple_htest t-test: paired label", {
  data(sleep)
  x_sleep <- sleep$extra[sleep$group == 1]
  y_sleep <- sleep$extra[sleep$group == 2]
  res <- hush(simple_htest(x_sleep, y_sleep, paired = TRUE, mu = 0))

  expect_equal(length(res$estimate), 1)
  expect_equal(names(res$estimate), "mean of the differences (z = x - y)")
})

test_that("simple_htest t-test: one-sample label unchanged", {
  res <- hush(simple_htest(1:20, mu = 0))
  expect_equal(length(res$estimate), 1)
  expect_equal(names(res$estimate), "mean of x")
})

test_that("simple_htest wilcox: estimate labels", {
  # One-sample
  res_one <- hush(simple_htest(1:20, test = "w", mu = 0))
  expect_equal(names(res_one$estimate), "(pseudo)median of x")

  # Paired
  data(sleep)
  x_sleep <- sleep$extra[sleep$group == 1]
  y_sleep <- sleep$extra[sleep$group == 2]
  res_paired <- hush(simple_htest(x_sleep, y_sleep, test = "w",
                                  paired = TRUE, mu = 0))
  expect_equal(names(res_paired$estimate), "Hodges-Lehmann estimate (z = x - y)")

  # Two-sample
  res_two <- hush(simple_htest(1:10, y = c(7:20), test = "w", mu = 0))
  expect_equal(names(res_two$estimate), "Hodges-Lehmann estimate (x - y)")
})

test_that("simple_htest: formula interface substitutes group names with quoting", {
  data(sleep)
  res <- hush(simple_htest(extra ~ group, data = sleep, mu = 0))
  nms <- names(res$estimate)
  # sleep$group levels are "1" and "2" -> should be quoted as '1' and '2'
  expect_equal(nms[1], "mean of group '1'")
  expect_equal(nms[2], "mean of group '2'")
  expect_true(grepl("'1' - '2'", nms[3]))

  # Wilcoxon formula
  res_w <- hush(simple_htest(extra ~ group, data = sleep, test = "w", mu = 0))
  expect_true(grepl("'1' - '2'", names(res_w$estimate)))
})

test_that("simple_htest: formula with non-numeric factor levels", {
  # Create data with text factor levels
  df <- data.frame(
    value = c(rnorm(10, 20), rnorm(10, 25)),
    group = factor(rep(c("auto", "manual"), each = 10))
  )
  res <- hush(simple_htest(value ~ group, data = df, mu = 0))
  nms <- names(res$estimate)
  # Text levels should NOT be quoted
  expect_equal(nms[1], "mean of group auto")
  expect_equal(nms[2], "mean of group manual")
  expect_true(grepl("auto - manual", nms[3]))
})

test_that("simple_htest: sample_size correctness", {
  # One-sample
  res_one <- hush(simple_htest(1:20, mu = 0))
  expect_equal(res_one$sample_size, c(n = 20L))

  # Two-sample
  res_two <- hush(simple_htest(1:10, y = c(7:20), mu = 0))
  expect_equal(res_two$sample_size, c(nx = 10L, ny = 14L))

  # Paired
  data(sleep)
  x_sleep <- sleep$extra[sleep$group == 1]
  y_sleep <- sleep$extra[sleep$group == 2]
  res_paired <- hush(simple_htest(x_sleep, y_sleep, paired = TRUE, mu = 0))
  expect_equal(res_paired$sample_size, c(n = 10L))

  # Paired with NAs
  x_na <- c(1, 2, NA, 4, 5, 3, 7)
  y_na <- c(6, NA, 8, 9, 10, 2, 5)
  res_na <- hush(simple_htest(x_na, y_na, paired = TRUE, mu = 0))
  expect_equal(res_na$sample_size, c(n = 5L))  # 5 complete pairs (indices 1,4,5,6,7)

  # Formula interface: names should be group levels
  res_formula <- hush(simple_htest(extra ~ group, data = sleep, mu = 0))
  expect_equal(names(res_formula$sample_size), c("1", "2"))
  expect_equal(unname(res_formula$sample_size), c(10L, 10L))
})

test_that("simple_htest: equivalence/MET paths work with new estimate structure", {
  # Equivalence t-test two-sample
  res_eq <- hush(simple_htest(1:10, y = c(7:20), alternative = "e", mu = 3))
  expect_equal(length(res_eq$estimate), 3)
  expect_true(grepl("mean difference", names(res_eq$estimate)[3]))
  expect_equal(unname(res_eq$estimate[3]),
               unname(res_eq$estimate[1] - res_eq$estimate[2]))

  # Equivalence wilcox two-sample
  res_eq_w <- hush(simple_htest(1:10, y = c(7:20), test = "w",
                                alternative = "e", mu = 3))
  expect_equal(names(res_eq_w$estimate), "Hodges-Lehmann estimate (x - y)")

  # MET t-test paired
  data(sleep)
  x_sleep <- sleep$extra[sleep$group == 1]
  y_sleep <- sleep$extra[sleep$group == 2]
  res_met <- hush(simple_htest(x_sleep, y_sleep, paired = TRUE,
                               alternative = "m", mu = 2))
  expect_equal(names(res_met$estimate), "mean of the differences (z = x - y)")
})

# ===========================================================================
# boot_t_test tests
# ===========================================================================

test_that("boot_t_test with tr == 0 inherits from simple_htest correctly", {
  set.seed(123)
  res <- hush(boot_t_test(1:10, y = c(7:20), mu = 0, R = 99))

  # Should inherit 3-element estimate from simple_htest
  expect_equal(length(res$estimate), 3)
  expect_true(grepl("mean difference", names(res$estimate)[3]))
  expect_equal(unname(res$estimate[3]),
               unname(res$estimate[1] - res$estimate[2]))

  # sample_size should be inherited
  expect_equal(res$sample_size, c(nx = 10L, ny = 14L))
})

test_that("boot_t_test with tr > 0: two-sample estimate has 3 elements", {
  set.seed(123)
  tr_val <- 0.1
  res <- hush(boot_t_test(1:10, y = c(7:20), mu = 0, tr = tr_val, R = 99))

  expect_equal(length(res$estimate), 3)
  expect_equal(names(res$estimate)[1], "trimmed mean of x")
  expect_equal(names(res$estimate)[2], "trimmed mean of y")
  expect_true(grepl("trimmed mean difference", names(res$estimate)[3]))
  expect_true(grepl(as.character(tr_val), names(res$estimate)[3]))

  # Third element should equal first minus second
  expect_equal(unname(res$estimate[3]),
               unname(res$estimate[1] - res$estimate[2]))
})

test_that("boot_t_test with tr > 0: paired label includes tr", {
  set.seed(123)
  tr_val <- 0.2
  data(sleep)
  x_sleep <- sleep$extra[sleep$group == 1]
  y_sleep <- sleep$extra[sleep$group == 2]
  res <- hush(boot_t_test(x_sleep, y_sleep, paired = TRUE,
                          mu = 0, tr = tr_val, R = 99))

  expect_equal(length(res$estimate), 1)
  expect_true(grepl("x - y", names(res$estimate)))
  expect_true(grepl(as.character(tr_val), names(res$estimate)))
})

test_that("boot_t_test: sample_size for tr > 0", {
  set.seed(123)
  # Two-sample
  res_two <- hush(boot_t_test(1:10, y = c(7:20), mu = 0, tr = 0.1, R = 99))
  expect_equal(res_two$sample_size, c(nx = 10L, ny = 14L))

  # One-sample
  res_one <- hush(boot_t_test(1:20, mu = 0, tr = 0.1, R = 99))
  expect_equal(res_one$sample_size, c(n = 20L))

  # Paired
  data(sleep)
  x_sleep <- sleep$extra[sleep$group == 1]
  y_sleep <- sleep$extra[sleep$group == 2]
  res_paired <- hush(boot_t_test(x_sleep, y_sleep, paired = TRUE,
                                 mu = 0, tr = 0.1, R = 99))
  expect_equal(res_paired$sample_size, c(n = 10L))
})

test_that("boot_t_test: formula interface substitutes group names", {
  set.seed(123)
  data(sleep)
  res <- hush(boot_t_test(extra ~ group, data = sleep, mu = 0, R = 99))

  nms <- names(res$estimate)
  expect_true(grepl("'1'", nms[1]))
  expect_true(grepl("'2'", nms[2]))
  expect_true(grepl("'1' - '2'", nms[3]))
  expect_equal(names(res$sample_size), c("1", "2"))

  # With trimming
  set.seed(123)
  res_tr <- hush(boot_t_test(extra ~ group, data = sleep,
                             mu = 0, tr = 0.1, R = 99))
  nms_tr <- names(res_tr$estimate)
  expect_true(grepl("trimmed mean of '1'", nms_tr[1]))
  expect_true(grepl("trimmed mean of '2'", nms_tr[2]))
  expect_true(grepl("'1' - '2'", nms_tr[3]))
})

test_that("boot_t_test: null.value names include tr for tr > 0", {
  set.seed(123)
  tr_val <- 0.1
  res <- hush(boot_t_test(1:10, y = c(7:20),
                          mu = c(-3, 3), alternative = "equivalence",
                          tr = tr_val, R = 99))
  expect_true(grepl(as.character(tr_val), names(res$null.value)[1]))
})

# ===========================================================================
# perm_t_test tests
# ===========================================================================

test_that("perm_t_test: two-sample estimate has 3 elements", {
  set.seed(123)
  res <- hush(perm_t_test(1:10, y = c(7:20), mu = 0, R = 99))

  expect_equal(length(res$estimate), 3)
  expect_equal(names(res$estimate)[1], "mean of group x")
  expect_equal(names(res$estimate)[2], "mean of group y")
  expect_equal(names(res$estimate)[3], "mean difference (x - y)")
  expect_equal(unname(res$estimate[3]),
               unname(res$estimate[1] - res$estimate[2]))
})

test_that("perm_t_test: paired label with and without trimming", {
  set.seed(123)
  data(sleep)
  x_sleep <- sleep$extra[sleep$group == 1]
  y_sleep <- sleep$extra[sleep$group == 2]

  # Without trimming
  res <- hush(perm_t_test(x_sleep, y_sleep, paired = TRUE, mu = 0, R = 99))
  expect_equal(names(res$estimate), "mean of the differences (z = x - y)")

  # With trimming
  set.seed(123)
  tr_val <- 0.1
  res_tr <- hush(perm_t_test(x_sleep, y_sleep, paired = TRUE,
                             mu = 0, tr = tr_val, R = 99))
  expect_true(grepl("x - y", names(res_tr$estimate)))
  expect_true(grepl(as.character(tr_val), names(res_tr$estimate)))
})

test_that("perm_t_test: trimmed two-sample includes tr in label", {
  set.seed(123)
  tr_val <- 0.1
  res <- hush(perm_t_test(1:10, y = c(7:20), mu = 0, tr = tr_val, R = 99))

  expect_equal(length(res$estimate), 3)
  expect_true(grepl("trimmed mean of x", names(res$estimate)[1]))
  expect_true(grepl("trimmed mean of y", names(res$estimate)[2]))
  expect_true(grepl(as.character(tr_val), names(res$estimate)[3]))
})

test_that("perm_t_test: sample_size correctness", {
  set.seed(123)
  # Two-sample
  res <- hush(perm_t_test(1:10, y = c(7:20), mu = 0, R = 99))
  expect_equal(res$sample_size, c(nx = 10L, ny = 14L))

  # One-sample
  res_one <- hush(perm_t_test(1:20, mu = 0, R = 99))
  expect_equal(res_one$sample_size, c(n = 20L))

  # Paired
  data(sleep)
  x_sleep <- sleep$extra[sleep$group == 1]
  y_sleep <- sleep$extra[sleep$group == 2]
  res_paired <- hush(perm_t_test(x_sleep, y_sleep, paired = TRUE,
                                 mu = 0, R = 99))
  expect_equal(res_paired$sample_size, c(n = 10L))
})

test_that("perm_t_test: formula interface substitutes group names", {
  set.seed(123)
  data(sleep)
  res <- hush(perm_t_test(extra ~ group, data = sleep, mu = 0, R = 99))

  nms <- names(res$estimate)
  expect_true(grepl("'1'", nms[1]))
  expect_true(grepl("'2'", nms[2]))
  expect_true(grepl("'1' - '2'", nms[3]))
  expect_equal(names(res$sample_size), c("1", "2"))
})

# ===========================================================================
# hodges_lehmann tests
# ===========================================================================

test_that("hodges_lehmann: estimate labels", {
  # One-sample
  res_one <- hush(hodges_lehmann(1:20, mu = 0))
  expect_equal(names(res_one$estimate), "(pseudo)median of x")

  # Paired
  before <- c(5.1, 4.8, 6.2, 5.7, 6.0, 5.5, 4.9, 5.8)
  after <- c(5.6, 5.2, 6.7, 6.1, 6.5, 5.8, 5.3, 6.2)
  res_paired <- hush(hodges_lehmann(before, after, paired = TRUE, mu = 0))
  expect_equal(names(res_paired$estimate), "Hodges-Lehmann estimate (z = x - y)")

  # Two-sample
  set.seed(123)
  x <- rnorm(20)
  y <- rnorm(20, mean = 0.5)
  res_two <- hush(hodges_lehmann(x, y, mu = 0))
  expect_equal(names(res_two$estimate), "Hodges-Lehmann estimate (x - y)")
  expect_equal(length(res_two$estimate), 1)  # always single value
})

test_that("hodges_lehmann: sample_size correctness", {
  # One-sample
  res_one <- hush(hodges_lehmann(1:20, mu = 0))
  expect_equal(res_one$sample_size, c(n = 20L))

  # Two-sample
  set.seed(123)
  x <- rnorm(15)
  y <- rnorm(20)
  res_two <- hush(hodges_lehmann(x, y, mu = 0))
  expect_equal(res_two$sample_size, c(nx = 15L, ny = 20L))

  # Paired
  before <- c(5.1, 4.8, 6.2, 5.7, 6.0, 5.5, 4.9, 5.8)
  after <- c(5.6, 5.2, 6.7, 6.1, 6.5, 5.8, 5.3, 6.2)
  res_paired <- hush(hodges_lehmann(before, after, paired = TRUE, mu = 0))
  expect_equal(res_paired$sample_size, c(n = 8L))
})

test_that("hodges_lehmann: formula interface substitutes group names with quoting", {
  data(sleep)
  res <- hush(hodges_lehmann(extra ~ group, data = sleep, mu = 0))

  # sleep$group levels "1"/"2" should be quoted
  expect_true(grepl("'1' - '2'", names(res$estimate)))
  expect_equal(names(res$sample_size), c("1", "2"))
})

test_that("hodges_lehmann: equivalence/MET paths work with updated labels", {
  set.seed(123)
  x <- rnorm(20)
  y <- rnorm(20, mean = 0.5)

  # Equivalence (asymptotic only — permutation not supported for equivalence/MET)
  res_eq <- hush(hodges_lehmann(x, y, alternative = "equivalence",
                                mu = c(-1, 1)))
  expect_equal(names(res_eq$estimate), "Hodges-Lehmann estimate (x - y)")
  expect_true(!is.null(res_eq$sample_size))

  # MET (asymptotic only)
  res_met <- hush(hodges_lehmann(x, y, alternative = "minimal.effect",
                                 mu = c(-1, 1)))
  expect_equal(names(res_met$estimate), "Hodges-Lehmann estimate (x - y)")
})

# ===========================================================================
# Cross-function consistency tests
# ===========================================================================

test_that("Formula interface group names are consistent across functions", {
  data(sleep)

  res_simple <- hush(simple_htest(extra ~ group, data = sleep, mu = 0))
  set.seed(1)
  res_boot <- hush(boot_t_test(extra ~ group, data = sleep, mu = 0, R = 99))
  set.seed(1)
  res_perm <- hush(perm_t_test(extra ~ group, data = sleep, mu = 0, R = 99))
  res_hl <- hush(hodges_lehmann(extra ~ group, data = sleep, mu = 0))

  # All should have group-named sample sizes
  for (res in list(res_simple, res_boot, res_perm, res_hl)) {
    expect_equal(names(res$sample_size), c("1", "2"),
                 info = paste("Failed for", res$method))
  }

  # Mean-based tests should all have 3 estimates with quoted group names
  for (res in list(res_simple, res_boot, res_perm)) {
    expect_equal(length(res$estimate), 3,
                 info = paste("Failed for", res$method))
    expect_true(grepl("'1' - '2'", names(res$estimate)[3]),
                info = paste("Failed for", res$method))
  }

  # HL should have 1 estimate with quoted group names
  expect_equal(length(res_hl$estimate), 1)
  expect_true(grepl("'1' - '2'", names(res_hl$estimate)))
})

test_that("quote_if_numeric quotes numeric names and leaves text names alone", {
  # Numeric names get quoted
  expect_equal(TOSTER:::quote_if_numeric("0"), "'0'")
  expect_equal(TOSTER:::quote_if_numeric("1"), "'1'")
  expect_equal(TOSTER:::quote_if_numeric("3.14"), "'3.14'")
  expect_equal(TOSTER:::quote_if_numeric("-2"), "'-2'")

  # Non-numeric names are unchanged
  expect_equal(TOSTER:::quote_if_numeric("auto"), "auto")
  expect_equal(TOSTER:::quote_if_numeric("group_A"), "group_A")
  expect_equal(TOSTER:::quote_if_numeric("1a"), "1a")
})

test_that("relabel_for_formula is idempotent on non-matching patterns", {
  # Create a fake result with labels that don't match patterns
  fake <- list(estimate = c(a = 1, b = 2), sample_size = c(n = 10))
  result <- relabel_for_formula(fake, c("groupA", "groupB"))
  # Should be unchanged since "of x"/"of y" patterns don't appear

  expect_equal(names(result$estimate), c("a", "b"))
  # sample_size with length 1 should not be renamed
  expect_equal(names(result$sample_size), "n")
})
