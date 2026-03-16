# Tests for denom parameter in smd_calc and boot_smd_calc

hush = function(code) {
  sink(nullfile())
  tmp = code
  sink()
  return(tmp)
}

# Test data
set.seed(123)
x_ind <- rnorm(30, mean = 100, sd = 15)
y_ind <- rnorm(30, mean = 110, sd = 18)
x_pair <- rnorm(20, mean = 5, sd = 2)
y_pair <- x_pair + rnorm(20, mean = 1, sd = 1)
x_one <- rnorm(25, mean = 3, sd = 2)

# --- 1. Equivalence with existing interface ---

test_that("denom = 'pooled' matches var.equal = TRUE", {
  r1 <- smd_calc(x = x_ind, y = y_ind, denom = "pooled")
  r2 <- smd_calc(x = x_ind, y = y_ind, var.equal = TRUE)
  expect_equal(unname(r1$estimate), unname(r2$estimate))
  expect_equal(r1$conf.int, r2$conf.int)
  expect_equal(r1$stderr, r2$stderr)
})

test_that("denom = 'avg' matches var.equal = FALSE", {
  r1 <- smd_calc(x = x_ind, y = y_ind, denom = "avg")
  r2 <- smd_calc(x = x_ind, y = y_ind, var.equal = FALSE)
  expect_equal(unname(r1$estimate), unname(r2$estimate))
  expect_equal(r1$conf.int, r2$conf.int)
  expect_equal(r1$stderr, r2$stderr)
})

test_that("denom = 'rm' matches rm_correction = TRUE for paired", {
  r1 <- smd_calc(x = x_pair, y = y_pair, paired = TRUE, denom = "rm")
  r2 <- smd_calc(x = x_pair, y = y_pair, paired = TRUE, rm_correction = TRUE)
  expect_equal(unname(r1$estimate), unname(r2$estimate))
  expect_equal(r1$conf.int, r2$conf.int)
  expect_equal(r1$stderr, r2$stderr)
})

test_that("denom = 'z' matches rm_correction = FALSE for paired", {
  r1 <- smd_calc(x = x_pair, y = y_pair, paired = TRUE, denom = "z")
  r2 <- smd_calc(x = x_pair, y = y_pair, paired = TRUE, rm_correction = FALSE)
  expect_equal(unname(r1$estimate), unname(r2$estimate))
  expect_equal(r1$conf.int, r2$conf.int)
  expect_equal(r1$stderr, r2$stderr)
})

test_that("denom = 'glass1' matches glass = 'glass1'", {
  r1 <- smd_calc(x = x_ind, y = y_ind, denom = "glass1")
  r2 <- smd_calc(x = x_ind, y = y_ind, glass = "glass1")
  expect_equal(unname(r1$estimate), unname(r2$estimate))
  expect_equal(r1$conf.int, r2$conf.int)
  expect_equal(r1$stderr, r2$stderr)
})

test_that("denom = 'glass2' matches glass = 'glass2'", {
  r1 <- smd_calc(x = x_ind, y = y_ind, denom = "glass2")
  r2 <- smd_calc(x = x_ind, y = y_ind, glass = "glass2")
  expect_equal(unname(r1$estimate), unname(r2$estimate))
  expect_equal(r1$conf.int, r2$conf.int)
  expect_equal(r1$stderr, r2$stderr)
})

test_that("denom = 'glass1' matches glass = 'glass1' for paired", {
  r1 <- smd_calc(x = x_pair, y = y_pair, paired = TRUE, denom = "glass1")
  r2 <- smd_calc(x = x_pair, y = y_pair, paired = TRUE, glass = "glass1")
  expect_equal(unname(r1$estimate), unname(r2$estimate))
  expect_equal(r1$conf.int, r2$conf.int)
  expect_equal(r1$stderr, r2$stderr)
})

test_that("denom = 'z' works for one-sample", {
  r1 <- smd_calc(x = x_one, denom = "z")
  r2 <- smd_calc(x = x_one)
  expect_equal(unname(r1$estimate), unname(r2$estimate))
  expect_equal(r1$conf.int, r2$conf.int)
})

# --- 2. Design-validity errors ---

test_that("denom = 'z' errors for independent samples", {
  expect_error(
    smd_calc(x = x_ind, y = y_ind, denom = "z"),
    "denom = 'z' is not valid for independent samples designs."
  )
})

test_that("denom = 'rm' errors for independent samples", {
  expect_error(
    smd_calc(x = x_ind, y = y_ind, denom = "rm"),
    "denom = 'rm' is only valid for paired samples designs."
  )
})

test_that("denom = 'pooled' errors for paired samples", {
  expect_error(
    smd_calc(x = x_pair, y = y_pair, paired = TRUE, denom = "pooled"),
    "denom = 'pooled' is only valid for independent samples designs."
  )
})

test_that("denom = 'avg' errors for paired samples", {
  expect_error(
    smd_calc(x = x_pair, y = y_pair, paired = TRUE, denom = "avg"),
    "denom = 'avg' is only valid for independent samples designs."
  )
})

test_that("denom = 'rm' errors for one-sample", {
  expect_error(
    smd_calc(x = x_one, denom = "rm"),
    "denom = 'rm' is only valid for paired samples designs."
  )
})

test_that("denom = 'glass1' errors for one-sample", {
  expect_error(
    smd_calc(x = x_one, denom = "glass1"),
    "denom = 'glass1' is not valid for one-sample designs."
  )
})

test_that("denom = 'glass2' errors for one-sample", {
  expect_error(
    smd_calc(x = x_one, denom = "glass2"),
    "denom = 'glass2' is not valid for one-sample designs."
  )
})

test_that("denom = 'pooled' errors for one-sample", {
  expect_error(
    smd_calc(x = x_one, denom = "pooled"),
    "denom = 'pooled' is only valid for independent samples designs."
  )
})

test_that("denom = 'avg' errors for one-sample", {
  expect_error(
    smd_calc(x = x_one, denom = "avg"),
    "denom = 'avg' is only valid for independent samples designs."
  )
})

# --- 3. Override messages ---

test_that("denom = 'pooled' with var.equal = FALSE emits message", {
  expect_message(
    smd_calc(x = x_ind, y = y_ind, denom = "pooled", var.equal = FALSE),
    "denom = 'pooled' overrides var.equal to TRUE."
  )
})

test_that("denom = 'rm' with rm_correction = FALSE emits message", {
  expect_message(
    smd_calc(x = x_pair, y = y_pair, paired = TRUE, denom = "rm", rm_correction = FALSE),
    "denom = 'rm' overrides rm_correction to TRUE."
  )
})

test_that("denom = 'glass1' with glass = 'glass2' emits message", {
  expect_message(
    smd_calc(x = x_ind, y = y_ind, denom = "glass1", glass = "glass2"),
    "denom = 'glass1' overrides glass argument."
  )
})

test_that("denom = 'avg' with glass = 'glass1' emits message", {
  expect_message(
    smd_calc(x = x_ind, y = y_ind, denom = "avg", glass = "glass1"),
    "denom = 'avg' overrides glass argument."
  )
})

test_that("denom = 'z' with rm_correction = TRUE emits message", {
  expect_message(
    smd_calc(x = x_pair, y = y_pair, paired = TRUE, denom = "z", rm_correction = TRUE),
    "denom = 'z' overrides rm_correction to FALSE."
  )
})

test_that("denom = 'avg' with var.equal = TRUE emits message", {
  expect_message(
    smd_calc(x = x_ind, y = y_ind, denom = "avg", var.equal = TRUE),
    "denom = 'avg' overrides var.equal to FALSE."
  )
})

# --- No message when no conflict ---

test_that("denom = 'pooled' with var.equal = TRUE does NOT emit message", {
  expect_no_message(
    smd_calc(x = x_ind, y = y_ind, denom = "pooled", var.equal = TRUE)
  )
})

test_that("denom = 'pooled' without explicit var.equal does NOT emit message", {
  expect_no_message(
    smd_calc(x = x_ind, y = y_ind, denom = "pooled")
  )
})

test_that("denom = 'rm' with rm_correction = TRUE does NOT emit message", {
  expect_no_message(
    smd_calc(x = x_pair, y = y_pair, paired = TRUE, denom = "rm", rm_correction = TRUE)
  )
})

# --- 4. Auto preserves existing behavior ---

test_that("denom = 'auto' is identical to omitting denom (independent)", {
  r1 <- smd_calc(x = x_ind, y = y_ind, denom = "auto", var.equal = TRUE)
  r2 <- smd_calc(x = x_ind, y = y_ind, var.equal = TRUE)
  expect_equal(unname(r1$estimate), unname(r2$estimate))
  expect_equal(r1$conf.int, r2$conf.int)
})

test_that("denom = 'auto' is identical to omitting denom (paired)", {
  r1 <- smd_calc(x = x_pair, y = y_pair, paired = TRUE, denom = "auto", rm_correction = TRUE)
  r2 <- smd_calc(x = x_pair, y = y_pair, paired = TRUE, rm_correction = TRUE)
  expect_equal(unname(r1$estimate), unname(r2$estimate))
  expect_equal(r1$conf.int, r2$conf.int)
})

test_that("denom = 'auto' is identical to omitting denom (one-sample)", {
  r1 <- smd_calc(x = x_one, denom = "auto")
  r2 <- smd_calc(x = x_one)
  expect_equal(unname(r1$estimate), unname(r2$estimate))
  expect_equal(r1$conf.int, r2$conf.int)
})

# --- 5. bias_correction independence ---

test_that("denom does not affect bias_correction", {
  r1 <- smd_calc(x = x_ind, y = y_ind, denom = "pooled", bias_correction = TRUE)
  r2 <- smd_calc(x = x_ind, y = y_ind, denom = "pooled", bias_correction = FALSE)
  # Hedges' g and Cohen's d should differ

  expect_false(unname(r1$estimate) == unname(r2$estimate))
})

# --- 6. Bootstrap: denom equivalence ---

test_that("boot_smd_calc denom = 'pooled' matches var.equal = TRUE", {
  set.seed(42)
  r1 <- boot_smd_calc(x = x_ind, y = y_ind, denom = "pooled", R = 199)
  set.seed(42)
  r2 <- boot_smd_calc(x = x_ind, y = y_ind, var.equal = TRUE, R = 199)
  expect_equal(unname(r1$estimate), unname(r2$estimate))
  expect_equal(r1$conf.int, r2$conf.int)
})

test_that("boot_smd_calc denom = 'avg' matches var.equal = FALSE", {
  set.seed(42)
  r1 <- boot_smd_calc(x = x_ind, y = y_ind, denom = "avg", R = 199)
  set.seed(42)
  r2 <- boot_smd_calc(x = x_ind, y = y_ind, var.equal = FALSE, R = 199)
  expect_equal(unname(r1$estimate), unname(r2$estimate))
  expect_equal(r1$conf.int, r2$conf.int)
})

test_that("boot_smd_calc denom = 'rm' matches rm_correction = TRUE (paired)", {
  set.seed(42)
  r1 <- boot_smd_calc(x = x_pair, y = y_pair, paired = TRUE, denom = "rm", R = 199)
  set.seed(42)
  r2 <- boot_smd_calc(x = x_pair, y = y_pair, paired = TRUE, rm_correction = TRUE, R = 199)
  expect_equal(unname(r1$estimate), unname(r2$estimate))
  expect_equal(r1$conf.int, r2$conf.int)
})

# --- 7. Bootstrap: design-validity errors ---

test_that("boot_smd_calc denom errors for invalid designs", {
  expect_error(
    boot_smd_calc(x = x_ind, y = y_ind, denom = "z", R = 99),
    "denom = 'z' is not valid for independent samples designs."
  )
  expect_error(
    boot_smd_calc(x = x_pair, y = y_pair, paired = TRUE, denom = "pooled", R = 99),
    "denom = 'pooled' is only valid for independent samples designs."
  )
  expect_error(
    boot_smd_calc(x = x_one, denom = "rm", R = 99),
    "denom = 'rm' is only valid for paired samples designs."
  )
})

# --- 8. Bootstrap: override message emitted once ---

test_that("boot_smd_calc emits message only once for override", {
  msgs <- character(0)
  withCallingHandlers(
    boot_smd_calc(x = x_ind, y = y_ind, denom = "pooled", var.equal = FALSE, R = 99),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  override_msgs <- grep("overrides", msgs, value = TRUE)
  expect_length(override_msgs, 1)
})

# --- 9. data.frame output also works with denom ---

test_that("smd_calc data.frame output works with denom", {
  r1 <- smd_calc(x = x_ind, y = y_ind, denom = "pooled", output = "data.frame")
  r2 <- smd_calc(x = x_ind, y = y_ind, var.equal = TRUE, output = "data.frame")
  expect_equal(r1$estimate, r2$estimate)
  expect_equal(r1$SE, r2$SE)
  expect_equal(r1$lower.ci, r2$lower.ci)
  expect_equal(r1$upper.ci, r2$upper.ci)
})

# --- 10. Formula interface passes denom through ---

test_that("smd_calc formula interface works with denom", {
  df <- data.frame(
    value = c(x_ind, y_ind),
    group = factor(rep(c("A", "B"), each = 30))
  )
  r1 <- smd_calc(formula = value ~ group, data = df, denom = "pooled")
  r2 <- smd_calc(formula = value ~ group, data = df, var.equal = TRUE)
  expect_equal(unname(r1$estimate), unname(r2$estimate))
})
