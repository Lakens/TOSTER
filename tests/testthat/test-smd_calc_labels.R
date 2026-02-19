# Tests for smd_calc label logic (estimate names, method strings, null.value names)

hush = function(code) {
  sink("NUL") # use /dev/null in UNIX
  tmp = code
  sink()
  return(tmp)
}

set.seed(42)
group1 <- rnorm(30, mean = 100, sd = 15)
group2 <- rnorm(30, mean = 110, sd = 18)
before <- c(5.1, 4.8, 6.2, 5.7, 6.0, 5.5, 4.9, 5.8)
after  <- c(5.6, 5.2, 6.7, 6.1, 6.5, 5.8, 5.3, 6.2)

# --- Estimate name tests (SMD (d/g[subscript]) format) ---

test_that("Two-sample var.equal=FALSE bias_correction=FALSE => SMD (d[av])", {
  res <- smd_calc(x = group1, y = group2,
                  var.equal = FALSE, bias_correction = FALSE)
  expect_equal(names(res$estimate), "SMD (d[av])")
})

test_that("Two-sample var.equal=TRUE bias_correction=FALSE => SMD (d[s])", {
  res <- smd_calc(x = group1, y = group2,
                  var.equal = TRUE, bias_correction = FALSE)
  expect_equal(names(res$estimate), "SMD (d[s])")
})

test_that("Two-sample var.equal=FALSE bias_correction=TRUE => SMD (g[av])", {
  res <- smd_calc(x = group1, y = group2,
                  var.equal = FALSE, bias_correction = TRUE)
  expect_equal(names(res$estimate), "SMD (g[av])")
})

test_that("Two-sample var.equal=TRUE bias_correction=TRUE => SMD (g[s])", {
  res <- smd_calc(x = group1, y = group2,
                  var.equal = TRUE, bias_correction = TRUE)
  expect_equal(names(res$estimate), "SMD (g[s])")
})

test_that("Two-sample denom='pooled' bias_correction=FALSE => SMD (d[s])", {
  res <- hush(smd_calc(x = group1, y = group2,
                       denom = "pooled", bias_correction = FALSE))
  expect_equal(names(res$estimate), "SMD (d[s])")
})

test_that("Two-sample denom='avg' bias_correction=FALSE => SMD (d[av])", {
  res <- hush(smd_calc(x = group1, y = group2,
                       denom = "avg", bias_correction = FALSE))
  expect_equal(names(res$estimate), "SMD (d[av])")
})

test_that("Two-sample glass='glass1' default => SMD (g[x])", {
  res <- smd_calc(x = group1, y = group2, glass = "glass1")
  expect_equal(names(res$estimate), "SMD (g[x])")
})

test_that("Two-sample glass='glass2' default => SMD (g[y])", {
  res <- smd_calc(x = group1, y = group2, glass = "glass2")
  expect_equal(names(res$estimate), "SMD (g[y])")
})

test_that("Two-sample glass='glass1' bias_correction=FALSE => SMD (d[x])", {
  res <- smd_calc(x = group1, y = group2, glass = "glass1",
                  bias_correction = FALSE)
  expect_equal(names(res$estimate), "SMD (d[x])")
})

test_that("Paired default => SMD (g[z])", {
  res <- smd_calc(x = before, y = after, paired = TRUE)
  expect_equal(names(res$estimate), "SMD (g[z])")
})

test_that("Paired rm_correction=TRUE => SMD (g[rm])", {
  res <- smd_calc(x = before, y = after, paired = TRUE, rm_correction = TRUE)
  expect_equal(names(res$estimate), "SMD (g[rm])")
})

test_that("One-sample default => SMD (g[z])", {
  res <- smd_calc(x = group1)
  expect_equal(names(res$estimate), "SMD (g[z])")
})

# --- Method string SD notation tests ---

test_that("Two-sample var.equal=FALSE method contains SD_avg", {
  res <- smd_calc(x = group1, y = group2, var.equal = FALSE)
  expect_true(grepl("SD_avg", res$method, fixed = TRUE))
})

test_that("Two-sample var.equal=TRUE method contains SD_pooled", {
  res <- smd_calc(x = group1, y = group2, var.equal = TRUE)
  expect_true(grepl("SD_pooled", res$method, fixed = TRUE))
})

test_that("Two-sample denom='pooled' method contains SD_pooled", {
  res <- hush(smd_calc(x = group1, y = group2, denom = "pooled"))
  expect_true(grepl("SD_pooled", res$method, fixed = TRUE))
})

test_that("Two-sample denom='avg' method contains SD_avg", {
  res <- hush(smd_calc(x = group1, y = group2, denom = "avg"))
  expect_true(grepl("SD_avg", res$method, fixed = TRUE))
})

test_that("Paired default method contains SD_z", {
  res <- smd_calc(x = before, y = after, paired = TRUE)
  expect_true(grepl("SD_z", res$method, fixed = TRUE))
})

test_that("Paired rm_correction=TRUE method contains SD_rm", {
  res <- smd_calc(x = before, y = after, paired = TRUE, rm_correction = TRUE)
  expect_true(grepl("SD_rm", res$method, fixed = TRUE))
})

test_that("One-sample default method contains SD_z", {
  res <- smd_calc(x = group1)
  expect_true(grepl("SD_z", res$method, fixed = TRUE))
})

test_that("Two-sample glass='glass1' method contains SD_x", {
  res <- smd_calc(x = group1, y = group2, glass = "glass1")
  expect_true(grepl("SD_x", res$method, fixed = TRUE))
})

# --- Method string contains SMD label ---

test_that("Method string contains merged SMD label with notation", {
  res <- smd_calc(x = group1, y = group2, var.equal = FALSE,
                  bias_correction = FALSE)
  expect_true(grepl("SMD (d[av]=(x-y)/SD_avg)", res$method, fixed = TRUE))

  res2 <- smd_calc(x = group1, y = group2, var.equal = TRUE,
                   bias_correction = TRUE)
  expect_true(grepl("SMD (g[s]=(x-y)/SD_pooled)", res2$method, fixed = TRUE))
})

# --- Formula method group name substitution tests ---

test_that("Formula glass='glass1' substitutes group names in estimate and method", {
  df <- data.frame(
    value = c(group1, group2),
    group = factor(rep(c("A", "B"), each = 30))
  )
  res <- smd_calc(formula = value ~ group, data = df, glass = "glass1")
  expect_equal(names(res$estimate), "SMD (g[A])")
  expect_true(grepl("g[A]", res$method, fixed = TRUE))
  expect_true(grepl("SD_A", res$method, fixed = TRUE))
})

test_that("Formula glass='glass2' substitutes group names in estimate and method", {
  df <- data.frame(
    value = c(group1, group2),
    group = factor(rep(c("A", "B"), each = 30))
  )
  res <- smd_calc(formula = value ~ group, data = df, glass = "glass2")
  expect_equal(names(res$estimate), "SMD (g[B])")
  expect_true(grepl("g[B]", res$method, fixed = TRUE))
  expect_true(grepl("SD_B", res$method, fixed = TRUE))
})

test_that("Formula glass='glass1' bias_correction=FALSE shows d not g", {
  df <- data.frame(
    value = c(group1, group2),
    group = factor(rep(c("A", "B"), each = 30))
  )
  res <- smd_calc(formula = value ~ group, data = df, glass = "glass1",
                  bias_correction = FALSE)
  expect_equal(names(res$estimate), "SMD (d[A])")
  expect_true(grepl("d[A]", res$method, fixed = TRUE))
})

test_that("Formula var.equal=FALSE keeps bracket subscript (no group name swap)", {
  df <- data.frame(
    value = c(group1, group2),
    group = factor(rep(c("A", "B"), each = 30))
  )
  res <- smd_calc(formula = value ~ group, data = df, var.equal = FALSE)
  # Subscript [av] should remain unchanged (not a group name)
  expect_equal(names(res$estimate), "SMD (g[av])")
  expect_true(grepl("SD_avg", res$method, fixed = TRUE))
  expect_true(grepl("(A-B)", res$method, fixed = TRUE))
})

# --- null.value name tests ---

test_that("Standard alternative null.value names are 'SMD'", {
  res <- smd_calc(x = group1, y = group2,
                  alternative = "two.sided", null.value = 0)
  expect_equal(names(res$null.value), "SMD")
})

test_that("Equivalence alternative null.value names are bounds", {
  res <- smd_calc(x = group1, y = group2,
                  alternative = "equivalence", null.value = c(-0.5, 0.5))
  expect_equal(names(res$null.value), c("lower bound", "upper bound"))
})

# --- Trimming note preserved ---

test_that("Trimmed estimate label includes trim note with subscript", {
  res <- hush(smd_calc(x = group1, y = group2,
                       tr = 0.1, var.equal = FALSE, bias_correction = FALSE))
  expect_equal(names(res$estimate), "SMD (d[av], 10% trimmed)")
})

# --- Consistency between smd_calc and boot_smd_calc ---

test_that("boot_smd_calc labels match smd_calc labels", {
  skip_on_cran()
  set.seed(42)
  x <- rnorm(20, mean = 5, sd = 2)
  y <- rnorm(20, mean = 6, sd = 2)

  res_smd <- smd_calc(x = x, y = y, bias_correction = FALSE)
  res_boot <- boot_smd_calc(x = x, y = y, bias_correction = FALSE, R = 99)
  expect_equal(names(res_smd$estimate), names(res_boot$estimate))

  res_smd2 <- smd_calc(x = x, y = y, bias_correction = TRUE)
  res_boot2 <- boot_smd_calc(x = x, y = y, bias_correction = TRUE, R = 99)
  expect_equal(names(res_smd2$estimate), names(res_boot2$estimate))

  res_smd3 <- smd_calc(x = x, y = y, paired = TRUE)
  res_boot3 <- boot_smd_calc(x = x, y = y, paired = TRUE, R = 99)
  expect_equal(names(res_smd3$estimate), names(res_boot3$estimate))
})
