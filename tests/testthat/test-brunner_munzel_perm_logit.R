# brunner_munzel(test_method = "perm_logit") ----

bm_pl <- function(args, mu = NULL, seed = 1) {
  if (!is.null(mu)) args$mu <- mu
  if (is.null(args$test_method)) args$test_method <- "perm_logit"
  set.seed(seed)
  suppressWarnings(suppressMessages(do.call(brunner_munzel, args)))
}

# Independent recomputation of the logit-scale permutation test ----
# (fully enumerated, so the permutation set is identical to the package's)

logit_stat <- function(pd, se, eps) {
  pl <- min(max(pd, eps), 1 - eps)
  c(L = qlogis(pl), seL = se / (pl * (1 - pl)))
}

two_stats <- function(x, y) {
  nx <- length(x); ny <- length(y); N <- nx + ny
  rxy <- rank(c(x, y)); rx <- rank(x); ry <- rank(y)
  pl2 <- (rxy[1:nx] - rx) / ny
  pl1 <- (rxy[(nx + 1):N] - ry) / nx
  V <- N * (var(pl2) / nx + var(pl1) / ny)
  if (V == 0) V <- N * 0.5 / (nx * ny)^2
  c(pd = mean(pl2), se = sqrt(V / N))
}

paired_stats <- function(x, y) {
  n <- length(x)
  rx <- rank(c(y, x))
  BM2 <- (rx[(n + 1):(2 * n)] - rank(x)) / n
  BM3 <- (rx[1:n] - rank(y)) / n - BM2
  v <- (sum(BM3^2) - n * mean(BM3)^2) / (n - 1)
  if (v == 0) v <- 1 / n
  c(pd = mean(BM2), se = sqrt(v / n))
}

manual_perm_logit <- function(obs, perm_stats, eps, mu, alpha = 0.05) {
  o <- logit_stat(obs[["pd"]], obs[["se"]], eps)
  Tperm <- apply(perm_stats, 1, function(s) {
    l <- logit_stat(s[["pd"]], s[["se"]], eps)
    l[["L"]] / l[["seL"]]
  })
  t_obs <- (o[["L"]] - qlogis(mu)) / o[["seL"]]
  b <- sum(abs(Tperm) >= abs(t_obs) - 1e-8 * max(1, abs(t_obs)))
  crit <- perm_crit(Tperm, alpha, "exact", "abs")
  list(p = b / length(Tperm),
       ci = plogis(o[["L"]] + c(-1, 1) * o[["seL"]] * crit))
}

test_that("perm_logit matches an independent calculation (two-sample, enumerated)", {
  skip_on_cran()
  set.seed(10)
  x <- rnorm(4, 1)
  y <- rnorm(5)
  z <- c(x, y)
  idx <- utils::combn(9, 4)   # 126 permutations, enumerated by the package
  ps <- t(apply(idx, 2, function(i) two_stats(z[i], z[-i])))
  man <- manual_perm_logit(two_stats(x, y), ps, eps = 0.5 / 20, mu = 0.5)

  res <- bm_pl(list(x = x, y = y, R = 1000))
  expect_equal(res$parameter[[1]], 126)
  expect_equal(res$p.value, man$p)
  expect_equal(as.numeric(res$conf.int), man$ci, tolerance = 1e-6)
})

test_that("perm_logit matches an independent calculation (paired, enumerated)", {
  skip_on_cran()
  set.seed(11)
  x <- rnorm(7, 0.5)
  y <- x + rnorm(7, -0.3)
  S <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), 7)))  # 128 swaps
  ps <- t(apply(S, 1, function(s) paired_stats(ifelse(s, y, x), ifelse(s, x, y))))
  man <- manual_perm_logit(paired_stats(x, y), ps, eps = 0.5 / 49, mu = 0.5)

  res <- bm_pl(list(x = x, y = y, paired = TRUE))
  expect_equal(res$parameter[[1]], 128)
  expect_equal(res$p.value, man$p)
  expect_equal(as.numeric(res$conf.int), man$ci, tolerance = 1e-6)
})

# CI / p-value duality ----

expect_pl_duality <- function(args, seed = 1) {
  res <- bm_pl(args, mu = 0.5, seed = seed)
  ci <- res$conf.int
  eps <- 1e-7
  checked <- 0
  if (ci[1] > 1e-6) {
    expect_gt(bm_pl(args, ci[1] + eps, seed)$p.value, 0.05)
    expect_lte(bm_pl(args, ci[1] - eps, seed)$p.value, 0.05)
    checked <- checked + 1
  }
  if (ci[2] < 1 - 1e-6) {
    expect_gt(bm_pl(args, ci[2] - eps, seed)$p.value, 0.05)
    expect_lte(bm_pl(args, ci[2] + eps, seed)$p.value, 0.05)
    checked <- checked + 1
  }
  expect_gt(checked, 0)
}

set.seed(1)
pl_x <- rnorm(8, 0.8, 2)
pl_y <- rnorm(20, 0, 1)
pl_px <- rnorm(15)
pl_py <- pl_px + rnorm(15, -0.4, 1.2)

test_that("perm_logit CI agrees with its p-value", {
  skip_on_cran()
  for (alt in c("two.sided", "less", "greater")) {
    expect_pl_duality(list(x = pl_x, y = pl_y, alternative = alt, R = 999))
    expect_pl_duality(list(x = pl_px, y = pl_py, paired = TRUE,
                           alternative = alt, R = 999))
    expect_pl_duality(list(x = pl_px[1:10], y = pl_py[1:10], paired = TRUE,
                           alternative = alt))
  }
})

test_that("perm_logit equivalence/MET decisions agree with the 1 - 2*alpha CI", {
  skip_on_cran()
  eps <- 1e-7
  for (paired in c(FALSE, TRUE)) {
    args <- if (paired) list(x = pl_px, y = pl_py, paired = TRUE, R = 999) else
      list(x = pl_x, y = pl_y, R = 999)
    run <- function(alt, mu) {
      args$alternative <- alt
      bm_pl(args, mu)
    }
    ci <- run("equivalence", c(0.05, 0.95))$conf.int
    expect_equal(attr(ci, "conf.level"), 0.90)
    expect_equal(run("minimal.effect", c(0.05, 0.95))$conf.int, ci)

    expect_lte(run("equivalence", c(ci[1] - eps, ci[2] + eps))$p.value, 0.05)
    expect_gt(run("equivalence", c(ci[1] + eps, ci[2] + eps))$p.value, 0.05)
    expect_gt(run("equivalence", c(ci[1] - eps, ci[2] - eps))$p.value, 0.05)

    expect_lte(run("minimal.effect", c(ci[2] + eps, min(ci[2] + 0.1, 0.999)))$p.value, 0.05)
    expect_gt(run("minimal.effect", c(ci[2] - eps, min(ci[2] + 0.1, 0.999)))$p.value, 0.05)
    expect_lte(run("minimal.effect", c(max(ci[1] - 0.1, 0.001), ci[1] - eps))$p.value, 0.05)
    expect_gt(run("minimal.effect", c(max(ci[1] - 0.1, 0.001), ci[1] + eps))$p.value, 0.05)
  }
})

# Range preservation and labelling ----

test_that("perm_logit intervals stay strictly inside (0, 1) under complete separation", {
  skip_on_cran()
  x <- c(10, 11, 12, 13, 14, 15)
  y <- c(1, 2, 3, 4, 5, 6)
  res <- bm_pl(list(x = x, y = y, R = 999))
  expect_equal(unname(res$estimate), 1)
  expect_gt(res$conf.int[1], 0)
  expect_lt(res$conf.int[2], 1)
  expect_lt(res$p.value, 0.05)
  # the identity-scale perm interval is clamped at 1 instead
  res_perm <- bm_pl(list(x = x, y = y, R = 999, test_method = "perm"))
  expect_equal(res_perm$conf.int[2], 1)
})

test_that("perm_logit method label and statistic names", {
  res <- bm_pl(list(x = pl_x, y = pl_y, R = 199))
  expect_match(res$method, "\\(logit\\)$")
  expect_equal(names(res$parameter), "N-permutations")
  res_p <- bm_pl(list(x = pl_px, y = pl_py, paired = TRUE, R = 199))
  expect_match(res_p$method, "Paired.*\\(logit\\)$")
})

test_that("perm_logit requires null values strictly inside (0, 1)", {
  expect_error(brunner_munzel(pl_x, pl_y, test_method = "perm_logit",
                              alternative = "less", mu = 1),
               "strictly between 0 and 1")
  expect_error(brunner_munzel(pl_x, pl_y, test_method = "perm_logit",
                              alternative = "equivalence", mu = c(0, 0.6)),
               "strictly between 0 and 1")
})

# Warnings for t and perm with nulls other than 0.5 ----

test_that("t and perm warn for minimal effect tests and nulls other than 0.5", {
  set.seed(2)
  x <- rnorm(20); y <- rnorm(20)
  for (tm in c("t", "perm")) {
    expect_warning(suppressMessages(brunner_munzel(
      x, y, test_method = tm, alternative = "minimal.effect", mu = c(0.3, 0.7),
      R = 199)), "Minimal effect tests")
    expect_warning(suppressMessages(brunner_munzel(
      x, y, test_method = tm, alternative = "equivalence", mu = c(0.3, 0.7),
      R = 199)), "null value other than 0.5")
    expect_warning(suppressMessages(brunner_munzel(
      x, y, test_method = tm, mu = 0.7, R = 199)), "null value other than 0.5")
  }
})

test_that("no null-value warning for logit, perm_logit, or the standard null", {
  set.seed(3)
  x <- rnorm(20); y <- rnorm(20)
  for (tm in c("logit", "perm_logit")) {
    expect_no_warning(suppressMessages(brunner_munzel(
      x, y, test_method = tm, alternative = "minimal.effect", mu = c(0.3, 0.7),
      R = 199)))
    expect_no_warning(suppressMessages(brunner_munzel(
      x, y, test_method = tm, alternative = "equivalence", mu = c(0.3, 0.7),
      R = 199)))
  }
  for (tm in c("t", "perm")) {
    expect_no_warning(suppressMessages(brunner_munzel(
      x, y, test_method = tm, R = 199)))
  }
})
