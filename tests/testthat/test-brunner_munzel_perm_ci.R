# Permutation CI / p-value duality for brunner_munzel (#120) ----

# Re-run at a new null value with the same permutations (same seed)
bm_run <- function(args, mu, seed = 1) {
  args$mu <- mu
  args$test_method <- "perm"
  set.seed(seed)
  suppressWarnings(suppressMessages(do.call(brunner_munzel, args)))
}

# mu just inside the CI gives p > alpha; just outside gives p <= alpha
# (endpoints clamped to 0 or 1 are skipped)
expect_bm_duality <- function(args, seed = 1) {
  alpha <- if (is.null(args$alpha)) 0.05 else args$alpha
  res <- bm_run(args, 0.5, seed)
  ci <- res$conf.int
  eps <- 1e-7
  checked <- 0

  if (ci[1] > 1e-6) {
    expect_gt(bm_run(args, ci[1] + eps, seed)$p.value, alpha)
    expect_lte(bm_run(args, ci[1] - eps, seed)$p.value, alpha)
    checked <- checked + 1
  }
  if (ci[2] < 1 - 1e-6) {
    expect_gt(bm_run(args, ci[2] - eps, seed)$p.value, alpha)
    expect_lte(bm_run(args, ci[2] + eps, seed)$p.value, alpha)
    checked <- checked + 1
  }
  expect_gt(checked, 0)
  invisible(res)
}

set.seed(1)
bm_x <- rnorm(8, mean = 0, sd = 2)
bm_y <- rnorm(20, mean = 0, sd = 1)
bm_x_small <- rnorm(6, mean = 0, sd = 2)
bm_y_small <- rnorm(8, mean = 0, sd = 1)
bm_px <- rnorm(20)
bm_py <- bm_px + rnorm(20, -0.2, 1.5)

test_that("BM permutation CI agrees with p-value: two-sample randomized", {
  skip_on_cran()
  for (alt in c("two.sided", "less", "greater")) {
    expect_bm_duality(list(x = bm_x, y = bm_y, alternative = alt, R = 999))
    expect_bm_duality(list(x = bm_x, y = bm_y, alternative = alt, R = 999,
                           p_method = "exact"))
  }
})

test_that("BM permutation CI agrees with p-value: two-sample exact enumeration", {
  skip_on_cran()
  # choose(14, 6) = 3003 permutations
  for (alt in c("two.sided", "less", "greater")) {
    for (pm in c("exact", "plusone")) {
      expect_bm_duality(list(x = bm_x_small, y = bm_y_small, alternative = alt,
                             R = 5000, p_method = pm))
    }
  }
})

test_that("BM permutation CI agrees with p-value: paired", {
  skip_on_cran()
  for (alt in c("two.sided", "less", "greater")) {
    # n > 13: randomized
    expect_bm_duality(list(x = bm_px, y = bm_py, paired = TRUE,
                           alternative = alt, R = 999))
    # n <= 13: exact enumeration of 2^12 swaps
    expect_bm_duality(list(x = bm_px[1:12], y = bm_py[1:12], paired = TRUE,
                           alternative = alt, p_method = "plusone"))
  }
})

test_that("BM permutation equivalence/MET decisions agree with the 1 - 2*alpha CI", {
  skip_on_cran()
  eps <- 1e-7
  for (paired in c(FALSE, TRUE)) {
    args <- if (paired) {
      list(x = bm_px, y = bm_py, paired = TRUE, R = 999)
    } else {
      list(x = bm_x, y = bm_y, R = 999)
    }
    run <- function(alt, mu) {
      args$alternative <- alt
      bm_run(args, mu)
    }
    ci <- run("equivalence", c(0.05, 0.95))$conf.int
    expect_equal(attr(ci, "conf.level"), 0.90)
    expect_true(ci[1] > 0.15 && ci[2] < 0.85)
    expect_equal(run("minimal.effect", c(0.05, 0.95))$conf.int, ci)

    # Equivalence: significant iff CI lies within the bounds
    expect_lte(run("equivalence", c(ci[1] - eps, ci[2] + eps))$p.value, 0.05)
    expect_gt(run("equivalence", c(ci[1] + eps, ci[2] + eps))$p.value, 0.05)
    expect_gt(run("equivalence", c(ci[1] - eps, ci[2] - eps))$p.value, 0.05)

    # Minimal effect: significant iff CI lies entirely outside the bounds
    expect_lte(run("minimal.effect", c(ci[2] + eps, ci[2] + 0.1))$p.value, 0.05)
    expect_gt(run("minimal.effect", c(ci[2] - eps, ci[2] + 0.1))$p.value, 0.05)
    expect_lte(run("minimal.effect", c(ci[1] - 0.1, ci[1] - eps))$p.value, 0.05)
    expect_gt(run("minimal.effect", c(ci[1] - 0.1, ci[1] + eps))$p.value, 0.05)
  }
})

test_that("BM two-sided permutation CI is centered on the estimate when unclamped", {
  skip_on_cran()
  res <- bm_run(list(x = bm_x, y = bm_y, R = 999), 0.5)
  ci <- res$conf.int
  skip_if(ci[1] <= 0 || ci[2] >= 1)
  est <- unname(res$estimate)
  expect_equal(est - ci[1], ci[2] - est)
})
