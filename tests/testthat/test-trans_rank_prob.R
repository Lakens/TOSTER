# test-trans_rank_prob.R

# Identity --------
test_that("trans_rank_prob returns inputs unchanged when from == to", {
  for (s in c("probability", "difference", "logodds", "odds")) {
    out <- trans_rank_prob(0.7, se = 0.05, ci = c(0.6, 0.8),
                           null = 0.5, from = s, to = s)
    expect_equal(out$estimate, 0.7)
    expect_equal(out$se, 0.05)
    expect_equal(out$ci, c(0.6, 0.8))
    expect_equal(out$null, 0.5)
  }
})

# Point estimates --------
test_that("trans_rank_prob computes correct point estimates from probability", {
  p <- 0.7
  expect_equal(trans_rank_prob(p, from = "probability", to = "difference")$estimate,
               2 * p - 1)
  expect_equal(trans_rank_prob(p, from = "probability", to = "odds")$estimate,
               p / (1 - p))
  expect_equal(trans_rank_prob(p, from = "probability", to = "logodds")$estimate,
               log(p / (1 - p)))
})

test_that("trans_rank_prob computes correct point estimates from difference", {
  d <- 0.4  # corresponds to p = 0.7
  expect_equal(trans_rank_prob(d, from = "difference", to = "probability")$estimate,
               0.7)
  expect_equal(trans_rank_prob(d, from = "difference", to = "odds")$estimate,
               0.7 / 0.3, tolerance = 1e-10)
})

test_that("trans_rank_prob computes correct point estimates from logodds", {
  eta <- qlogis(0.7)
  expect_equal(trans_rank_prob(eta, from = "logodds", to = "probability")$estimate,
               0.7, tolerance = 1e-10)
  expect_equal(trans_rank_prob(eta, from = "logodds", to = "difference")$estimate,
               0.4, tolerance = 1e-10)
})

test_that("trans_rank_prob computes correct point estimates from odds", {
  alpha <- 0.7 / 0.3
  expect_equal(trans_rank_prob(alpha, from = "odds", to = "probability")$estimate,
               0.7, tolerance = 1e-10)
})

# Round-trip SE --------
test_that("trans_rank_prob SE round-trips through all scales", {
  for (s in c("difference", "logodds", "odds")) {
    out <- trans_rank_prob(0.7, se = 0.05, from = "probability", to = s)
    back <- trans_rank_prob(out$estimate, se = out$se, from = s, to = "probability")
    expect_equal(back$se, 0.05, tolerance = 1e-10)
  }
})

# Round-trip CI --------
test_that("trans_rank_prob CI round-trips through all scales", {
  for (s in c("difference", "logodds", "odds")) {
    out <- trans_rank_prob(0.7, ci = c(0.6, 0.8), from = "probability", to = s)
    back <- trans_rank_prob(out$estimate, ci = out$ci, from = s, to = "probability")
    expect_equal(back$ci, c(0.6, 0.8), tolerance = 1e-10)
  }
})

# Cross-scale without going through probability --------
test_that("trans_rank_prob handles non-probability from/to pairs", {
  # difference -> logodds (should route through probability internally)
  d <- 0.4
  out <- trans_rank_prob(d, se = 0.1, ci = c(0.2, 0.6),
                         from = "difference", to = "logodds")
  # Verify against manual: d=0.4 -> p=0.7 -> logodds = qlogis(0.7)
  expect_equal(out$estimate, qlogis(0.7), tolerance = 1e-10)

  # odds -> difference
  alpha <- 0.7 / 0.3
  out2 <- trans_rank_prob(alpha, from = "odds", to = "difference")
  expect_equal(out2$estimate, 0.4, tolerance = 1e-10)
})

# Null values --------
test_that("trans_rank_prob transforms null values", {
  out <- trans_rank_prob(0.7, null = 0.5, from = "probability", to = "difference")
  expect_equal(out$null, 0)

  out2 <- trans_rank_prob(0.7, null = 0.5, from = "probability", to = "logodds")
  expect_equal(out2$null, 0)

  out3 <- trans_rank_prob(0.7, null = 0.5, from = "probability", to = "odds")
  expect_equal(out3$null, 1)
})

test_that("trans_rank_prob transforms equivalence bounds", {
  out <- trans_rank_prob(0.7, null = c(0.35, 0.65),
                         from = "probability", to = "difference")
  expect_equal(out$null, c(-0.3, 0.3))
})

# NULL handling --------
test_that("trans_rank_prob handles NULL optional arguments", {
  out <- trans_rank_prob(0.7, from = "probability", to = "difference")
  expect_null(out$se)
  expect_null(out$ci)
  expect_null(out$null)
  expect_equal(out$estimate, 0.4)
})

# Boundary behavior --------
test_that("trans_rank_prob handles boundary probability values", {
  out <- trans_rank_prob(0, from = "probability", to = "logodds")
  expect_equal(out$estimate, -Inf)

  out2 <- trans_rank_prob(1, from = "probability", to = "logodds")
  expect_equal(out2$estimate, Inf)

  out3 <- trans_rank_prob(1, from = "probability", to = "odds")
  expect_equal(out3$estimate, Inf)
})

# Consistency with ses_calc internals --------
test_that("trans_rank_prob matches ses_calc transformation functions", {
  # rb = 2*cstat - 1 (i.e., difference = 2*probability - 1)
  p <- 0.7
  expect_equal(
    trans_rank_prob(p, from = "probability", to = "difference")$estimate,
    2 * p - 1
  )
  # odds = p / (1-p)
  expect_equal(
    trans_rank_prob(p, from = "probability", to = "odds")$estimate,
    p / (1 - p)
  )
  # logodds = log(p / (1-p))
  expect_equal(
    trans_rank_prob(p, from = "probability", to = "logodds")$estimate,
    log(p / (1 - p))
  )
  # SE for difference = 2 * SE_p
  se_p <- 0.05
  expect_equal(
    trans_rank_prob(p, se = se_p, from = "probability", to = "difference")$se,
    2 * se_p
  )
  # SE for logodds = SE_p / (p * (1-p))
  expect_equal(
    trans_rank_prob(p, se = se_p, from = "probability", to = "logodds")$se,
    se_p / (p * (1 - p))
  )
  # SE for odds = SE_p / (1-p)^2
  expect_equal(
    trans_rank_prob(p, se = se_p, from = "probability", to = "odds")$se,
    se_p / (1 - p)^2
  )
})
