# Implementing NPC Alpha Calibration for Permutation Equivalence Tests

## Background

This document describes how to implement the alpha calibration procedure recommended by Arboretti, Pesarin, & Salmaso (2021) for permutation-based equivalence testing. The calibration corrects for the conservatism of the naive (uncalibrated) TOST procedure, particularly at small sample sizes and narrow equivalence margins.

### Mapping Between Paper Terminology and TOSTER

| Paper Term | TOSTER Argument | Null Hypothesis | Alternative | Combination Rule |
|---|---|---|---|---|
| **IU-NPC** (Intersection-Union) | `alternative = "equivalence"` | Non-equivalence (NEq): d ≤ -ε_I OR d ≥ ε_S | Equivalence (Eq): -ε_I < d < ε_S | `max(p_low, p_high)` compared to α_c |
| **UI-NPC** (Union-Intersection) | `alternative = "minimal.effect"` | Equivalence (Eq): -ε_I ≤ d ≤ ε_S | Non-equivalence (NEq): d < -ε_I OR d > ε_S | `min(p_low, p_high)` compared to α̃_c |

Both directions are already implemented in `perm_t_test` and `brunner_munzel`. What is currently missing is the calibration of the critical value used to evaluate the global p-value. Both functions currently compare the global p-value to the nominal α (the "naive" approach). The calibrated approach replaces α with a data-informed critical value (α_c or α̃_c) that ensures the true global type I error rate equals the nominal α.

### Why Calibration Matters

**For `alternative = "equivalence"` (IU-NPC):**
The two partial tests (lower bound, upper bound) are negatively dependent. Using uncalibrated α for each partial test can make the global procedure dramatically conservative. The calibrated α_c lies in the interval [α, (1+α)/2). For α = 0.05, this means α_c ∈ [0.05, 0.525). At small n and narrow margins, the naive procedure can have maximum power near zero, meaning it would fail to detect that a treatment is equivalent to itself.

**For `alternative = "minimal.effect"` (UI-NPC):**
The calibrated α̃_c lies in [α/2, α]. For α = 0.05, this means α̃_c ∈ [0.025, 0.05]. The gap is much smaller, the convergence to α is faster, and the practical impact of not calibrating is less severe. Still, for small samples with tight margins, calibration improves accuracy.

## Implementation Plan

### Overview

The calibration is a Monte Carlo procedure that must be run before the main test. It estimates the critical value (α_c or α̃_c) such that the type I error rate at the boundary of the null hypothesis equals the nominal α. This is computationally expensive (a double loop: MC outer simulations × R inner permutations per simulation), so it should be implemented as:

1. A standalone helper function (`calibrate_alpha_perm`) that returns α_c or α̃_c
2. An optional `calibrate` argument on `perm_t_test` (and potentially `brunner_munzel`) that triggers calibration internally

### Algorithm

The algorithm below follows Section 4.2 of the paper's supplementary material. The key idea: simulate data at the boundary of H₀ (the worst case for type I error), run the permutation test at each simulation, and find the partial critical value that yields a global rejection rate of exactly α.

#### For `perm_t_test` (location-based tests)

**Step 0: Define the boundary.**

For the IU direction (`alternative = "equivalence"`), H₀ is true when d = -ε_I OR d = ε_S. The calibrated α_c must control type I error at both boundaries. In practice, for symmetric margins (ε_I = ε_S = ε), both boundaries yield the same α_c by symmetry, so only one boundary needs to be simulated.

For the UI direction (`alternative = "minimal.effect"`), H̃₀ is true when -ε_I ≤ d ≤ ε_S. The worst case for type I error is at the boundary d = -ε_I or d = ε_S (the edges of the equivalence region).

**Step 1: For a candidate critical value α_candidate, estimate the global rejection rate.**

```
function estimate_rejection_rate(n1, n2, eps_I, eps_S, alpha_candidate,
                                  direction, MC, R, sigma_hat,
                                  boundary = "upper"):

    rejections = 0

    for m in 1:MC:
        # Generate data at the boundary of H0
        if boundary == "upper":
            # d = eps_S: group 1 has mean 0, group 2 has mean -eps_S
            # (so d = mean_1 - mean_2 = eps_S, right at the boundary)
            x_sim = rnorm(n1, mean = 0, sd = sigma_hat)
            y_sim = rnorm(n2, mean = -eps_S, sd = sigma_hat)
        else:  # boundary == "lower"
            # d = -eps_I: group 1 has mean 0, group 2 has mean eps_I
            x_sim = rnorm(n1, mean = 0, sd = sigma_hat)
            y_sim = rnorm(n2, mean = eps_I, sd = sigma_hat)

        # Run the permutation test (the same perm_t_test internals)
        # but compare the global p-value to alpha_candidate instead of alpha
        result = perm_t_test(x_sim, y_sim,
                             alternative = direction,
                             mu = c(-eps_I, eps_S),
                             R = R)

        if direction == "equivalence":
            # IU: reject H0 (conclude Eq) when p.value <= alpha_candidate
            if result$p.value <= alpha_candidate:
                rejections += 1
        else:  # "minimal.effect"
            # UI: reject H̃0 (conclude NEq) when p.value <= alpha_candidate
            if result$p.value <= alpha_candidate:
                rejections += 1

    return rejections / MC
```

**Step 2: Search for the α_candidate that yields a rejection rate of α.**

Use a one-dimensional root-finding algorithm (e.g., `uniroot` in R) to solve:

```
estimate_rejection_rate(alpha_candidate) - alpha = 0
```

Search within the theoretical bounds:
- IU (`"equivalence"`): search α_candidate ∈ [α, (1 + α)/2)
- UI (`"minimal.effect"`): search α_candidate ∈ [α/2, α]

### Detailed R Pseudocode

```r
#' Calibrate Alpha for Permutation Equivalence/Minimal Effect Tests
#'
#' Estimates the calibrated critical value (alpha_c or alpha_tilde_c) via
#' Monte Carlo simulation at the boundary of the null hypothesis.
#'
#' @param n1 Sample size for group 1 (or total n for paired/one-sample)
#' @param n2 Sample size for group 2 (NULL for paired/one-sample)
#' @param eqb Equivalence bounds: a vector of length 2, c(lower, upper).
#'     For symmetric bounds, c(-eps, eps).
#' @param alpha Nominal significance level (default 0.05)
#' @param alternative Either "equivalence" (IU-NPC) or "minimal.effect" (UI-NPC)
#' @param paired Logical: paired or one-sample design?
#' @param var.equal Logical: assume equal variances? (passed to perm_t_test)
#' @param sigma_hat Estimated standard deviation of the data. If NULL, defaults to 1
#'     (equivalent to standardized margins). Can be estimated from pilot data or
#'     from the pooled SD of the observed data.
#' @param R Number of permutations per simulated test (default 2500, per paper)
#' @param MC Number of Monte Carlo simulation runs (default 5000, per paper)
#' @param tr Trimming fraction (passed to perm_t_test, default 0)
#' @param seed Random seed for reproducibility
#'
#' @return A list with:
#'   - alpha_c: the calibrated critical value
#'   - alpha_nominal: the input alpha
#'   - alternative: which direction was calibrated
#'   - search_interval: the theoretical bounds searched
#'   - rejection_rate_at_boundary: estimated type I error at the calibrated alpha_c
#'   - MC: number of Monte Carlo runs used
#'   - R: number of permutations per run
#'
#' @details
#' For alternative = "equivalence" (IU-NPC from Arboretti et al., 2021):
#'   - Searches alpha_c in [alpha, (1 + alpha)/2)
#'   - Data are simulated at d = eps_S (upper boundary of H0)
#'   - The calibrated alpha_c replaces the nominal alpha when evaluating the
#'     global p-value from perm_t_test(..., alternative = "equivalence")
#'
#' For alternative = "minimal.effect" (UI-NPC from Arboretti et al., 2021):
#'   - Searches alpha_tilde_c in [alpha/2, alpha]
#'   - Data are simulated at d = eps_S (boundary of H-tilde-0)
#'   - The calibrated alpha_tilde_c replaces the nominal alpha when evaluating
#'     the global p-value from perm_t_test(..., alternative = "minimal.effect")
#'
#' This function is computationally expensive (MC * R permutations total).
#' With defaults MC = 5000 and R = 2500, this is 12.5 million permutations.
#' Consider reducing MC and R for exploratory use, and increasing for final
#' analysis.
#'
#' @references
#' Arboretti, R., Pesarin, F., & Salmaso, L. (2021). A unified approach to
#' permutation testing for equivalence. Statistical Methods & Applications,
#' 30, 1033-1052.

calibrate_alpha_perm <- function(n1, n2 = NULL,
                                  eqb,
                                  alpha = 0.05,
                                  alternative = c("equivalence",
                                                  "minimal.effect"),
                                  paired = FALSE,
                                  var.equal = FALSE,
                                  sigma_hat = 1,
                                  R = 2500,
                                  MC = 5000,
                                  tr = 0,
                                  seed = NULL) {

  alternative <- match.arg(alternative)

  # Validate and parse bounds
  if (length(eqb) == 1) {
    eqb <- c(-abs(eqb), abs(eqb))
  }
  stopifnot(length(eqb) == 2)
  eqb <- sort(eqb)
  eps_I <- abs(eqb[1])  # lower margin magnitude
  eps_S <- eqb[2]       # upper margin

  # Define search interval
  if (alternative == "equivalence") {
    search_lower <- alpha
    search_upper <- (1 + alpha) / 2 - 1e-6  # open upper bound
  } else {
    search_lower <- alpha / 2
    search_upper <- alpha
  }

  if (!is.null(seed)) set.seed(seed)

  # Inner function: estimate rejection rate at a candidate alpha_c
  est_rejection <- function(alpha_c) {
    rejections <- 0L

    for (m in seq_len(MC)) {
      # Simulate data at the upper boundary: d = eps_S
      # x ~ N(0, sigma), y ~ N(-eps_S, sigma) => mean(x) - mean(y) ≈ eps_S
      if (paired || is.null(n2)) {
        # One-sample/paired: simulate differences with mean = eps_S
        n <- n1
        z_sim <- rnorm(n, mean = eps_S, sd = sigma_hat)
        result <- perm_t_test(x = z_sim,
                              alternative = alternative,
                              mu = eqb,
                              R = R,
                              tr = tr,
                              keep_perm = FALSE)
      } else {
        # Two-sample
        x_sim <- rnorm(n1, mean = 0, sd = sigma_hat)
        y_sim <- rnorm(n2, mean = -eps_S, sd = sigma_hat)
        result <- perm_t_test(x = x_sim, y = y_sim,
                              alternative = alternative,
                              mu = eqb,
                              var.equal = var.equal,
                              R = R,
                              tr = tr,
                              keep_perm = FALSE)
      }

      if (result$p.value <= alpha_c) {
        rejections <- rejections + 1L
      }
    }

    rejections / MC
  }

  # Use uniroot to find alpha_c where rejection_rate = alpha
  # The function to zero: f(alpha_c) = est_rejection(alpha_c) - alpha
  #
  # Note: est_rejection is stochastic, so uniroot may need a tolerance
  # compatible with MC precision. With MC = 5000, SE ≈ sqrt(0.05*0.95/5000) ≈ 0.003.

  root_result <- uniroot(
    f = function(ac) est_rejection(ac) - alpha,
    interval = c(search_lower, search_upper),
    tol = 0.005,      # tolerance compatible with MC noise
    maxiter = 20       # limit iterations given stochastic function
  )

  alpha_c <- root_result$root

  # Verify with a final estimate
  final_rate <- est_rejection(alpha_c)

  list(
    alpha_c = alpha_c,
    alpha_nominal = alpha,
    alternative = alternative,
    search_interval = c(search_lower, search_upper),
    rejection_rate_at_boundary = final_rate,
    MC = MC,
    R = R,
    eqb = eqb,
    sigma_hat = sigma_hat
  )
}
```

### Important Implementation Notes

#### 1. The sigma_hat problem

The calibration requires simulating data from a known distribution F. The paper (IU.3, Section 7.1) acknowledges that F is never fully known in practice, and recommends substituting the unknown σ with the pooled sample estimate σ̂. This introduces approximation, but the paper notes:

- For UI (minimal.effect): the approximation error is bounded by α/2 and is generally negligible in practice.
- For IU (equivalence): the approximation can be more consequential, particularly when n and margins are small.

In the implementation, `sigma_hat` should default to the pooled SD from the observed data when called internally. When exposed as a standalone function, let the user provide it.

#### 2. Computational cost

With the paper's recommended defaults (MC = 5000, R = 2500), each evaluation of the objective function requires 12.5 million permutation test computations. The root-finding algorithm typically needs 5–15 evaluations, so total cost is on the order of 60–190 million permutations. This is substantial. Practical considerations:

- **Progress reporting**: Add a progress indicator or message for each uniroot iteration.
- **Reduced defaults for exploration**: Consider MC = 1000, R = 1000 as a "fast" mode.
- **Caching**: If the user's data hasn't changed, the calibrated value can be reused.
- **Parallelization**: The MC loop is embarrassingly parallel. Consider `future.apply` or `parallel::mclapply` as optional speedups.

#### 3. Symmetry shortcut

When ε_I = ε_S (symmetric bounds), the calibration only needs to simulate at one boundary. The α_c found at d = ε_S will be the same as at d = -ε_I by symmetry. For asymmetric bounds, simulate at both boundaries and take the α_c that yields the maximum rejection rate (the conservative choice).

#### 4. uniroot with stochastic functions

`uniroot` expects a deterministic function. Since `est_rejection` is stochastic (Monte Carlo), there are a few options to handle this:

- **Option A (simple)**: Use large enough MC that noise is small relative to tolerance. With MC = 5000 and α ≈ 0.05, the SE of the rejection rate estimate is about 0.003, so `tol = 0.005` is reasonable.
- **Option B (robust)**: Replace `uniroot` with a grid search over the interval, evaluate each grid point, and interpolate or select the closest. This is less efficient but immune to the non-monotonicity caused by simulation noise.
- **Option C (stochastic approximation)**: Use the Robbins-Monro algorithm or similar. More complex to implement but theoretically sound.

Option A is recommended for initial implementation. If users report instability, upgrade to Option B.

#### 5. Mid-rank alternative

The paper (IU.3, UI.1) suggests that when no distributional assumption is tenable, mid-rank transformation of data and margins can provide reliable calibration. This corresponds to using the Brunner-Munzel test rather than the t-test. For `brunner_munzel`, the calibration would simulate from a distribution where the relative effect p equals the bound value (e.g., p = 0.3 or p = 0.7), which is more complex to set up because you need to construct two distributions with a specific probability of superiority. One approach: use normal distributions with a location shift chosen to achieve the target relative effect.

#### 6. When NOT to calibrate

The paper shows that calibrated and uncalibrated α converge when the standardized equivalence interval is large enough. Specifically (from Section 7.2, point IU.2), when:

(ε_I + ε_S) × sqrt(n1 × n2 / (n × σ²)) > ~5.4

the naive and calibrated values approximately coincide. The function could check this condition and skip calibration with an informative message when it's met.

### Integration with Existing Functions

#### Option A: Standalone function + manual workflow

The user runs calibration separately and passes the result:

```r
# Step 1: Calibrate (computationally expensive, do once)
cal <- calibrate_alpha_perm(
  n1 = 20, n2 = 20,
  eqb = c(-0.5, 0.5),
  sigma_hat = pooled_sd,
  alternative = "equivalence"
)

# Step 2: Run the test, comparing p-value to calibrated alpha
result <- perm_t_test(x, y, alternative = "equivalence", mu = c(-0.5, 0.5))

# Step 3: Evaluate using calibrated threshold
if (result$p.value <= cal$alpha_c) {
  # Reject H0 (conclude equivalence) at calibrated alpha
}
```

#### Option B: Integrated argument on perm_t_test

Add an optional `calibrate` argument:

```r
result <- perm_t_test(x, y,
                      alternative = "equivalence",
                      mu = c(-0.5, 0.5),
                      calibrate = TRUE,        # triggers calibration
                      cal_MC = 5000,            # MC runs for calibration
                      cal_R = 2500)             # permutations per calibration run
```

When `calibrate = TRUE`, the function:
1. Estimates σ̂ from the pooled data
2. Runs `calibrate_alpha_perm` internally
3. Compares the global p-value to α_c instead of α
4. Reports α_c in the output (e.g., as an additional list element)

#### Recommendation

Start with **Option A** (standalone function). This keeps the main test functions clean, gives users transparency into the calibration process, and avoids surprising users with long computation times. Option B can be added later as a convenience wrapper. The standalone function is also useful for power analysis and study design, independent of any particular dataset.

### Unit Tests for Calibration

#### Test 1: Search interval bounds are respected

```r
test_that("calibrated alpha_c is within theoretical bounds", {
  cal_iu <- calibrate_alpha_perm(n1 = 12, n2 = 12,
                                  eqb = c(-0.5, 0.5),
                                  alternative = "equivalence",
                                  MC = 500, R = 500, seed = 42)
  expect_gte(cal_iu$alpha_c, 0.05)
  expect_lt(cal_iu$alpha_c, (1 + 0.05) / 2)

  cal_ui <- calibrate_alpha_perm(n1 = 12, n2 = 12,
                                  eqb = c(-0.5, 0.5),
                                  alternative = "minimal.effect",
                                  MC = 500, R = 500, seed = 42)
  expect_gte(cal_ui$alpha_c, 0.05 / 2)
  expect_lte(cal_ui$alpha_c, 0.05)
})
```

#### Test 2: IU calibrated alpha_c > nominal alpha for small margins

```r
test_that("IU alpha_c exceeds nominal alpha for small margins and moderate n", {
  cal <- calibrate_alpha_perm(n1 = 12, n2 = 12,
                               eqb = c(-0.25, 0.25),
                               alternative = "equivalence",
                               MC = 1000, R = 1000, seed = 123)
  # With small margins relative to sigma, alpha_c should be notably above 0.05
  expect_gt(cal$alpha_c, 0.05 + 0.01)
})
```

#### Test 3: Large margins/samples yield alpha_c ≈ alpha

```r
test_that("alpha_c converges to alpha for large standardized margins", {
  cal <- calibrate_alpha_perm(n1 = 50, n2 = 50,
                               eqb = c(-1.0, 1.0),
                               alternative = "equivalence",
                               MC = 1000, R = 1000, seed = 42)
  expect_equal(cal$alpha_c, 0.05, tolerance = 0.015)
})
```

#### Test 4: UI alpha_tilde_c is close to alpha

```r
test_that("UI alpha_tilde_c is close to nominal alpha", {
  cal <- calibrate_alpha_perm(n1 = 15, n2 = 15,
                               eqb = c(-0.5, 0.5),
                               alternative = "minimal.effect",
                               MC = 1000, R = 1000, seed = 42)
  # UI calibrated value should be close to alpha (within [alpha/2, alpha])
  expect_equal(cal$alpha_c, 0.05, tolerance = 0.03)
})
```

## References

- Arboretti, R., Pesarin, F., & Salmaso, L. (2021). A unified approach to permutation testing for equivalence. *Statistical Methods & Applications*, 30, 1033–1052.
- Arboretti, R., Carrozzo, E., Pesarin, F., & Salmaso, L. (2018). Testing for equivalence: an intersection-union permutation solution. *Statistics in Biopharmaceutical Research*, 10, 130–138.
- Pesarin, F., & Salmaso, L. (2010). *Permutation Tests for Complex Data: Theory, Applications and Software*. Wiley.
- Wellek, S. (2010). *Testing Statistical Hypotheses of Equivalence and Noninferiority*. Chapman & Hall/CRC.
