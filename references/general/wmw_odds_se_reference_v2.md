# Reference: Updating Standard Error and CI Calculations for Non-Parametric Effect Sizes in TOSTER

## Overview

This document provides comprehensive guidance for updating the `ses_calc()` and related functions in the TOSTER package to compute proper standard errors and confidence intervals for all non-parametric effect sizes: rank-biserial correlation, concordance probability (c-statistic), WMW odds, and log-odds.

**Key recommendation:** Compute variance on the log-odds scale (which has the best asymptotic properties), then back-transform CIs to all other scales including the rank-biserial correlation.

---

## Part 1: Effect Size Definitions and Transformations

### Effect Size Measures

| Measure | Symbol | Definition | Range | Null Value |
|---------|--------|------------|-------|------------|
| Concordance probability | p | Pr(Y > X) + 0.5·Pr(Y = X) | [0, 1] | 0.5 |
| Rank-biserial correlation | r | 2p - 1 | [-1, 1] | 0 |
| WMW odds | α | p / (1 - p) | (0, ∞) | 1 |
| WMW log-odds | η | log(p / (1 - p)) = logit(p) | (-∞, ∞) | 0 |

### Transformation Functions (Bidirectional)

```r
# === Concordance (p) as the hub ===

# Concordance <-> Rank-biserial
cstat_to_rb <- function(p) 2 * p - 1
rb_to_cstat <- function(r) (r + 1) / 2

# Concordance <-> Odds
cstat_to_odds <- function(p) p / (1 - p)
odds_to_cstat <- function(alpha) alpha / (1 + alpha)

# Concordance <-> Log-odds
cstat_to_logodds <- function(p) log(p / (1 - p))  # = qlogis(p)
logodds_to_cstat <- function(eta) exp(eta) / (1 + exp(eta))  # = plogis(eta)

# === Direct transformations ===

# Rank-biserial <-> Odds
rb_to_odds <- function(r) (1 + r) / (1 - r)
odds_to_rb <- function(alpha) (alpha - 1) / (alpha + 1)

# Rank-biserial <-> Log-odds  
rb_to_logodds <- function(r) log((1 + r) / (1 - r))  # = atanh(r) * 2
logodds_to_rb <- function(eta) (exp(eta) - 1) / (exp(eta) + 1)  # = tanh(eta/2)

# Odds <-> Log-odds
odds_to_logodds <- function(alpha) log(alpha)
logodds_to_odds <- function(eta) exp(eta)
```

---

## Part 2: Two Independent Samples

### 2.1 Notation

- x: observations from group 1 (n₁ observations)
- y: observations from group 2 (n₂ observations)  
- n = n₁ + n₂: total sample size
- w₁ = n₁/n, w₂ = n₂/n: sample proportions

### 2.2 Point Estimates

The concordance probability is estimated by:

```r
# Using Mann-Whitney U
U <- sum(sapply(y, function(yj) sum(x < yj) + 0.5 * sum(x == yj)))
p_hat <- U / (n1 * n2)

# Equivalently, using placement values
W <- sapply(y, function(yj) mean(x < yj) + 0.5 * mean(x == yj))
p_hat <- mean(W)
```

### 2.3 Variance of Concordance Probability

From Bamber (1975) and Agresti (1980), using Lehmann (1975, eq. 2.21):

```
Var(p̂) = (1/n₁) · (P₂₂₁ - p²) + (1/n₂) · (P₂₁₁ - p²)
```

Where:
- **P₂₂₁** = Pr(X < Y₁, X < Y₂) where Y₁, Y₂ are independent draws from group 2
  - Interpretation: probability that a random X is less than TWO independent Y's
- **P₂₁₁** = Pr(X₁ < Y, X₂ < Y) where X₁, X₂ are independent draws from group 1
  - Interpretation: probability that TWO independent X's are both less than a random Y

### 2.4 Estimating the Variance Components

**Step 1: Compute placement values**

```r
# For each x_i: what proportion of y's are greater (plus half ties)?
# This estimates Pr(Y > x_i | X = x_i)
V <- sapply(x, function(xi) mean(y > xi) + 0.5 * mean(y == xi))

# For each y_j: what proportion of x's are less (plus half ties)?
# This estimates Pr(X < y_j | Y = y_j)
W <- sapply(y, function(yj) mean(x < yj) + 0.5 * mean(x == yj))
```

**Step 2: Estimate variance components**

```r
# P_221: Pr(X < Y1, X < Y2) - estimated by E[W^2] with finite-sample correction
# Unbiased estimator:
P_221 <- (n1 / (n1 - 1)) * (mean(W^2) - p_hat^2 / n1)
# Or simplified (slightly biased but simpler):
P_221_simple <- mean(W^2)

# P_211: Pr(X1 < Y, X2 < Y) - estimated by E[V^2] with finite-sample correction
# Note: V estimates Pr(Y > X), so V corresponds to concordance from Y's perspective
P_211 <- (n2 / (n2 - 1)) * (mean(V^2) - p_hat^2 / n2)
# Or simplified:
P_211_simple <- mean(V^2)
```

**Step 3: Compute variance of p̂**

```r
var_p <- (P_221 - p_hat^2) / n1 + (P_211 - p_hat^2) / n2
se_p <- sqrt(max(var_p, 0))
```

### 2.5 Standard Errors for All Effect Sizes (Two-Sample)

Using the delta method, we can derive SEs for all effect sizes from SE(p̂):

| Effect Size | SE Formula | Derivation |
|-------------|------------|------------|
| Concordance (p) | SE(p̂) | Direct |
| Rank-biserial (r) | 2 · SE(p̂) | r = 2p - 1, so dr/dp = 2 |
| Odds (α) | SE(p̂) / (1-p̂)² | α = p/(1-p), so dα/dp = 1/(1-p)² |
| Log-odds (η) | SE(p̂) / [p̂(1-p̂)] | η = logit(p), so dη/dp = 1/[p(1-p)] |

```r
se_cstat <- se_p
se_rb <- 2 * se_p
se_odds <- se_p / (1 - p_hat)^2
se_logodds <- se_p / (p_hat * (1 - p_hat))
```

### 2.6 Confidence Interval Construction (Two-Sample)

**Recommended approach: Compute CI on log-odds scale, back-transform to all others**

The log-odds (η = logit(p)) has the fastest convergence to normality (Agresti, 1980).

```r
# 1. Compute CI on log-odds scale
eta_hat <- cstat_to_logodds(p_hat)
se_eta <- se_p / (p_hat * (1 - p_hat))
z <- qnorm(1 - alpha/2)

ci_logodds <- eta_hat + c(-1, 1) * z * se_eta

# 2. Back-transform to all other scales
ci_odds <- exp(ci_logodds)
ci_cstat <- logodds_to_cstat(ci_logodds)  # = plogis(ci_logodds)
ci_rb <- cstat_to_rb(ci_cstat)            # = 2 * ci_cstat - 1
```

**Why this approach is superior:**
1. Log-odds is unbounded (-∞, ∞), avoiding boundary issues
2. Faster convergence to normality than other scales
3. CIs are guaranteed to respect bounds after back-transformation
4. Consistent methodology across all effect sizes

---

## Part 3: Paired Samples (One-Sample on Differences)

### 3.1 Setup

For paired data (x₁, y₁), (x₂, y₂), ..., (xₙ, yₙ), we analyze the differences:
- dᵢ = yᵢ - xᵢ - μ₀ (where μ₀ is typically 0)

The matched-pairs rank-biserial correlation is based on signed ranks of |dᵢ|.

### 3.2 Point Estimates (Paired)

Let:
- n = number of non-zero differences
- R⁺ = sum of ranks for positive differences  
- R⁻ = sum of ranks for negative differences
- T = min(R⁺, R⁻) = Wilcoxon signed-rank statistic
- S = R⁺ + R⁻ = n(n+1)/2 = total rank sum

The matched-pairs rank-biserial correlation:
```r
r_hat <- (R_plus - R_minus) / S
# Equivalently:
r_hat <- 1 - (4 * T) / (n * (n + 1))
```

The paired concordance probability:
```r
p_hat <- (r_hat + 1) / 2  # = R_plus / S
```

### 3.3 Variance for Paired Samples

**From Agresti (1980, eq. 4.4-4.5) for matched pairs:**

For matched pairs where qᵢⱼ denotes the proportion of pairs with first observation in category i and second in category j:

```
α' = (Σᵢ<ⱼ qᵢⱼ) / (Σᵢ>ⱼ qᵢⱼ)
```

The estimated asymptotic variance of α'√n for a random sample of n pairs is:

```
Var(α') · n = α'² · (1 - Σᵢ qᵢᵢ) / (Σᵢ<ⱼ qᵢⱼ)²
```

**For continuous paired differences (practical implementation):**

When working with signed ranks, we can estimate the variance using a different approach based on the signed-rank statistic's variance.

The variance of the Wilcoxon signed-rank statistic W = R⁺ under H₀ is:
```
Var(W) = n(n+1)(2n+1) / 24
```

With ties at t values with counts c₁, c₂, ..., cₜ:
```
Var(W) = [n(n+1)(2n+1) - Σⱼ cⱼ(cⱼ-1)(cⱼ+1)/2] / 24
```

**Deriving SE for paired concordance:**

Since p̂ = R⁺/S where S = n(n+1)/2:

```r
# Variance of W (= R_plus) with tie correction
tie_correction <- sum(sapply(table(abs(d[d != 0])), function(t) t*(t-1)*(t+1)/2))
var_W <- (n * (n+1) * (2*n+1) - tie_correction) / 24

# SE of concordance
S <- n * (n + 1) / 2
se_p_paired <- sqrt(var_W) / S
```

### 3.4 Standard Errors for All Effect Sizes (Paired)

```r
se_cstat_paired <- se_p_paired
se_rb_paired <- 2 * se_p_paired
se_odds_paired <- se_p_paired / (1 - p_hat)^2
se_logodds_paired <- se_p_paired / (p_hat * (1 - p_hat))
```

### 3.5 Confidence Interval Construction (Paired)

**Same approach as two-sample: compute on log-odds scale, back-transform**

```r
eta_hat <- cstat_to_logodds(p_hat)
se_eta <- se_p_paired / (p_hat * (1 - p_hat))
z <- qnorm(1 - alpha/2)

ci_logodds <- eta_hat + c(-1, 1) * z * se_eta

ci_odds <- exp(ci_logodds)
ci_cstat <- logodds_to_cstat(ci_logodds)
ci_rb <- cstat_to_rb(ci_cstat)
```

---

## Part 4: Summary of SE and CI Formulas

### 4.1 Standard Errors

| Effect Size | Two-Sample SE | Paired SE |
|-------------|---------------|-----------|
| Concordance (p) | √[(P₂₂₁-p²)/n₁ + (P₂₁₁-p²)/n₂] | √Var(W) / S |
| Rank-biserial (r) | 2 · SE(p) | 2 · SE(p) |
| Odds (α) | SE(p) / (1-p)² | SE(p) / (1-p)² |
| Log-odds (η) | SE(p) / [p(1-p)] | SE(p) / [p(1-p)] |

### 4.2 Confidence Intervals (All via log-odds)

For both two-sample and paired:

```r
# Step 1: Get SE of concordance (method differs by design)
se_p <- ...  # two-sample or paired formula

# Step 2: Transform to log-odds scale
eta_hat <- qlogis(p_hat)  # = log(p_hat / (1 - p_hat))
se_eta <- se_p / (p_hat * (1 - p_hat))

# Step 3: Form CI on log-odds scale
z <- qnorm(1 - alpha/2)
ci_eta <- eta_hat + c(-1, 1) * z * se_eta

# Step 4: Back-transform to desired scale
ci_logodds <- ci_eta
ci_odds <- exp(ci_eta)
ci_cstat <- plogis(ci_eta)
ci_rb <- 2 * plogis(ci_eta) - 1
```

---

## Part 5: Complete Implementation

### 5.1 Helper Functions

```r
#' Compute placement values for two-sample comparison
#' @param x numeric vector for group 1
#' @param y numeric vector for group 2
#' @return list with V (placements for x), W (placements for y), p_hat
compute_placements <- function(x, y) {
  n1 <- length(x)
  n2 <- length(y)
  
  # For each y: proportion of x's that are less (concordant pairs)
  W <- sapply(y, function(yj) {
    mean(x < yj) + 0.5 * mean(x == yj)
  })
  
  # For each x: proportion of y's that are greater
  V <- sapply(x, function(xi) {
    mean(y > xi) + 0.5 * mean(y == xi)
  })
  
  p_hat <- mean(W)
  
  list(V = V, W = W, p_hat = p_hat, n1 = n1, n2 = n2)
}

#' Compute variance of concordance probability (two-sample)
#' @param placements output from compute_placements
#' @return variance of p_hat
var_concordance_twosample <- function(placements) {
  with(placements, {
    # Variance components
    P_221 <- mean(W^2)
    P_211 <- mean(V^2)
    
    # Variance of p_hat (Lehmann/Agresti formula)
    var_p <- (P_221 - p_hat^2) / n1 + (P_211 - p_hat^2) / n2
    
    max(var_p, 0)  # Ensure non-negative
  })
}

#' Compute variance of concordance probability (paired)
#' @param d vector of differences (y - x)
#' @return variance of p_hat for paired design
var_concordance_paired <- function(d) {
  # Remove zeros
  d_nonzero <- d[d != 0]
  n <- length(d_nonzero)
  
  if (n < 2) return(NA)
  
  # Tie correction
  tie_counts <- table(abs(d_nonzero))
  tie_correction <- sum(tie_counts * (tie_counts - 1) * (tie_counts + 1) / 2)
  
  # Variance of W (Wilcoxon statistic)
  var_W <- (n * (n + 1) * (2 * n + 1) - tie_correction) / 24
  
  # Total rank sum
  S <- n * (n + 1) / 2
  
  # Variance of p_hat = W / S
  var_p <- var_W / S^2
  
  max(var_p, 0)
}
```

### 5.2 Main SE Calculation Function

```r
#' Compute SEs for all non-parametric effect sizes
#' @param x numeric vector (group 1 or pre-treatment)
#' @param y numeric vector (group 2 or post-treatment), NULL for one-sample
#' @param paired logical, TRUE for paired samples
#' @param mu hypothesized difference (default 0)
#' @return list with point estimates and SEs for all effect sizes
ses_compute <- function(x, y = NULL, paired = FALSE, mu = 0) {
  
  if (is.null(y)) {
    # One-sample: compare x to mu
    y <- rep(mu, length(x))
    paired <- TRUE
  }
  
  if (paired) {
    # === PAIRED SAMPLES ===
    d <- y - x - mu
    d_nonzero <- d[d != 0]
    n <- length(d_nonzero)
    
    if (n < 2) {
      warning("Fewer than 2 non-zero differences")
      return(NULL)
    }
    
    # Signed ranks
    ranks <- rank(abs(d_nonzero))
    R_plus <- sum(ranks[d_nonzero > 0])
    R_minus <- sum(ranks[d_nonzero < 0])
    S <- n * (n + 1) / 2
    
    # Point estimates
    p_hat <- R_plus / S
    r_hat <- (R_plus - R_minus) / S  # = 2 * p_hat - 1
    
    # Variance
    var_p <- var_concordance_paired(d)
    se_p <- sqrt(var_p)
    
  } else {
    # === TWO INDEPENDENT SAMPLES ===
    x <- na.omit(x)
    y <- na.omit(y)
    
    # Compute placements and variance
    placements <- compute_placements(x - mu, y)
    p_hat <- placements$p_hat
    r_hat <- 2 * p_hat - 1
    
    var_p <- var_concordance_twosample(placements)
    se_p <- sqrt(var_p)
  }
  
  # Point estimates for all scales
  alpha_hat <- p_hat / (1 - p_hat)
  eta_hat <- log(alpha_hat)  # = qlogis(p_hat)
  
  # Standard errors (delta method)
  se_cstat <- se_p
  se_rb <- 2 * se_p
  se_odds <- se_p / (1 - p_hat)^2
  se_logodds <- se_p / (p_hat * (1 - p_hat))
  
  list(
    # Point estimates
    cstat = p_hat,
    rb = r_hat,
    odds = alpha_hat,
    logodds = eta_hat,
    # Standard errors
    se_cstat = se_cstat,
    se_rb = se_rb,
    se_odds = se_odds,
    se_logodds = se_logodds,
    # Design info
    paired = paired
  )
}
```

### 5.3 CI Calculation Function

```r
#' Compute CIs for non-parametric effect sizes via log-odds transformation
#' @param ses_results output from ses_compute
#' @param conf.level confidence level (default 0.95)
#' @return list with CIs for all effect sizes
ses_ci <- function(ses_results, conf.level = 0.95) {
  
  alpha <- 1 - conf.level
  z <- qnorm(1 - alpha / 2)
  
  with(ses_results, {
    # CI on log-odds scale (best asymptotic properties)
    ci_logodds <- logodds + c(-1, 1) * z * se_logodds
    
    # Back-transform to all other scales
    ci_odds <- exp(ci_logodds)
    ci_cstat <- plogis(ci_logodds)
    ci_rb <- 2 * ci_cstat - 1
    
    list(
      ci_cstat = ci_cstat,
      ci_rb = ci_rb,
      ci_odds = ci_odds,
      ci_logodds = ci_logodds,
      conf.level = conf.level
    )
  })
}
```

### 5.4 Updated np_ses Function

```r
#' Non-parametric standardized effect sizes with proper variance estimation
#' @param x numeric vector (group 1 or pre-treatment)
#' @param y numeric vector (group 2 or post-treatment), NULL for one-sample
#' @param mu hypothesized value/difference (default 0)
#' @param conf.level confidence level (default 0.95)
#' @param paired logical for paired samples
#' @param ses effect size type: "rb", "cstat", "odds", or "logodds"
#' @return list with effect size estimate, CI, and metadata
np_ses_updated <- function(x, y = NULL, mu = 0, conf.level = 0.95,
                           paired = FALSE, ses = c("rb", "odds", "logodds", "cstat")) {
  
  ses <- match.arg(ses)
  
  # Compute all estimates and SEs
  est <- ses_compute(x = x, y = y, paired = paired, mu = mu)
  
  if (is.null(est)) {
    return(list(type = ses, est = NA, conf.int = c(NA, NA), 
                se = NA, paired = paired, mu = mu))
  }
  
  # Compute CIs via log-odds transformation
  cis <- ses_ci(est, conf.level = conf.level)
  
  # Extract requested effect size
  result <- switch(ses,
    "rb" = list(
      est = est$rb,
      se = est$se_rb,
      conf.int = cis$ci_rb
    ),
    "cstat" = list(
      est = est$cstat,
      se = est$se_cstat,
      conf.int = cis$ci_cstat
    ),
    "odds" = list(
      est = est$odds,
      se = est$se_odds,
      conf.int = cis$ci_odds
    ),
    "logodds" = list(
      est = est$logodds,
      se = est$se_logodds,
      conf.int = cis$ci_logodds
    )
  )
  
  list(
    type = ses,
    est = result$est,
    se = result$se,
    conf.int = result$conf.int,
    conf.level = conf.level,
    paired = paired,
    mu = mu
  )
}
```

---

## Part 6: Comparison with Current Implementation

### Current Approach (rbs.R)

The current implementation:
1. Computes SE for rank-biserial using a formula derived from the Wilcoxon statistic variance
2. Uses Fisher's z-transformation (atanh) for CIs on rank-biserial
3. Simply transforms the CI bounds for odds/logodds/cstat

```r
# Current: Fisher z-transform for rank-biserial
rf <- atanh(r_rbs)
rfSE <- ...  # SE on Fisher z scale
confint <- tanh(rf + c(-1, 1) * qnorm(1 - alpha/2) * rfSE)

# Current: Simple transformation of CI bounds
confint_odds <- rb_to_odds(confint)  # NOT proper propagation
```

### New Approach

1. Computes SE for concordance probability using Lehmann/Agresti variance formula
2. Uses log-odds transformation for CIs (better asymptotic properties than Fisher z)
3. Back-transforms CIs to all scales including rank-biserial

```r
# New: Log-odds transform for all effect sizes
eta <- qlogis(p_hat)
se_eta <- se_p / (p_hat * (1 - p_hat))
ci_eta <- eta + c(-1, 1) * qnorm(1 - alpha/2) * se_eta

# New: Proper back-transformation
ci_rb <- 2 * plogis(ci_eta) - 1  # Proper CI for rank-biserial
ci_odds <- exp(ci_eta)           # Proper CI for odds
```

### Key Differences

| Aspect | Current | New |
|--------|---------|-----|
| Primary scale | Rank-biserial | Log-odds |
| Transformation | Fisher z (atanh) | Logit |
| Variance basis | Wilcoxon W variance | Concordance probability variance |
| CI for odds | Transform bounds | Proper back-transform |
| Asymptotic properties | Good | Better (faster convergence) |

---

## Part 7: Edge Cases and Robustness

### 7.1 Boundary Cases

When p̂ is near 0 or 1:
- Log-odds approaches ±∞
- SE formula SE(η) = SE(p) / [p(1-p)] can become very large
- CIs will be wide but remain valid after back-transformation

```r
# Safeguard for extreme values
p_hat_safe <- pmin(pmax(p_hat, 1e-10), 1 - 1e-10)
```

### 7.2 Small Samples

For very small samples:
- Variance estimates may be unstable
- Consider using exact methods or bootstrap as alternatives
- Add warnings when n₁ < 5 or n₂ < 5 (two-sample) or n < 5 (paired)

### 7.3 Heavy Ties

Ties are handled through:
- Half-credit in placement values (0.5 for each tie)
- Tie correction in paired variance formula

The variance formulas remain valid with ties, though coverage may be slightly affected with very heavy tying.

---

## Part 8: Testing Recommendations

### 8.1 Verification Tests

1. **Point estimates unchanged**: All point estimates should match the current implementation exactly.

2. **Known values**: Test against published examples (e.g., Agresti 1980 examples).

3. **Transformation consistency**: Verify that:
   ```r
   all.equal(2 * ci_cstat - 1, ci_rb)
   all.equal(exp(ci_logodds), ci_odds)
   all.equal(plogis(ci_logodds), ci_cstat)
   ```

### 8.2 Coverage Simulations

Run simulations varying:
- Sample sizes: (10,10), (20,20), (50,50), (10,30), (30,10)
- True effect sizes: p = 0.5, 0.6, 0.7, 0.8, 0.9
- Distributions: Normal, Exponential, Heavy-tailed

Target: 95% CIs should achieve close to 95% coverage.

### 8.3 Comparison Tests

Compare new CIs to:
- Current TOSTER implementation
- Bootstrap percentile CIs (as reference)
- Other software (e.g., SAS PROC NPAR1WAY)

---

## References

1. **Agresti, A. (1980).** Generalized odds ratios for ordinal data. *Biometrics*, 36, 59-67.

2. **Bamber, D. (1975).** The area above the ordinal dominance graph and the area below the receiver operating characteristic graph. *Journal of Mathematical Psychology*, 12, 387-415.

3. **O'Brien, R.G. & Castelloe, J.M. (2006).** Exploiting the link between the Wilcoxon-Mann-Whitney test and a simple odds statistic. *SUGI 31 Proceedings*, Paper 209-31.

4. **Lehmann, E.L. (1975).** *Nonparametrics: Statistical Methods Based on Ranks*. Holden-Day.

5. **Kerby, D.S. (2014).** The simple difference formula: An approach to teaching nonparametric correlation. *Comprehensive Psychology*, 3, 11-IT.

6. **Divine, G.W., Norton, H.J., Barón, A.E., & Juarez-Colunga, E. (2018).** The Wilcoxon-Mann-Whitney procedure fails as a test of medians. *The American Statistician*, 72(3), 278-286.
