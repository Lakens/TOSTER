# Wilcoxon Rank-Based Tests: What They Actually Test

## Overview

This reference clarifies what the Wilcoxon-Mann-Whitney (two-sample), Wilcoxon signed-rank (paired/one-sample), and sign test procedures actually test, the assumptions required for various interpretations, and how to extend these to equivalence testing frameworks.

**Key insight:** The perception that the WMW procedure tests equality of medians is pervasive and frequently encountered. Unfortunately, this perception is mostly wrong (Divine et al., 2018).

---

## The Two-Sample Case: Wilcoxon-Mann-Whitney (WMW)

### The Test Statistic

The Mann-Whitney U statistic counts the number of pairs (Xᵢ, Yⱼ) where Xᵢ > Yⱼ. Equivalently, it can be expressed via ranks (Wilcoxon rank-sum formulation).

### What It Actually Tests (Distribution-Free)

**Without any assumptions**, the WMW test is a test of:

$$H_0: p = \Pr(X_1 < X_2) + \Pr(X_1 = X_2)/2 = 0.5$$

where X₁ and X₂ are random observations from the two groups being compared (Divine et al., 2018; O'Brien & Castelloe, 2006).

The sample estimate is:

$$\hat{p} = \frac{U}{n_1 n_2}$$

where U is the Mann-Whitney U statistic.

**This formulation:**
- Has nothing directly to do with means, medians, or even the shapes of the distributions
- Is valid for tied data (contrary to some textbook claims)
- Is the only interpretation that holds without additional assumptions

### What It Tests (Assumption-Dependent)

| Assumption Level | Null Hypothesis | Parameter Tested |
|------------------|-----------------|------------------|
| **None (distribution-free)** | p = Pr(X < Y) + Pr(X = Y)/2 = 0.5 | Stochastic equality |
| **Identical shape, different location (shift alternative)** | Δ = 0 | Hodges-Lehmann pseudomedian of pairwise differences |
| **Shift alternative + symmetric distributions** | median(X) - median(Y) = 0 | Difference in medians |
| **Shift alternative + symmetric distributions** | mean(X) - mean(Y) = 0 | Difference in means |

### The WMWodds: An Interpretable Effect Size

O'Brien and Castelloe (2006) recommend converting p̂ to an odds measure:

$$\text{WMWodds} = \frac{\hat{p}}{1 - \hat{p}}$$

**Interpretation:** If WMWodds = 2.0, the odds are 2:1 that a randomly selected observation from group 1 is less than a randomly selected observation from group 2 (splitting ties evenly).

**Advantages:**
- Null value of 1.0 is intuitive
- Comparable across studies regardless of scale
- Can be displayed in forest plots like odds ratios
- Confidence intervals available via Agresti's (1980) generalized odds ratio formulas

### Connection to ROC Curve Area (c-statistic)

The quantity p̂ = Pr(X < Y) + Pr(X = Y)/2 equals the area under the ROC curve (AUC) when viewing the outcome as a classifier for group membership (Bamber, 1975; Hanley & McNeil, 1982). This provides:

- An intuitive interpretation: the probability that a randomly selected observation from group 2 exceeds one from group 1
- A graphical representation via ROC curves
- Connection to discrimination/classification frameworks

### The Hodges-Lehmann Estimator

**Only under the location-shift model** F_X(x) = F_Y(x - Δ) does the Hodges-Lehmann estimator have meaning:

$$\hat{\Delta} = \text{median}\{X_i - Y_j : \text{all } i, j\}$$

This is the median of all n₁ × n₂ pairwise differences. Importantly:

- It estimates the **pseudomedian** of the distribution of X - Y
- The pseudomedian equals the median only for symmetric distributions
- The pseudomedian equals the mean only for symmetric distributions

### Counterexamples: Why WMW Fails as a Test of Medians

Divine et al. (2018) provide compelling counterexamples:

**1. Equal medians, significant WMW test:**
In an aromatherapy trial (Hunt et al., 2013), both groups had median PON scores of -1, yet WMW p < 0.001.

**2. Very different medians, non-significant WMW test:**
Constructed samples with medians of 9 vs. 99, but p̂ = 0.502, WMWodds = 1.01, p ≈ 1.0.

**3. Medians in the wrong direction:**
Samples where median A = 99 >> median B = 9, but p̂ = 0.716 with p = 0.046, indicating observations from A tend to be *lower* than B.

**4. Global intransitivity:**
It's possible to have A < B < C < A by WMW tests, which is impossible for any measure of central tendency.

### Validity with Tied Data

Despite claims in some textbooks, the WMW test **does not require continuous data** to be valid. Lehmann (1975) established the asymptotic normality of the WMW test statistic for tied data, with only a mild condition that no single point accounts for nearly all the probability.

When ties are present:
- Tied observations receive their average rank
- The variance estimator is adjusted for ties
- The null hypothesis p = 0.5 remains valid

### The Behrens-Fisher Problem

When variances are unequal between groups:
- The Brunner-Munzel (2000) variation performs well when minimum n ≥ 30
- The Fligner-Policello (1981) variation assumes continuous data
- For smaller samples or many ties, exact/permutation tests are recommended

### R Implementation

```r
# Basic test (tests p = 0.5)
wilcox.test(x = x, y = y)

# With location shift (tests H₀: Δ = μ₀, requires shift assumption)
wilcox.test(x = x, y = y, mu = 0)

# One-sided with confidence interval
wilcox.test(x = x, y = y, mu = 0, alternative = "greater", conf.int = TRUE)

# Compute WMWodds manually
U <- wilcox.test(x, y)$statistic
n1 <- length(x)
n2 <- length(y)
p_hat <- U / (n1 * n2)
wmw_odds <- p_hat / (1 - p_hat)
```

---

## The Paired/One-Sample Case: Wilcoxon Signed-Rank Test

### The Test Statistic

For differences Dᵢ = Xᵢ - Yᵢ (paired) or observations Xᵢ (one-sample), the test:

1. Computes absolute values |Dᵢ|
2. Ranks the absolute values
3. Sums ranks of positive differences (W⁺) and negative differences (W⁻)

### What It Tests (Assumption-Dependent)

| Assumption Level | Null Hypothesis | Parameter Tested |
|------------------|-----------------|------------------|
| **Symmetric distribution around θ** | θ = 0 | Pseudomedian of differences |
| **Symmetric distribution** | median(D) = 0 | Median of differences |
| **Symmetric distribution** | mean(D) = 0 | Mean of differences |

### Critical Assumption: Symmetry

The signed-rank test **requires** the assumption that the distribution of differences is symmetric around some value θ. Without this assumption, the test does not have a clear interpretation.

**Note:** The signed-rank test has a similarly poor connection to sample medians as the WMW test (Divine et al., 2018). The quantity it tests, Pr[X₁ + X₂ < 0], is not as interpretable as Pr(X < Y) is for the WMW test.

### The Hodges-Lehmann Estimator (One-Sample/Paired)

For the signed-rank setting:

$$\hat{\theta} = \text{median}\left\{\frac{D_i + D_j}{2} : i \leq j\right\}$$

This is the median of all n(n+1)/2 Walsh averages (pairwise averages including self-pairs).

### R Implementation

```r
# Paired test
wilcox.test(x = x, y = y, paired = TRUE)

# One-sample test against μ₀
wilcox.test(x = x, mu = 0)

# With confidence interval
wilcox.test(x = x, y = y, paired = TRUE, conf.int = TRUE)
```

---

## The Sign Test

### When to Use

When the symmetry assumption for the signed-rank test is untenable, the sign test provides a valid alternative.

### What It Tests

- **H₀:** median(D) = 0 (or P(D > 0) = 0.5)
- Uses only the signs of differences, ignoring magnitudes
- Valid without symmetry assumption
- Less powerful than signed-rank when symmetry holds

### R Implementation

```r
# Using binom.test on signs
positive <- sum(d > 0)
n_nonzero <- sum(d != 0)
binom.test(positive, n_nonzero)
```

---

## Equivalence and Non-Inferiority Testing with TOST

### Framework

Two One-Sided Tests (TOST) can be applied to rank-based tests by specifying equivalence bounds on the Hodges-Lehmann scale.

### Non-Inferiority (One-Sided)

For a non-inferiority margin of δ (e.g., X is non-inferior to Y if Δ > -δ):

```r
# H₀: Δ ≤ -δ vs H₁: Δ > -δ
wilcox.test(x = x, y = y, mu = -delta, alternative = "greater", conf.int = TRUE)
```

**Interpretation:** If p < α, conclude that X is non-inferior to Y on the Hodges-Lehmann scale, with margin δ.

### Equivalence (Two One-Sided Tests)

For equivalence bounds (-δ, δ):

```r
# Test 1: H₀: Δ ≤ -δ vs H₁: Δ > -δ
test_lower <- wilcox.test(x = x, y = y, mu = -delta, alternative = "greater")

# Test 2: H₀: Δ ≥ δ vs H₁: Δ < δ
test_upper <- wilcox.test(x = x, y = y, mu = delta, alternative = "less")

# Equivalence p-value
p_equiv <- max(test_lower$p.value, test_upper$p.value)
```

**Interpretation:** If p_equiv < α, conclude that the Hodges-Lehmann location shift lies within (-δ, δ).

### Using TOSTER Package

```r
library(TOSTER)

# Equivalence test for two independent samples
wilcox_TOST(x = x, y = y, low_eqbound = -delta, high_eqbound = delta)

# Paired equivalence test
wilcox_TOST(x = x, y = y, paired = TRUE, low_eqbound = -delta, high_eqbound = delta)
```

### Concordance of P-Values and Confidence Intervals

For TOST on the Hodges-Lehmann scale:

- The (1 - 2α) confidence interval should be contained within (-δ, δ) if and only if p_equiv < α
- Minor discrepancies at boundaries may occur due to the discrete nature of rank statistics

**Caution:** Divine et al. (2018) note that the Hodges-Lehmann confidence interval may perform poorly when the location-shift assumption is violated. In the aromatherapy example, the exact Hodges-Lehmann CI was [-1, 0], which seems inconsistent with p < 0.001.

---

## Graphical Representations of p̂

Divine et al. (2018) describe several ways to visualize p̂ = Pr(X < Y) + Pr(X = Y)/2:

### Bubble Plot
Plot all (X, Y) pairs with bubble size proportional to frequency. The proportion of bubble area above the identity line equals p̂.

### ROC Curve
The area under the ROC curve (treating group as "disease status") equals p̂. Useful when students are familiar with diagnostic testing.

### Dominance Diagram
A grid displaying the direction of difference for all combinations of ordered X and Y values (Newson, 2002). The proportion of "X < Y" cells plus half the tied cells equals p̂.

---

## Reporting Guidelines

### Methods Section Language

**Two-sample (distribution-free framing):**

> "Groups were compared using the Wilcoxon-Mann-Whitney test, which tests the null hypothesis that a randomly selected observation from group A is equally likely to be greater or less than a randomly selected observation from group B (i.e., Pr(X < Y) + Pr(X = Y)/2 = 0.5). The WMWodds effect size and its 95% CI are reported."

**Two-sample (location-shift framing, if justified):**

> "Non-inferiority was assessed using the Wilcoxon-Mann-Whitney test with a margin of [δ] units on the Hodges-Lehmann scale. Under the assumption of a pure location shift between groups, this tests whether the pseudomedian of pairwise differences exceeds -[δ]."

**Paired samples:**

> "Equivalence was assessed using TOST with Wilcoxon signed-rank tests and bounds of ±[δ] units. Under the assumption of symmetric differences, this tests whether the pseudomedian of paired differences lies within the equivalence region."

### What to Report

1. The p̂ estimate (or equivalently, WMWodds)
2. The Hodges-Lehmann point estimate (if location-shift assumption is reasonable)
3. The confidence interval (specify confidence level)
4. For equivalence/non-inferiority: the margin with justification
5. Clear statement of which assumptions are being invoked

---

## Summary Table: Assumptions and Interpretations

| Test | Minimum Assumption | Parameter | Additional Assumption for Median/Mean |
|------|-------------------|-----------|--------------------------------------|
| WMW (two-sample) | None | p = Pr(X < Y) + Pr(X = Y)/2 | Location shift + symmetry |
| WMW with Hodges-Lehmann | Location shift | Pseudomedian of pairwise differences | Symmetry for median/mean |
| Signed-rank (paired/one-sample) | Symmetric differences | Pseudomedian of D | Automatic (symmetry assumed) |
| Sign test | None | Median of D | Not applicable |

---

## Key Takeaways

1. **Empirically, the WMW test should be regarded as a test of p = Pr(X < Y) + Pr(X = Y)/2 = 0.5** (Divine et al., 2018). This is the only interpretation that holds without additional assumptions.

2. **The WMW procedure is in no way a function of the observed sample medians.** Counterexamples demonstrate significant tests with equal medians, non-significant tests with different medians, and even significant results in the direction opposite the median difference.

3. **The location-shift assumption** (identical shapes) is required to interpret WMW as testing the Hodges-Lehmann estimator. This assumption is often implausible in practice.

4. **WMWodds = p̂/(1-p̂)** provides an interpretable effect size that can be compared across studies, analogous to an odds ratio.

5. **The WMW test is valid for tied data.** This is established by Lehmann (1975), despite claims to the contrary in some textbooks.

6. **The signed-rank test requires symmetry** of the difference distribution; without it, use the sign test.

7. **Equivalence testing via TOST** is valid on the Hodges-Lehmann scale, but the bounds should be interpreted as pseudomedian differences, not necessarily mean or median differences.

8. **Confidence intervals and p-values** from `wilcox.test()` are concordant (both based on Hodges-Lehmann), with minor boundary discrepancies due to discreteness.

---

## References

- Agresti, A. (1980). Generalized odds ratios for ordinal data. *Biometrics*, 36, 59-67.
- Bamber, D. (1975). The area above the ordinal dominance graph and the area below the receiver operating characteristic graph. *Journal of Mathematical Psychology*, 12, 387-415.
- Brunner, E., & Munzel, U. (2000). The nonparametric Behrens-Fisher problem: Asymptotic theory and a small-sample approximation. *Biometrical Journal*, 42, 17-25.
- Divine, G. W., Norton, H. J., Barón, A. E., & Juarez-Colunga, E. (2018). The Wilcoxon-Mann-Whitney procedure fails as a test of medians. *The American Statistician*, 72(3), 278-286.
- Divine, G., Norton, H. J., Hunt, R., & Dienemann, J. (2013). A review of analysis and sample size calculation considerations for Wilcoxon tests. *Anesthesia & Analgesia*, 117(3), 699-710.
- Fligner, M. A., & Policello, G. E. (1981). Robust rank procedures for the Behrens-Fisher problem. *Journal of the American Statistical Association*, 76, 162-168.
- Hanley, J., & McNeil, B. (1982). The meaning and use of the area under a receiver operating characteristic (ROC) curve. *Radiology*, 143, 29-36.
- Hodges, J. L., & Lehmann, E. L. (1963). Estimates of location based on rank tests. *Annals of Mathematical Statistics*, 34(2), 598-611.
- Hunt, R., Dienemann, J., Norton, H., Hartley, W., Hudgens, A., Stern, T., & Divine, G. (2013). Aromatherapy as treatment for postoperative nausea. *Anesthesia & Analgesia*, 117, 597-604.
- Lehmann, E. L. (1975). *Nonparametrics: Statistical Methods Based on Ranks*. San Francisco: Holden-Day.
- Newson, R. (2002). Parameters behind "nonparametric" statistics: Kendall's tau, Somers' D and median differences. *Stata Journal*, 2, 45-64.
- O'Brien, R. G., & Castelloe, J. M. (2006). Exploiting the link between the Wilcoxon-Mann-Whitney test and a simple odds statistic. *Proceedings of the 31st Annual SAS Users Group International Conference*, Paper 209-31.
- Wellek, S. (2010). *Testing Statistical Hypotheses of Equivalence and Noninferiority* (2nd ed.). CRC Press.
