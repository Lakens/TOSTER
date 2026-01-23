# TOSTER Package: Conceptual and Methodological Reference

This document provides background on equivalence testing methodology that underlies the TOSTER R package. It is derived from the foundational papers by Lakens and colleagues.

---

## Core Concept: What Equivalence Testing Does

Equivalence testing determines whether an observed effect is small enough to be considered practically equivalent to zero (or another specified value). Unlike null-hypothesis significance testing (NHST), which can only reject a null hypothesis of zero difference, equivalence testing allows researchers to reject the presence of meaningful effects.

**Key distinction**: A non-significant NHST result does not demonstrate absence of an effect. Equivalence testing provides a formal statistical framework to conclude that effects are too small to matter.

---

## The Two One-Sided Tests (TOST) Procedure

### Mechanism

TOST tests whether an observed effect falls within pre-specified equivalence bounds (ΔL and ΔU). The procedure performs two one-sided tests:

1. Test whether the effect is significantly greater than the lower bound (ΔL)
2. Test whether the effect is significantly less than the upper bound (ΔU)

If both tests are significant, the effect is declared statistically equivalent to zero (within the specified bounds).

### Interpretation via Confidence Intervals

For alpha = 0.05, equivalence is concluded when the 90% confidence interval falls entirely within the equivalence bounds. This is mathematically equivalent to performing the two one-sided tests.

**Note**: The 90% CI is used (not 95%) because TOST consists of two one-sided tests at alpha = 0.05 each. The CI width is 1 − 2α.

### No Multiple Comparison Correction Needed

Both one-sided tests must be significant to conclude equivalence, so no correction for multiple comparisons is required. When reporting, it suffices to report the test with the larger p-value (i.e., the more conservative result).

---

## Types of Tests

### Equivalence Test (Two-Sided)
- **H0**: Effect is ≤ ΔL OR ≥ ΔU (effect is meaningfully different from zero)
- **H1**: Effect is between ΔL and ΔU (effect is practically equivalent to zero)
- Requires symmetric or asymmetric bounds around zero (or another reference value)

### Inferiority/Superiority Test (One-Sided)
- Tests whether an effect is reliably smaller (or larger) than a single bound
- Appropriate when only one direction matters (e.g., replication studies testing if effect is smaller than original)
- Uses only one equivalence bound (Δ)

### Minimal Effects Test
- The logical inverse of equivalence testing
- **H0**: Effect falls between bounds
- **H1**: Effect is outside bounds
- Less commonly used but conceptually important for understanding the framework

---

## Possible Outcomes When Combining NHST and Equivalence Tests

| NHST Result | Equivalence Result | Interpretation |
|-------------|-------------------|----------------|
| Not significant | Not equivalent | Inconclusive; insufficient data |
| Not significant | Equivalent | Effect is trivially small |
| Significant | Not equivalent | Meaningful effect exists |
| Significant | Equivalent | Statistically significant but practically insignificant effect |

The fourth outcome (significant AND equivalent) occurs with very large samples where even trivially small effects become statistically detectable.

---

## Smallest Effect Size of Interest (SESOI)

The SESOI defines the equivalence bounds and determines what question the test answers. Bounds can be symmetric (±Δ) or asymmetric (ΔL ≠ |ΔU|).

### Approaches to Justifying SESOI

**Objective justifications:**
- Just-noticeable differences (perceptual thresholds)
- Minimal clinically important differences (patient-reported thresholds)
- Quantitative theoretical predictions

**Subjective justifications:**
- *Based on previous studies*: Small telescopes approach (effect size original study had 33% power to detect), critical effect size (smallest detectable effect in original study), or meta-analytic estimates
- *Based on resources*: Effect size the planned study has adequate power to detect (answers a resource question rather than theoretical question)
- *Benchmarks*: Cohen's conventions (d = 0.1, 0.2, 0.5) — weakest justification, generally discouraged

### Raw vs. Standardized Bounds

- **Raw bounds** (e.g., 5 points on a scale): Independent of standard deviation; appropriate when raw units are meaningful
- **Standardized bounds** (e.g., d = 0.3): Depend on observed SD; useful for cross-study comparisons or when raw units lack inherent meaning

These ask slightly different questions and may yield different results if SDs vary substantially.

---

## Power and Sample Size Considerations

Equivalence tests require adequate sample size to obtain confidence intervals narrow enough to fall within the equivalence bounds. As the SESOI becomes smaller (narrower bounds), larger samples are needed.

Power for equivalence tests should be calculated a priori, similar to NHST. The study should have sufficient power both to:
1. Detect effects that exceed the SESOI (via NHST)
2. Demonstrate equivalence if the true effect is near zero (via TOST)

---

## Statistical Frameworks

### Frequentist (TOST)
- Controls long-run error rates
- Provides dichotomous decisions (equivalent vs. not equivalent)
- Based on Neyman-Pearson hypothesis testing framework
- Cannot conclude the true effect is exactly zero; only that effects larger than the SESOI can be rejected

### Bayesian Alternatives (for reference only)

**Region of Practical Equivalence (ROPE) procedure:**
- Uses Bayesian estimation to generate posterior distribution
- Compares 95% Highest Density Interval to equivalence bounds
- Similar decision rules to TOST but based on credible intervals rather than confidence intervals

**Bayes factors:**
- Compare relative evidence for null model vs. alternative model
- Require specification of prior distributions for the alternative hypothesis
- Provide continuous measure of evidence rather than dichotomous decision

---

## Common Applications

1. **Replication studies**: Testing whether a replication effect is too small to support the original finding
2. **Non-inferiority trials**: Testing whether a new treatment is not meaningfully worse than standard treatment
3. **Bioequivalence**: Testing whether generic drugs produce equivalent effects to brand-name drugs
4. **Ruling out confounds**: Demonstrating that groups do not differ meaningfully on potential confounding variables
5. **Interpreting null results**: Distinguishing between "no evidence of effect" and "evidence of no meaningful effect"

---

## Reporting Guidelines

When reporting equivalence tests:

1. State the SESOI and justify its selection (before seeing data)
2. Report the equivalence bounds in the same units as the effect (raw or standardized)
3. Report the t-statistic (or other test statistic) and p-value for the more conservative one-sided test
4. Report the relevant confidence interval (90% for alpha = 0.05)
5. State whether equivalence was concluded
6. If conducting both NHST and equivalence tests, report both results

**Example phrasing**: "The equivalence test was significant, t(97.78) = −1.90, p = .030, indicating that effects as large as d = 0.49 can be rejected."

---

## Key References

- Lakens, D. (2017). Equivalence Tests: A Practical Primer for t Tests, Correlations, and Meta-Analyses. *Social Psychological and Personality Science*, 8(4), 355–362.

- Lakens, D., Scheel, A. M., & Isager, P. M. (2018). Equivalence Testing for Psychological Research: A Tutorial. *Advances in Methods and Practices in Psychological Science*, 1(2), 259–269.

- Harms, C., & Lakens, D. (2018). Making 'Null Effects' Informative: Statistical Techniques and Inferential Frameworks. *Journal of Clinical and Translational Research*.

- Schuirmann, D. J. (1987). A comparison of the Two One-Sided Tests Procedure and the Power Approach for assessing the equivalence of average bioavailability. *Journal of Pharmacokinetics and Biopharmaceutics*, 15(6), 657–680.

- Rogers, J. L., Howard, K. I., & Vessey, J. T. (1993). Using significance tests to evaluate equivalence between two experimental groups. *Psychological Bulletin*, 113(3), 553–565.
