# Paired Brunner-Munzel: `t` vs `logit` vs `perm` by sample size

*Simulation study, 2026-09-28. Companion to `junk/bm_t_vs_perm_sample_size.md` (two-sample).
Scripts and results live in `junk/` (not part of the package).*

## Question

For `brunner_munzel(..., paired = TRUE)`, how do `test_method = "t"`, `"logit"`, and
`"perm"` behave in small samples? From what sample size does `"t"` give essentially the
same inference as `"perm"`? Does the two-sample answer (about 15 per group) carry over?

## Short answer

**No, the two-sample answer does not carry over.** For paired data the within-pair
correlation matters more than the sample size.

- **Uncorrelated pairs (ρ = 0):** similar to the two-sample case. `t` is practically
  equivalent to `perm` from about **13–20 pairs**.
- **Positively correlated pairs (ρ = 0.5, 0.8), the usual paired design:**
  - **`t` is conservative.** With identical marginals and ρ = 0.8 it rejects 1.8–3.8%
    for n = 7–50, while `perm` stays near 5%.
  - At ρ = 0.8 the gap is still about 1 percentage point at **n = 50**, with about 1%
    decision disagreement. At ρ = 0.5 it shrinks to near zero by n = 50 with identical
    marginals.
  - Conservative is not invalid, but it means `perm` is more powerful. **For correlated
    pairs, `perm` is preferable well beyond n = 15, at least up to 50 pairs.**
- **Five or fewer pairs:** `perm` **cannot reject at two-sided α = .05.** Full
  enumeration of 2^5 = 32 swaps with the symmetric |T| rule gives a smallest possible
  p-value of 2/32 = .0625. This is a property of the design, not a bug. Six pairs
  (64 swaps, smallest p = .031) is the minimum at which a rejection is possible.
- **`logit`:** very conservative in small samples (Type I error around 1–3% up to about
  15 pairs), and still below nominal at 50 pairs (.036–.045).
- **Two limits of `perm`:**
  - With **ties plus strong correlation** (ordinal, ρ = 0.8) it is conservative in small
    samples (.003–.03 up to n = 15). Many pairs have x = y, and swapping those changes
    nothing.
  - With **unequal marginals plus strong correlation** it is slightly liberal
    (.052–.062). Within-pair swapping assumes (X, Y) and (Y, X) have the same
    distribution, which fails when the marginals differ. This is the paired version of
    the weak-null issue.

## Design

- **Methods:**
  - `brunner_munzel(x, y, paired = TRUE, test_method = m)` for m in `"t"`, `"logit"`
    and `"perm"`, all applied to **the same datasets**.
  - For `perm`, R = 999. The paired branch fully enumerates all 2^n within-pair swaps
    when n ≤ 13, regardless of `R`, and uses `p_method = "exact"` there. It samples
    R swap patterns when n > 13 and uses `p_method = "plusone"`.
- **Null:** the relative effect p = P(X > Y) + 0.5 P(X = Y) is defined on the marginal
  distributions, so it does not depend on the within-pair correlation. Every scenario
  has p = 0.5, so every rejection is a Type I error. A Monte Carlo check with 2e5 draws
  gave 0.4985–0.4998.
- **Pairs:** generated with a Gaussian copula. Correlated normals with correlation ρ are
  mapped to the target marginals through their quantile functions. ρ is the latent
  (normal-scale) correlation.
- **Two-sided α:** 0.05. **Replications:** 4,000 per cell. The Monte Carlo SE of one
  rate near .05 is about .0035. The paired SE of the `t` − `perm` difference is
  √(disagreement / nsim).
- **Grid:**
  - n (pairs) ∈ {5, 7, 10, 13, 15, 20, 30, 50}
  - ρ ∈ {0, 0.5, 0.8}
  - four marginal pairs (x vs y):

| label | x marginal | y marginal | what it probes |
|---|---|---|---|
| identical_normal | N(0, 1) | N(0, 1) | exchangeable within pairs |
| identical_ordinal | 5-point, probs (.1,.2,.4,.2,.1) | same | exchangeable, heavy ties |
| normal_sd3_vs_sd1 | N(0, 3) | N(0, 1) | unequal spread |
| skewed_scale3_vs_1 | 3·Exp(1) − c | Exp(1) | skewed, unequal scale |

With identical marginals and a Gaussian copula, (X, Y) and (Y, X) have the same
distribution, so the paired permutation test is exact. With unequal marginals it is
not exact.

All results use the code after the floating point tie fix (tolerant permutation counts),
so no merging of runs was needed.

### Criterion for "practically equivalent"

A cell counts as equivalent if |t − perm| ≤ tol and disagreement ≤ 1%. The reported n
is the smallest n from which **every** larger n in the grid is also equivalent. NA means
the criterion is not met even at n = 50.

## Results

### Sample size at which `t` becomes practically equivalent to `perm`

| marginals | ρ | tol = 0.005 | tol = 0.010 |
|---|---|---|---|
| identical_normal | 0 | 15 | 13 |
| identical_normal | 0.5 | 50 | 50 |
| identical_normal | 0.8 | NA (> 50) | NA (> 50) |
| identical_ordinal | 0 | 30 | 20 |
| identical_ordinal | 0.5 | 20 | 20 |
| identical_ordinal | 0.8 | NA (> 50) | 30 |
| normal_sd3_vs_sd1 | 0 | 20 | 10 |
| normal_sd3_vs_sd1 | 0.5 | 50 | 20 |
| normal_sd3_vs_sd1 | 0.8 | NA (> 50) | 30 |
| skewed_scale3_vs_1 | 0 | 13 | 13 |
| skewed_scale3_vs_1 | 0.5 | 20 | 20 |
| skewed_scale3_vs_1 | 0.8 | NA (> 50) | NA (> 50) |

For the ρ = 0.8 cells that meet the 1-point tolerance at n = 30, `t` is still
conservative by .005–.01, just inside the tolerance.

### Main patterns

1. **n = 5.** `perm` rejects 0% everywhere, because its smallest attainable two-sided
   p-value is .0625. The large `t` − `perm` differences at n = 5 in the tables reflect
   this, not a failure of `t`.
2. **ρ = 0.**
   - `t` is mildly liberal at 5–15 pairs (.051–.058) and within about half a point of
     `perm` from 13–20 pairs. This mirrors the two-sample result.
   - With ordinal data at n = 15, both methods are somewhat high (t .067, perm .057).
3. **ρ > 0.**
   - `t` is consistently **conservative**, more so as ρ increases. Its Type I error
     climbs slowly toward .05 as n grows but is still around .038–.048 at n = 50 when
     ρ = 0.8.
   - `perm` stays close to .05 with identical marginals.
   - The paired z-statistics for `t` − `perm` are −5 to −11 through n = 50 at ρ = 0.8,
     so this is a clear systematic difference, not noise.
   - The cause is likely the placement-based variance estimator and the t reference
     distribution under strong within-pair dependence. It was not investigated further.
4. **Unequal marginals with correlation.**
   - `perm` is slightly liberal at ρ = 0.8: normal_sd3_vs_sd1 .052–.062, skewed
     .051–.058.
   - Swapping within a pair is only guaranteed valid when (X, Y) and (Y, X) have the
     same distribution. The studentized statistic makes this asymptotically valid, not
     exact.
   - `t` is conservative in the same cells, so the two methods err in opposite
     directions.
5. **Ties with strong correlation.** In ordinal ρ = 0.8, `perm` is .003 at n = 7 and
   .033 at n = 15, reaching about .04–.05 by n = 20–30. Pairs with x = y add nothing to
   the permutation distribution, which becomes very discrete. `t` is also conservative
   here (.017–.043).
6. **`logit`.** Conservative in every cell: about .006–.03 up to n = 15 and .036–.045 at
   n = 50. It trades considerable power for range-preserving intervals.

## Practical guidance (paired)

- **Five or fewer pairs:** `perm` cannot reject at two-sided α = .05. Use `t` with caution
  (its level is uncertain at this n), or treat the design as underpowered.
- **6–15 pairs:** prefer `perm`.
- **More than 15 pairs:**
  - **Uncorrelated pairs:** `t` is fine.
  - **Positively correlated pairs,** the typical case: `perm` still holds its level more
    accurately and is more powerful than the conservative `t`, at least up to 50 pairs.
    `t` remains valid (conservative) if speed matters.
- **Unequal marginals with strong correlation:** neither method is exact. `perm` may be
  slightly liberal (about .06), and `t` conservative.

## Caveats

- **Type I error only.** Power was not simulated. `t`'s conservatism suggests a power
  loss relative to `perm` for correlated pairs, but its size was not measured.
- **Two-sided α = .05 only.**
- **Latent correlation.** ρ is the latent (normal-scale) correlation of the Gaussian
  copula. The correlation on the observed scale is somewhat lower for skewed and ordinal
  marginals.
- **Monte Carlo floors.** For n > 13, `perm` uses sampled permutations with R = 999. That
  sets floors of about 0.5% decision disagreement and about .010 mean |p_t − p_perm|.

## Full results

Columns: Type I error for each method, paired `t` − `perm` difference with its
z-statistic, decision disagreement, and mean absolute p-value difference.

#### identical_normal 

| rho | n | t | logit | perm | t - perm (z) | disagree | mean \|p_t - p_perm\| |
|---|---|---|---|---|---|---|---|
| 0 | 5 | 0.058 | 0.011 | 0.000 | 0.058 (15.2) | 0.058 | 0.056 |
| 0 | 7 | 0.052 | 0.013 | 0.046 | 0.007 (3.1) | 0.019 | 0.029 |
| 0 | 10 | 0.053 | 0.022 | 0.046 | 0.006 (3.9) | 0.011 | 0.015 |
| 0 | 13 | 0.055 | 0.026 | 0.050 | 0.005 (3.5) | 0.010 | 0.011 |
| 0 | 15 | 0.051 | 0.030 | 0.050 | 0.001 (1.1) | 0.007 | 0.013 |
| 0 | 20 | 0.053 | 0.034 | 0.050 | 0.003 (2.0) | 0.007 | 0.011 |
| 0 | 30 | 0.054 | 0.041 | 0.052 | 0.002 (1.3) | 0.005 | 0.011 |
| 0 | 50 | 0.052 | 0.045 | 0.051 | 0.000 (0.5) | 0.004 | 0.010 |
| 0.5 | 5 | 0.038 | 0.015 | 0.000 | 0.038 (12.4) | 0.038 | 0.059 |
| 0.5 | 7 | 0.034 | 0.013 | 0.045 | -0.011 (-4.7) | 0.021 | 0.040 |
| 0.5 | 10 | 0.032 | 0.016 | 0.050 | -0.018 (-8.2) | 0.018 | 0.027 |
| 0.5 | 13 | 0.042 | 0.024 | 0.051 | -0.009 (-5.4) | 0.011 | 0.021 |
| 0.5 | 15 | 0.044 | 0.026 | 0.053 | -0.010 (-5.3) | 0.013 | 0.020 |
| 0.5 | 20 | 0.042 | 0.030 | 0.050 | -0.008 (-5.1) | 0.009 | 0.016 |
| 0.5 | 30 | 0.048 | 0.040 | 0.058 | -0.009 (-5.8) | 0.010 | 0.012 |
| 0.5 | 50 | 0.044 | 0.038 | 0.046 | -0.001 (-1.1) | 0.005 | 0.011 |
| 0.8 | 5 | 0.033 | 0.015 | 0.000 | 0.033 (11.4) | 0.033 | 0.065 |
| 0.8 | 7 | 0.025 | 0.016 | 0.050 | -0.025 (-9.5) | 0.028 | 0.058 |
| 0.8 | 10 | 0.018 | 0.013 | 0.042 | -0.024 (-9.6) | 0.024 | 0.047 |
| 0.8 | 13 | 0.028 | 0.021 | 0.054 | -0.027 (-10.4) | 0.027 | 0.040 |
| 0.8 | 15 | 0.025 | 0.019 | 0.054 | -0.028 (-10.6) | 0.028 | 0.036 |
| 0.8 | 20 | 0.029 | 0.024 | 0.048 | -0.018 (-8.5) | 0.018 | 0.027 |
| 0.8 | 30 | 0.038 | 0.034 | 0.056 | -0.018 (-8.5) | 0.018 | 0.020 |
| 0.8 | 50 | 0.038 | 0.036 | 0.049 | -0.010 (-6.3) | 0.011 | 0.014 |

#### identical_ordinal 

| rho | n | t | logit | perm | t - perm (z) | disagree | mean \|p_t - p_perm\| |
|---|---|---|---|---|---|---|---|
| 0 | 5 | 0.070 | 0.015 | 0.000 | 0.070 (16.7) | 0.070 | 0.120 |
| 0 | 7 | 0.064 | 0.022 | 0.014 | 0.050 (14.1) | 0.050 | 0.068 |
| 0 | 10 | 0.055 | 0.025 | 0.034 | 0.020 (8.6) | 0.022 | 0.034 |
| 0 | 13 | 0.049 | 0.028 | 0.038 | 0.011 (6.4) | 0.012 | 0.019 |
| 0 | 15 | 0.067 | 0.041 | 0.056 | 0.011 (5.9) | 0.013 | 0.018 |
| 0 | 20 | 0.058 | 0.041 | 0.051 | 0.006 (4.4) | 0.008 | 0.014 |
| 0 | 30 | 0.054 | 0.042 | 0.050 | 0.004 (3.3) | 0.007 | 0.012 |
| 0 | 50 | 0.051 | 0.043 | 0.047 | 0.004 (2.7) | 0.006 | 0.011 |
| 0.5 | 5 | 0.047 | 0.013 | 0.000 | 0.047 (13.7) | 0.047 | 0.154 |
| 0.5 | 7 | 0.037 | 0.015 | 0.006 | 0.032 (11.1) | 0.032 | 0.096 |
| 0.5 | 10 | 0.034 | 0.020 | 0.024 | 0.011 (5.3) | 0.018 | 0.049 |
| 0.5 | 13 | 0.049 | 0.035 | 0.044 | 0.005 (2.9) | 0.011 | 0.029 |
| 0.5 | 15 | 0.042 | 0.030 | 0.042 | 0.000 (0.3) | 0.011 | 0.026 |
| 0.5 | 20 | 0.047 | 0.040 | 0.048 | -0.001 (-1.0) | 0.007 | 0.017 |
| 0.5 | 30 | 0.051 | 0.044 | 0.052 | -0.001 (-0.7) | 0.005 | 0.013 |
| 0.5 | 50 | 0.046 | 0.042 | 0.048 | -0.002 (-2.0) | 0.005 | 0.011 |
| 0.8 | 5 | 0.017 | 0.006 | 0.000 | 0.017 (8.2) | 0.017 | 0.210 |
| 0.8 | 7 | 0.023 | 0.010 | 0.002 | 0.021 (9.1) | 0.021 | 0.149 |
| 0.8 | 10 | 0.028 | 0.018 | 0.011 | 0.018 (8.0) | 0.019 | 0.093 |
| 0.8 | 13 | 0.029 | 0.022 | 0.022 | 0.006 (3.6) | 0.013 | 0.061 |
| 0.8 | 15 | 0.035 | 0.031 | 0.032 | 0.002 (1.5) | 0.011 | 0.048 |
| 0.8 | 20 | 0.040 | 0.033 | 0.042 | -0.001 (-0.7) | 0.012 | 0.030 |
| 0.8 | 30 | 0.041 | 0.039 | 0.049 | -0.007 (-4.7) | 0.010 | 0.018 |
| 0.8 | 50 | 0.042 | 0.040 | 0.048 | -0.005 (-4.0) | 0.007 | 0.013 |

#### normal_sd3_vs_sd1 

| rho | n | t | logit | perm | t - perm (z) | disagree | mean \|p_t - p_perm\| |
|---|---|---|---|---|---|---|---|
| 0 | 5 | 0.054 | 0.008 | 0.000 | 0.054 (14.7) | 0.054 | 0.053 |
| 0 | 7 | 0.058 | 0.011 | 0.041 | 0.017 (6.6) | 0.026 | 0.027 |
| 0 | 10 | 0.052 | 0.014 | 0.045 | 0.007 (4.7) | 0.009 | 0.015 |
| 0 | 13 | 0.056 | 0.020 | 0.050 | 0.006 (4.6) | 0.007 | 0.012 |
| 0 | 15 | 0.056 | 0.026 | 0.051 | 0.005 (3.5) | 0.008 | 0.013 |
| 0 | 20 | 0.050 | 0.031 | 0.049 | 0.000 (0.6) | 0.003 | 0.012 |
| 0 | 30 | 0.051 | 0.036 | 0.050 | 0.001 (0.8) | 0.004 | 0.011 |
| 0 | 50 | 0.046 | 0.037 | 0.046 | 0.000 (0.0) | 0.005 | 0.010 |
| 0.5 | 5 | 0.051 | 0.015 | 0.000 | 0.051 (14.3) | 0.051 | 0.056 |
| 0.5 | 7 | 0.040 | 0.013 | 0.043 | -0.003 (-1.6) | 0.020 | 0.036 |
| 0.5 | 10 | 0.042 | 0.016 | 0.050 | -0.008 (-4.6) | 0.013 | 0.026 |
| 0.5 | 13 | 0.041 | 0.019 | 0.049 | -0.008 (-5.2) | 0.010 | 0.021 |
| 0.5 | 15 | 0.048 | 0.022 | 0.054 | -0.006 (-3.8) | 0.011 | 0.019 |
| 0.5 | 20 | 0.045 | 0.028 | 0.048 | -0.003 (-2.1) | 0.007 | 0.016 |
| 0.5 | 30 | 0.042 | 0.031 | 0.050 | -0.007 (-4.9) | 0.009 | 0.012 |
| 0.5 | 50 | 0.052 | 0.044 | 0.056 | -0.004 (-3.5) | 0.006 | 0.011 |
| 0.8 | 5 | 0.045 | 0.017 | 0.000 | 0.045 (13.5) | 0.045 | 0.066 |
| 0.8 | 7 | 0.022 | 0.009 | 0.056 | -0.033 (-11.1) | 0.036 | 0.052 |
| 0.8 | 10 | 0.037 | 0.014 | 0.058 | -0.021 (-9.1) | 0.021 | 0.041 |
| 0.8 | 13 | 0.035 | 0.018 | 0.058 | -0.023 (-9.6) | 0.023 | 0.034 |
| 0.8 | 15 | 0.043 | 0.023 | 0.062 | -0.019 (-8.7) | 0.019 | 0.031 |
| 0.8 | 20 | 0.038 | 0.024 | 0.052 | -0.014 (-7.6) | 0.014 | 0.024 |
| 0.8 | 30 | 0.040 | 0.033 | 0.048 | -0.008 (-5.8) | 0.008 | 0.017 |
| 0.8 | 50 | 0.043 | 0.040 | 0.052 | -0.008 (-5.8) | 0.008 | 0.013 |

#### skewed_scale3_vs_1 

| rho | n | t | logit | perm | t - perm (z) | disagree | mean \|p_t - p_perm\| |
|---|---|---|---|---|---|---|---|
| 0 | 5 | 0.055 | 0.006 | 0.000 | 0.055 (14.8) | 0.055 | 0.052 |
| 0 | 7 | 0.051 | 0.011 | 0.036 | 0.014 (5.7) | 0.026 | 0.027 |
| 0 | 10 | 0.052 | 0.012 | 0.044 | 0.009 (4.9) | 0.013 | 0.015 |
| 0 | 13 | 0.047 | 0.020 | 0.045 | 0.002 (1.2) | 0.006 | 0.012 |
| 0 | 15 | 0.054 | 0.029 | 0.051 | 0.002 (1.8) | 0.006 | 0.014 |
| 0 | 20 | 0.051 | 0.031 | 0.049 | 0.002 (1.5) | 0.005 | 0.012 |
| 0 | 30 | 0.054 | 0.042 | 0.053 | 0.001 (1.0) | 0.007 | 0.011 |
| 0 | 50 | 0.050 | 0.042 | 0.049 | 0.001 (0.6) | 0.006 | 0.010 |
| 0.5 | 5 | 0.045 | 0.010 | 0.000 | 0.045 (13.4) | 0.045 | 0.056 |
| 0.5 | 7 | 0.036 | 0.012 | 0.045 | -0.009 (-3.8) | 0.022 | 0.037 |
| 0.5 | 10 | 0.041 | 0.014 | 0.048 | -0.008 (-4.8) | 0.010 | 0.025 |
| 0.5 | 13 | 0.051 | 0.024 | 0.056 | -0.006 (-4.0) | 0.007 | 0.020 |
| 0.5 | 15 | 0.048 | 0.026 | 0.057 | -0.009 (-5.5) | 0.010 | 0.019 |
| 0.5 | 20 | 0.046 | 0.028 | 0.050 | -0.004 (-3.4) | 0.006 | 0.015 |
| 0.5 | 30 | 0.047 | 0.035 | 0.049 | -0.002 (-2.0) | 0.006 | 0.012 |
| 0.5 | 50 | 0.048 | 0.042 | 0.050 | -0.002 (-2.2) | 0.004 | 0.011 |
| 0.8 | 5 | 0.048 | 0.014 | 0.000 | 0.048 (13.9) | 0.048 | 0.066 |
| 0.8 | 7 | 0.029 | 0.009 | 0.051 | -0.022 (-8.2) | 0.028 | 0.052 |
| 0.8 | 10 | 0.034 | 0.012 | 0.053 | -0.019 (-8.6) | 0.020 | 0.041 |
| 0.8 | 13 | 0.037 | 0.019 | 0.052 | -0.016 (-7.9) | 0.016 | 0.034 |
| 0.8 | 15 | 0.036 | 0.022 | 0.057 | -0.020 (-9.1) | 0.020 | 0.030 |
| 0.8 | 20 | 0.043 | 0.030 | 0.058 | -0.015 (-7.6) | 0.014 | 0.024 |
| 0.8 | 30 | 0.043 | 0.035 | 0.056 | -0.013 (-7.1) | 0.013 | 0.017 |
| 0.8 | 50 | 0.048 | 0.042 | 0.058 | -0.010 (-6.3) | 0.011 | 0.013 |

## Reproduction

```r
Sys.setenv(SIM_NSIM = 4000)
source("junk/sim_brunner_munzel_paired_t_vs_perm_n.R")   # about 1 hour on 13 cores
```

Files:

- `junk/sim_brunner_munzel_paired_t_vs_perm_n.R`: simulation script
- `junk/sim_brunner_munzel_paired_t_vs_perm_n_results.rds`: results, including the
  equivalence summary (`equiv_n`) and the Monte Carlo check of p = 0.5 (`p_check`)
