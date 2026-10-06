# Brunner-Munzel: when does the `t` method match the permutation test?

*Simulation study, 2026-09-26. Scripts and results live in `junk/` (not part of the package).*

## Question

`brunner_munzel()` offers a t-approximation (`test_method = "t"`, Welch-Satterthwaite df)
and a studentized permutation test (`test_method = "perm"`, Neubert & Brunner, 2007).
In small samples the permutation test is recommended. **From what sample size does the
`t` method give essentially the same inference as `perm`?**

## Short answer

- With **15 or more observations in the smaller group**, `t` and `perm` are practically
  equivalent in every scenario studied:
  - Their Type I error rates differ by at most 1 percentage point.
  - They reach different decisions on at most 1% of datasets.
- For **continuous data**, a stricter tolerance of 0.5 percentage points is met from
  about **10–15 in the smaller group**.
- For **heavily tied (ordinal) data**, the stricter tolerance needs about **15–20**.
  One balanced cell sits right on the 0.5-point boundary up to n = 40.
- Below about 10 per group the `t` method is clearly **liberal**:
  - At 5 v 5 it rejects 7–9% of the time.
  - The permutation test stays near or below 5% when the groups are exchangeable.
- The package's existing advice fits these results. It messages users to prefer `perm`
  when either group has fewer than 15 observations.
- Above about 15–20 per group, `t` is an adequate and much faster substitute. Even at
  50 v 100, `t` stays slightly more liberal than `perm` in some designs, but only by
  0.1–0.4 percentage points (see Caveats).

## Design

- **Methods:** `brunner_munzel(x, y, test_method = "t")` and `test_method = "perm"`,
  with R = 999 and the default `p_method`. The permutation test fully enumerates when
  there are at most 999 permutations (only 5 v 5) and uses `plusone` otherwise.
  Both methods are applied to **the same datasets**, so their comparison is paired.
- **Null:** in every scenario the true relative effect
  p = P(X > Y) + 0.5 P(X = Y) is exactly 0.5, so every rejection of H0: p = 0.5
  is a Type I error. For the skewed pair, the shift that makes p = 0.5 was solved by
  numerical integration.
- **Two-sided α:** 0.05.
- **Replications:** 4,000 per cell. The Monte Carlo SE of a single Type I rate near
  .05 is about .0035. The paired SE of the difference is much smaller:
  √(disagreement rate / nsim), roughly .001–.002 for n ≥ 10.
- **Sample sizes:** smaller-group n ∈ {5, 7, 10, 15, 20, 25, 30, 40, 50}. Each is run
  balanced (n v n) and unbalanced in both directions: (n v 2n) and (2n v n).
  x is always listed first.
- **Distribution pairs (x vs y):**

| label | x | y | what it probes |
|---|---|---|---|
| identical_normal | N(0, 1) | N(0, 1) | exchangeable control |
| identical_ordinal | 5-point, probs (.1,.2,.4,.2,.1) | same | exchangeable, heavy ties |
| normal_sd3_vs_sd1 | N(0, 3) | N(0, 1) | unequal variance, same shape |
| skewed_scale3_vs_1 | 3·Exp(1) − c | Exp(1) | skewed, unequal scale |

In the 1:2 designs the smaller group is the more variable one (x). In the 2:1 designs the
larger group is the more variable one.

### Metrics per cell

- `t`, `perm`: two-sided Type I error.
- `t − perm (z)`: paired difference in Type I error, with its z-statistic.
- `disagree`: proportion of datasets on which the two methods reach different decisions.
- `mean |p_t − p_perm|`: average absolute difference in p-values.

### Criterion for "practically equivalent"

A cell counts as equivalent if |t − perm| ≤ tol **and** disagreement ≤ 1%. The reported
n is the smallest smaller-group n from which **every** larger n is also equivalent.

## Important: floating point fix during this study

The first run showed the *exact* permutation test above .05 in exchangeable 5 v 5 cells
(normal .058, ordinal .080). Full enumeration is exact under exchangeability, so that
pointed to a bug.

**Cause.** The observed statistic is computed from ranks. The permuted statistics are
computed by `perm_loop()` from pairwise comparisons. Mathematically tied values
therefore differed by about 1e-16, and strict comparisons such as `Tperm >= t_obs`
dropped them from the count. Among those dropped values was the observed arrangement's
own permutation, in about a third of datasets.

**Fix.** Permutation counts in `perm_t_test()`, `brunner_munzel()` and
`hodges_lehmann()` now treat values within a small relative tolerance as ties
(`perm_count()` / `perm_tol()` in `R/perm_t_test.R`). `perm_crit()` uses the same
tolerance, so CIs still agree with p-values. After the fix, exact 5 v 5 ordinal
rejects 1.4% of the time and normal 4.0%.

**All results below use the corrected code:**
- **Rerun with the corrected code:**
  - every cell with n ≤ 10
  - every `identical_ordinal` cell, because tied data are affected at any n
- **Kept from the original run:** continuous-data cells with n > 10. With continuous
  data and sampled permutations, exact ties essentially never occur, so the fix does
  not change them.
- **Merging:** `junk/bm_t_vs_perm_merge.R` combines the two runs.

## Results

### Sample size at which `t` becomes practically equivalent to `perm`

Smaller-group n from which all larger n meet the criterion:

| distribution | design | tol = 0.005 | tol = 0.010 |
|---|---|---|---|
| identical_normal | 1:1 | 10 | 10 |
| identical_normal | 1:2 | 10 | 7 |
| identical_normal | 2:1 | 25 | 10 |
| identical_ordinal | 1:1 | 50* | 15 |
| identical_ordinal | 1:2 | 20 | 15 |
| identical_ordinal | 2:1 | 15 | 10 |
| normal_sd3_vs_sd1 | 1:1 | 10 | 10 |
| normal_sd3_vs_sd1 | 1:2 | 10 | 10 |
| normal_sd3_vs_sd1 | 2:1 | 15 | 10 |
| skewed_scale3_vs_1 | 1:1 | 10 | 10 |
| skewed_scale3_vs_1 | 1:2 | 15 | 15 |
| skewed_scale3_vs_1 | 2:1 | 15 | 15 |

\* This cell is driven by single cells at 30 v 30 and 40 v 40 where the difference is
.004–.005, essentially on the boundary. From n = 20 onward the difference is ≤ .005 apart
from those two cells.

The identical_normal 2:1 value of 25 at the strict tolerance comes from one cell
(40 v 20, difference .005).

### Main patterns

1. **Very small samples (n ≤ 7).**
   - With exchangeable data, `t` is liberal: 5 v 5 rejects .083 (normal) and .067
     (ordinal), and 5 v 10 or 10 v 5 about .07.
   - `perm` is at or below .05 in these exchangeable cells. It is conservative with
     heavy ties (5 v 5 ordinal: .020), as expected for a discrete exact test.
   - The two methods disagree on 2–5% of datasets.
2. **Unequal variance or skew, with the smaller group more variable (1:2).**
   - Both methods are liberal at 5 v 10 (t .082/.079; perm .079/.074).
   - This is a weak-null (Behrens–Fisher-type) problem that affects both methods, not a
     difference between them.
   - From 7 v 14 on, `t` is actually slightly *less* liberal than `perm` (by about .005).
3. **Unequal variance or skew, with the larger group more variable (2:1).**
   - `perm` is conservative in small samples (.029–.043 up to 20 v 10), while `t` is
     near nominal.
   - The difference shrinks below .005 by 30 v 15.
4. **Large samples.**
   - Disagreement levels off near 0.5% and mean |p_t − p_perm| near .010. Both floors
     reflect Monte Carlo error in `perm`'s own p-value at R = 999 (SE about .01 near
     p = .05), not a real difference between the methods.
   - A small systematic offset remains detectable thanks to the paired design (|z| ≈ 2–4):
     - For n ≥ 25, `t` averages 0.1–0.4 points *above* `perm` for exchangeable and 2:1
       designs.
     - It averages 0.1–0.2 points *below* for 1:2 unequal-spread designs.
   - This is negligible in practice.

### Full results

Columns: Type I error for each method, paired difference with its z-statistic, decision
disagreement, and mean absolute p-value difference.

#### identical_normal

| ratio | n (x v y) | t | perm | t − perm (z) | disagree | mean \|p_t − p_perm\| |
|---|---|---|---|---|---|---|
| 1:1 | 5 v 5 | 0.083 | 0.040 | 0.043 (13.2) | 0.043 | 0.022 |
| 1:1 | 7 v 7 | 0.056 | 0.046 | 0.010 (5.7) | 0.012 | 0.019 |
| 1:1 | 10 v 10 | 0.052 | 0.048 | 0.004 (3.0) | 0.007 | 0.014 |
| 1:1 | 15 v 15 | 0.054 | 0.052 | 0.002 (1.8) | 0.006 | 0.013 |
| 1:1 | 20 v 20 | 0.050 | 0.049 | 0.002 (1.3) | 0.007 | 0.011 |
| 1:1 | 25 v 25 | 0.052 | 0.050 | 0.002 (1.5) | 0.006 | 0.011 |
| 1:1 | 30 v 30 | 0.052 | 0.049 | 0.004 (3.0) | 0.006 | 0.010 |
| 1:1 | 40 v 40 | 0.055 | 0.053 | 0.002 (2.0) | 0.006 | 0.011 |
| 1:1 | 50 v 50 | 0.054 | 0.052 | 0.002 (2.1) | 0.005 | 0.010 |
| 1:2 | 5 v 10 | 0.072 | 0.051 | 0.022 (9.3) | 0.022 | 0.021 |
| 1:2 | 7 v 14 | 0.057 | 0.050 | 0.007 (5.3) | 0.008 | 0.016 |
| 1:2 | 10 v 20 | 0.052 | 0.049 | 0.003 (2.7) | 0.004 | 0.013 |
| 1:2 | 15 v 30 | 0.043 | 0.041 | 0.003 (2.2) | 0.006 | 0.011 |
| 1:2 | 20 v 40 | 0.048 | 0.045 | 0.002 (2.0) | 0.006 | 0.011 |
| 1:2 | 25 v 50 | 0.041 | 0.039 | 0.002 (2.2) | 0.005 | 0.010 |
| 1:2 | 30 v 60 | 0.053 | 0.051 | 0.002 (1.8) | 0.005 | 0.010 |
| 1:2 | 40 v 80 | 0.050 | 0.047 | 0.003 (2.6) | 0.006 | 0.010 |
| 1:2 | 50 v 100 | 0.044 | 0.047 | -0.004 (-2.9) | 0.006 | 0.010 |
| 2:1 | 10 v 5 | 0.071 | 0.051 | 0.020 (8.9) | 0.021 | 0.021 |
| 2:1 | 14 v 7 | 0.062 | 0.051 | 0.010 (6.2) | 0.011 | 0.016 |
| 2:1 | 20 v 10 | 0.058 | 0.051 | 0.007 (4.1) | 0.010 | 0.013 |
| 2:1 | 30 v 15 | 0.058 | 0.056 | 0.002 (2.0) | 0.006 | 0.011 |
| 2:1 | 40 v 20 | 0.055 | 0.050 | 0.005 (4.2) | 0.006 | 0.011 |
| 2:1 | 50 v 25 | 0.053 | 0.050 | 0.003 (2.4) | 0.006 | 0.010 |
| 2:1 | 60 v 30 | 0.045 | 0.046 | -0.000 (-0.5) | 0.004 | 0.010 |
| 2:1 | 80 v 40 | 0.050 | 0.049 | 0.001 (0.6) | 0.006 | 0.010 |
| 2:1 | 100 v 50 | 0.052 | 0.049 | 0.002 (2.4) | 0.004 | 0.010 |

#### identical_ordinal

| ratio | n (x v y) | t | perm | t − perm (z) | disagree | mean \|p_t − p_perm\| |
|---|---|---|---|---|---|---|
| 1:1 | 5 v 5 | 0.067 | 0.020 | 0.046 (13.4) | 0.048 | 0.095 |
| 1:1 | 7 v 7 | 0.055 | 0.028 | 0.027 (10.0) | 0.029 | 0.054 |
| 1:1 | 10 v 10 | 0.059 | 0.041 | 0.018 (8.3) | 0.019 | 0.031 |
| 1:1 | 15 v 15 | 0.058 | 0.049 | 0.009 (5.8) | 0.009 | 0.017 |
| 1:1 | 20 v 20 | 0.050 | 0.045 | 0.005 (3.8) | 0.007 | 0.014 |
| 1:1 | 25 v 25 | 0.053 | 0.050 | 0.003 (2.9) | 0.006 | 0.012 |
| 1:1 | 30 v 30 | 0.050 | 0.046 | 0.004 (3.4) | 0.006 | 0.012 |
| 1:1 | 40 v 40 | 0.049 | 0.044 | 0.005 (4.0) | 0.007 | 0.011 |
| 1:1 | 50 v 50 | 0.054 | 0.053 | 0.001 (0.7) | 0.005 | 0.010 |
| 1:2 | 5 v 10 | 0.068 | 0.036 | 0.032 (11.2) | 0.032 | 0.039 |
| 1:2 | 7 v 14 | 0.059 | 0.041 | 0.018 (8.1) | 0.019 | 0.024 |
| 1:2 | 10 v 20 | 0.066 | 0.056 | 0.010 (5.6) | 0.012 | 0.016 |
| 1:2 | 15 v 30 | 0.054 | 0.049 | 0.005 (3.5) | 0.009 | 0.012 |
| 1:2 | 20 v 40 | 0.055 | 0.050 | 0.005 (4.1) | 0.006 | 0.011 |
| 1:2 | 25 v 50 | 0.052 | 0.048 | 0.003 (2.5) | 0.007 | 0.011 |
| 1:2 | 30 v 60 | 0.051 | 0.049 | 0.002 (1.5) | 0.005 | 0.010 |
| 1:2 | 40 v 80 | 0.053 | 0.050 | 0.003 (2.1) | 0.007 | 0.010 |
| 1:2 | 50 v 100 | 0.051 | 0.050 | 0.001 (0.6) | 0.006 | 0.010 |
| 2:1 | 10 v 5 | 0.064 | 0.039 | 0.025 (10.0) | 0.025 | 0.038 |
| 2:1 | 14 v 7 | 0.063 | 0.046 | 0.017 (7.7) | 0.019 | 0.023 |
| 2:1 | 20 v 10 | 0.056 | 0.048 | 0.008 (5.3) | 0.010 | 0.016 |
| 2:1 | 30 v 15 | 0.058 | 0.054 | 0.004 (3.3) | 0.007 | 0.012 |
| 2:1 | 40 v 20 | 0.049 | 0.047 | 0.002 (1.7) | 0.005 | 0.011 |
| 2:1 | 50 v 25 | 0.052 | 0.050 | 0.002 (1.5) | 0.006 | 0.011 |
| 2:1 | 60 v 30 | 0.054 | 0.050 | 0.004 (3.4) | 0.006 | 0.011 |
| 2:1 | 80 v 40 | 0.048 | 0.045 | 0.003 (2.3) | 0.007 | 0.011 |
| 2:1 | 100 v 50 | 0.053 | 0.052 | 0.001 (0.7) | 0.005 | 0.010 |

#### normal_sd3_vs_sd1

| ratio | n (x v y) | t | perm | t − perm (z) | disagree | mean \|p_t − p_perm\| |
|---|---|---|---|---|---|---|
| 1:1 | 5 v 5 | 0.085 | 0.054 | 0.032 (11.2) | 0.032 | 0.025 |
| 1:1 | 7 v 7 | 0.054 | 0.058 | -0.004 (-1.9) | 0.013 | 0.022 |
| 1:1 | 10 v 10 | 0.059 | 0.060 | -0.001 (-0.7) | 0.008 | 0.017 |
| 1:1 | 15 v 15 | 0.057 | 0.056 | 0.001 (0.6) | 0.007 | 0.014 |
| 1:1 | 20 v 20 | 0.050 | 0.050 | 0.000 (0.0) | 0.007 | 0.012 |
| 1:1 | 25 v 25 | 0.046 | 0.048 | -0.001 (-0.9) | 0.008 | 0.011 |
| 1:1 | 30 v 30 | 0.056 | 0.058 | -0.002 (-1.5) | 0.005 | 0.011 |
| 1:1 | 40 v 40 | 0.045 | 0.045 | -0.000 (-0.2) | 0.005 | 0.010 |
| 1:1 | 50 v 50 | 0.045 | 0.047 | -0.002 (-1.8) | 0.006 | 0.010 |
| 1:2 | 5 v 10 | 0.082 | 0.079 | 0.004 (2.8) | 0.007 | 0.029 |
| 1:2 | 7 v 14 | 0.049 | 0.056 | -0.007 (-4.3) | 0.012 | 0.021 |
| 1:2 | 10 v 20 | 0.057 | 0.062 | -0.005 (-3.2) | 0.009 | 0.016 |
| 1:2 | 15 v 30 | 0.054 | 0.059 | -0.005 (-3.7) | 0.007 | 0.013 |
| 1:2 | 20 v 40 | 0.046 | 0.048 | -0.002 (-1.7) | 0.005 | 0.012 |
| 1:2 | 25 v 50 | 0.052 | 0.055 | -0.003 (-2.6) | 0.006 | 0.011 |
| 1:2 | 30 v 60 | 0.050 | 0.052 | -0.002 (-2.2) | 0.004 | 0.010 |
| 1:2 | 40 v 80 | 0.054 | 0.058 | -0.004 (-2.9) | 0.006 | 0.010 |
| 1:2 | 50 v 100 | 0.051 | 0.051 | -0.001 (-0.6) | 0.006 | 0.010 |
| 2:1 | 10 v 5 | 0.058 | 0.034 | 0.024 (9.8) | 0.024 | 0.022 |
| 2:1 | 14 v 7 | 0.049 | 0.038 | 0.011 (6.6) | 0.011 | 0.016 |
| 2:1 | 20 v 10 | 0.048 | 0.039 | 0.010 (6.0) | 0.010 | 0.013 |
| 2:1 | 30 v 15 | 0.054 | 0.049 | 0.005 (3.9) | 0.006 | 0.011 |
| 2:1 | 40 v 20 | 0.054 | 0.050 | 0.004 (3.4) | 0.006 | 0.010 |
| 2:1 | 50 v 25 | 0.043 | 0.042 | 0.001 (1.3) | 0.005 | 0.010 |
| 2:1 | 60 v 30 | 0.049 | 0.048 | 0.001 (0.9) | 0.005 | 0.010 |
| 2:1 | 80 v 40 | 0.048 | 0.046 | 0.002 (1.9) | 0.006 | 0.010 |
| 2:1 | 100 v 50 | 0.048 | 0.048 | -0.000 (-0.2) | 0.006 | 0.010 |

#### skewed_scale3_vs_1

| ratio | n (x v y) | t | perm | t − perm (z) | disagree | mean \|p_t − p_perm\| |
|---|---|---|---|---|---|---|
| 1:1 | 5 v 5 | 0.091 | 0.059 | 0.033 (11.4) | 0.033 | 0.025 |
| 1:1 | 7 v 7 | 0.061 | 0.058 | 0.004 (2.0) | 0.014 | 0.022 |
| 1:1 | 10 v 10 | 0.051 | 0.052 | -0.000 (-0.4) | 0.005 | 0.017 |
| 1:1 | 15 v 15 | 0.060 | 0.059 | 0.002 (1.5) | 0.006 | 0.014 |
| 1:1 | 20 v 20 | 0.052 | 0.054 | -0.002 (-1.5) | 0.007 | 0.012 |
| 1:1 | 25 v 25 | 0.049 | 0.050 | -0.000 (-0.4) | 0.006 | 0.011 |
| 1:1 | 30 v 30 | 0.049 | 0.050 | -0.000 (-0.5) | 0.004 | 0.011 |
| 1:1 | 40 v 40 | 0.048 | 0.050 | -0.002 (-2.0) | 0.005 | 0.011 |
| 1:1 | 50 v 50 | 0.050 | 0.050 | 0.000 (0.2) | 0.004 | 0.010 |
| 1:2 | 5 v 10 | 0.079 | 0.074 | 0.005 (2.8) | 0.013 | 0.029 |
| 1:2 | 7 v 14 | 0.058 | 0.061 | -0.003 (-2.6) | 0.007 | 0.021 |
| 1:2 | 10 v 20 | 0.057 | 0.065 | -0.007 (-4.5) | 0.011 | 0.016 |
| 1:2 | 15 v 30 | 0.054 | 0.056 | -0.002 (-1.6) | 0.006 | 0.013 |
| 1:2 | 20 v 40 | 0.045 | 0.046 | -0.000 (-0.2) | 0.005 | 0.012 |
| 1:2 | 25 v 50 | 0.055 | 0.057 | -0.002 (-1.8) | 0.006 | 0.011 |
| 1:2 | 30 v 60 | 0.052 | 0.055 | -0.002 (-2.4) | 0.004 | 0.011 |
| 1:2 | 40 v 80 | 0.051 | 0.051 | -0.001 (-0.6) | 0.006 | 0.010 |
| 1:2 | 50 v 100 | 0.050 | 0.051 | -0.001 (-1.1) | 0.005 | 0.010 |
| 2:1 | 10 v 5 | 0.050 | 0.029 | 0.021 (9.2) | 0.021 | 0.022 |
| 2:1 | 14 v 7 | 0.055 | 0.043 | 0.012 (6.6) | 0.012 | 0.016 |
| 2:1 | 20 v 10 | 0.051 | 0.040 | 0.011 (6.5) | 0.011 | 0.013 |
| 2:1 | 30 v 15 | 0.052 | 0.048 | 0.004 (3.4) | 0.007 | 0.011 |
| 2:1 | 40 v 20 | 0.052 | 0.048 | 0.003 (3.3) | 0.004 | 0.011 |
| 2:1 | 50 v 25 | 0.051 | 0.048 | 0.004 (3.3) | 0.005 | 0.010 |
| 2:1 | 60 v 30 | 0.048 | 0.044 | 0.004 (3.6) | 0.005 | 0.010 |
| 2:1 | 80 v 40 | 0.058 | 0.055 | 0.003 (2.6) | 0.005 | 0.010 |
| 2:1 | 100 v 50 | 0.048 | 0.045 | 0.004 (2.8) | 0.007 | 0.010 |

## Caveats

- **Type I error only.** Power was not studied. With equivalent Type I error and about
  0.5% decision disagreement, power should also be nearly identical at larger n, but this
  was not checked.
- **Two-sample only.** The paired Brunner-Munzel test and the `logit` method were not
  studied here. An earlier run (`junk/sim_brunner_munzel_small.R`) found `logit` very
  conservative (.004–.030) at n ≤ 12. That run predates the floating point fix, and its
  `perm` ordinal cells should be rerun before being cited.
- **Two-sided α = .05 only.** Tail behaviour at smaller α (e.g., .005) typically needs
  larger n for the t-approximation to be adequate.
- **Floors from Monte Carlo error.** The disagreement (about 0.5%) and mean
  |p_t − p_perm| (about .010) floors are set by R = 999 in `perm`. Larger R would lower
  them.
- **Merged runs.** Continuous-data cells with n > 10 come from the run before the fix
  (see above). That is justified because exact ties essentially never occur there, but it
  is not identical code.
- **The "fairly large" message.** `brunner_munzel()` currently suggests that `perm` is
  unnecessary only when both groups are large (min n > 250 and total > 500). These
  results suggest a much lower point, about 20–30 per group, would be defensible. The
  current threshold is conservative rather than wrong.

## Reproduction

```r
# Full grid (original run; about 1 hour on 13 cores)
Sys.setenv(SIM_NSIM = 4000)
source("junk/sim_brunner_munzel_t_vs_perm_n.R")

# Reruns after the floating point fix
Sys.setenv(SIM_NSIM = 4000, SIM_NMAX = 10,
           SIM_OUT = "junk/sim_bm_t_vs_perm_n_rerun_small.rds")
source("junk/sim_brunner_munzel_t_vs_perm_n.R")
Sys.unsetenv("SIM_NMAX")
Sys.setenv(SIM_DISTS = "identical_ordinal",
           SIM_OUT = "junk/sim_bm_t_vs_perm_n_rerun_ordinal.rds")
source("junk/sim_brunner_munzel_t_vs_perm_n.R")

# Merge and print tables
source("junk/bm_t_vs_perm_merge.R")
```

Files:

- `junk/sim_brunner_munzel_t_vs_perm_n.R`: simulation script
- `junk/bm_t_vs_perm_merge.R`: merge and tables
- `junk/sim_brunner_munzel_t_vs_perm_n_results.rds`: original run
- `junk/sim_bm_t_vs_perm_n_rerun_small.rds`, `junk/sim_bm_t_vs_perm_n_rerun_ordinal.rds`:
  reruns
- `junk/sim_brunner_munzel_t_vs_perm_n_merged.rds`: merged results used above

## References

Neubert, K., & Brunner, E. (2007). A studentized permutation test for the non-parametric
Behrens–Fisher problem. *Computational Statistics & Data Analysis*, 51(10), 5192–5204.
