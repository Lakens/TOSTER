# perm_t_test: is a shifted permutation null worth adding?

*Simulation study, 2026-09-29. Script: `junk/sim_perm_t_test_shifted_null.R`.
Results: `junk/sim_perm_t_test_shifted_null_results.rds`.*

## Question

`perm_t_test()` tests a null value θ (e.g., an equivalence bound) by comparing the
observed statistic (d̂ − θ)/SE with the permutation distribution of the **unshifted**
data. A **shifted** null (as in Arboretti, Pesarin & Salmaso, 2021) instead permutes data
shifted so that θ becomes zero:
- two-sample: permute pooled x − θ and y
- paired / one-sample: flip the signs of d − θ

The shifted test is exact at θ under a pure location shift (two-sample), or under
symmetry of the differences about θ (paired). The unshifted test is only asymptotically
valid there.

Would a `perm_null = "shifted"` option improve TOST and minimal-effect testing enough to
justify the cost? The main cost is that CIs would need numerical root-finding to stay
consistent with the p-values.

## Short answer: no, the juice isn't worth the squeeze

1. **Where shifting is exact (equal shapes), the current unshifted test is already at
   nominal.** This holds even at n = 6 per group, with skewed or heavy-tailed data, and
   with bounds as wide as 2 SD. Over all equal-shape two-sample cells, TOST Type I error
   averaged .049 (unshifted) vs .049 (shifted), and minimal-effect .050 vs .050. The
   studentized statistic is close enough to pivotal under a location shift that the
   unshifted reference distribution works.
2. **Where the current test fails (unequal variances, small n), shifting doesn't help.**
   It is slightly worse for TOST (.065–.081 vs .050–.076) and slightly better for
   minimal-effect tests (.061–.076 vs .065–.081). Both are liberal, while Welch stays near
   .05.
3. **Paired data:**
   - Symmetric (normal) differences: both versions are nominal.
   - Skewed differences: all methods, including Welch, over-reject (.09–.16), and shifting
     changes nothing. The sign-flip needs symmetry about θ, which a skewed distribution
     doesn't have.
4. **Power is unchanged.** Apart from the heteroskedastic cells, where "shifted" rejects
   more only because it is more liberal, TOST power differs by at most about 0.01–0.02.

So the shifted null would buy a theoretical guarantee in a setting where the current test
already performs at nominal. It would do nothing, or slightly harm, where the current test
has a real problem. It would also cost the closed-form, always-consistent CI we just
implemented for issue #120.

## Design

- **Methods:** the three are applied to the same datasets, and both permutation versions
  use the same permutation indices.
  - unshifted: the current `perm_t_test`. The prototype matched the package's p-values
    exactly in fully enumerated cases; maximum difference 0.
  - shifted: the proposed option.
  - Welch t-test: reference.
- **Truth placed exactly at a bound (θ = −B or +B):**
  - TOST Type I error at that bound = the one-sided test in the equivalence direction
  - minimal-effect Type I error = the one-sided test in the opposite direction
  - their sum = 1 − coverage of the 90% CI at the bound
  
  Truth at 0 gives TOST power. Tables report the worse (larger) of the two bounds.
- **Two-sample, equal shape (location shift; shifting is exact here):**
  - distributions: normal, skewed (Exp(1) − 1) and heavy-tailed (t₃/√3), all with SD 1
  - sample sizes: 6v6 (all 924 permutations enumerated), 10v10 and 20v20
  - bounds: B = 0.5, 1 and 2 SD
- **Two-sample, heteroskedastic (not exact for either method):** normal, the smaller group
  has SD 3, at 8v16 and 16v32, with B = 0.5, 1 and 2.
- **Paired (one-sample on differences):** normal and skewed differences (SD 1), n = 8 (all
  256 sign patterns enumerated), 15 and 30, with B = 0.5 and 1.
- **Replications:** 4,000 per cell, R = 999 permutations, α = 0.05 (one-sided tests,
  i.e. the 90% CI). MC SE ≈ .0034 per cell. Taking the larger of two bounds biases the
  table values slightly upward.

## Results

### Two-sample, equal shape: both at nominal

Mean over the 27 equal-shape cells (each bound separately; SE of a mean ≈ .0005–.0008):

| cells | TOST unshifted | TOST shifted | TOST Welch | MET unshifted | MET shifted | MET Welch |
|---|---|---|---|---|---|---|
| 6v6 (enumerated) | .050 | .050 | .043 | .049 | .050 | .043 |
| 10v10, 20v20 | .049 | .049 | .047 | .050 | .050 | .048 |

- Worst single cell (larger of the two bounds):
  - TOST: unshifted .0565 (heavy, 10v10, B = 0.5); shifted .0575 (skewed, 6v6, B = 2)
  - minimal-effect: unshifted .058; shifted .0565
  
  All are consistent with MC noise around .05 after taking the maximum of two bounds.
- Welch is slightly conservative with skewed or heavy-tailed data in small samples
  (down to .036–.044).

### Two-sample, heteroskedastic (small group SD 3): shifting doesn't help

| n (x v y) | B | TOST unshifted | TOST shifted | TOST Welch | MET unshifted | MET shifted | MET Welch |
|---|---|---|---|---|---|---|---|
| 8 v 16 | 0.5 | .067 | .071 | .049 | .078 | .076 | .052 |
| 8 v 16 | 1 | .076 | .081 | .057 | .078 | .075 | .053 |
| 8 v 16 | 2 | .060 | .074 | .053 | .081 | .076 | .053 |
| 16 v 32 | 0.5 | .062 | .067 | .051 | .076 | .074 | .061 |
| 16 v 32 | 1 | .063 | .072 | .054 | .068 | .064 | .049 |
| 16 v 32 | 2 | .050 | .065 | .053 | .065 | .061 | .051 |

Both permutation versions are liberal here. This is the weak-null problem: pooling the
groups under permutation mixes the variances, and shifting moves locations but not
variances. Welch is the best of the three in this setting.

### Paired / one-sample differences

| differences | n | TOST unshifted | TOST shifted | TOST Welch |
|---|---|---|---|---|
| normal | 8 | .052 | .049 | .052 |
| normal | 15 | .048 | .048 | .049 |
| normal | 30 | .056 | .057 | .056 |
| skewed | 8 | .156 | .147 | .153 |
| skewed | 15 | .130 | .127 | .127 |
| skewed | 30 | .097 | .096 | .096 |

(Worst over B = 0.5 and 1. Minimal-effect results follow the same pattern.)

With skewed differences, every mean-based method over-rejects at the bounds in small
samples. Shifting doesn't help, because the shifted sign-flip requires symmetry about θ.

### Power (truth at 0)

TOST power is essentially identical for unshifted and shifted, within about 0.01–0.02 in
all equal-shape and paired cells. For example, 10v10 normal with B = 1 gives .394 vs .392,
and paired n = 15 normal with B = 1 gives .952 vs .951. In the heteroskedastic cells the
shifted version rejects slightly more often (e.g. .262 vs .190 at 8v16 with B = 2), which
only reflects its higher Type I error.

## Recommendation

- **Don't add `perm_null`.** In the only setting where it is exact, it improves nothing
  measurable. In the problem setting (heteroskedastic small samples), it is no better and
  sometimes worse. It would also require replacing the closed-form, always-consistent CI
  with numerical root-finding over the null value, which is slower and needs handling for
  non-monotone p-value functions.
- **Good news for the docs and the Avocado paper.** The current `perm_t_test` TOST is
  already well calibrated at the bounds under location-shift models, including skewed and
  heavy-tailed data at n = 6 per group. The "approximate rather than exact" caveat in
  `?perm_t_test` is still technically correct, but in practice the approximation is
  excellent in these settings. That could be stated.
- **The real weaknesses, which shifting doesn't address and which are worth documenting
  for TOST specifically:**
  1. Heteroskedastic small samples where the smaller group is more variable: permutation
     TOST and minimal-effect Type I error are about .06–.08. Welch (`t_TOST`) holds about
     .05.
  2. Skewed paired differences in small samples: all mean-based methods reach about
     .10–.16 at the bounds. Untested here: whether `boot_t_test` does better, which the
     assumption-check supplement already recommends when differences are asymmetric.
     That would be a worthwhile follow-up check.
- **If Arboretti-style calibration is ever revisited:** it targets a different issue,
  TOST's inherent conservativeness with narrow margins (visible here as very low power,
  e.g. .001–.03 for B = 0.5 with n ≤ 10). It should be evaluated separately, and its
  conflict with the CI–TOST duality weighed first.
