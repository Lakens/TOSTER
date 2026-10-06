# Type I error: boot_cor_test and boot_smd_calc (by `boot_ci`) vs analytic methods

Scripts: `junk/sim_boot_cor_test_type1.R`, `junk/sim_boot_smd_calc_type1.R`
Results: `junk/sim_boot_cor_test_type1_results.rds`, `junk/sim_boot_smd_calc_type1_results.rds`

All runs use `alpha = 0.05` and bootstrap `R = 999`. Within each dataset, every `boot_ci` type is computed from the same bootstrap resamples (the RNG state is restored before each call).

## Correlations (Pearson): `boot_cor_test` vs `z_cor_test`

> **Superseded (2026-10-06):** these results use the pre-`d25c85c` `boot_cor_test`, whose Pearson "stud" used a constant Fisher-z SE (basic-on-z). Upstream commit `d25c85c` replaced it with a per-resample sandwich SE and made "stud" the Pearson default via `boot_ci = "auto"`. The "keep bca" recommendation below no longer applies.

nsim = 5000 per cell, so the MC SE is ~0.003 near 0.05; a cell is notable when it falls outside about 0.044–0.056.

The DGP is `y = b*x + e`, with x and e standardized:

- **normal**: x, e ~ N(0,1)
- **skewed**: x ~ Exp(1) − 1, e ~ t(4)/√2
- **hetero**: x ~ N(0,1), e = z·√((1+x²)/2). At rho = 0, x and y are dependent but uncorrelated, which is the classic case where the Fisher-z SE is wrong.

### Two-sided test at the true value (rho = 0)

|  dist  |  n |  bca  |  stud | basic |  perc | fisher_z |
|--------|---:|------:|------:|------:|------:|---------:|
| normal | 10 | 0.056 | 0.036 | 0.217 | 0.091 | 0.048 |
| normal | 20 | 0.076 | 0.067 | 0.159 | 0.085 | 0.061 |
| normal | 50 | 0.063 | 0.059 | 0.095 | 0.068 | 0.053 |
| skewed | 10 | 0.051 | 0.043 | 0.233 | 0.085 | 0.058 |
| skewed | 20 | 0.061 | 0.062 | 0.155 | 0.080 | 0.049 |
| skewed | 50 | 0.072 | 0.069 | 0.106 | 0.069 | 0.050 |
| hetero | 10 | 0.094 | 0.080 | 0.334 | 0.115 | 0.127 |
| hetero | 20 | 0.102 | 0.107 | 0.245 | 0.098 | 0.150 |
| hetero | 50 | 0.082 | 0.085 | 0.147 | 0.076 | 0.155 |

### Equivalence at the bound (rho = 0.3, bounds ±0.3)

These are one-sided rejection rates at the near bound (90% CI upper limit < 0.3), which should be about 0.05. The TOST decision itself is near 0 for n ≤ 20 with these bounds for every method, because the far bound can't be rejected.

|  dist  |  n |  bca  |  stud | basic |  perc | fisher_z |
|--------|---:|------:|------:|------:|------:|---------:|
| normal | 10 | 0.053 | 0.050 | 0.088 | 0.064 | 0.048 |
| normal | 20 | 0.058 | 0.054 | 0.063 | 0.058 | 0.048 |
| normal | 50 | 0.059 | 0.059 | 0.053 | 0.057 | 0.051 |
| skewed | 10 | 0.049 | 0.059 | 0.093 | 0.056 | 0.047 |
| skewed | 20 | 0.066 | 0.075 | 0.081 | 0.061 | 0.058 |
| skewed | 50 | 0.074 | 0.083 | 0.074 | 0.060 | 0.062 |
| hetero | 10 | 0.066 | 0.075 | 0.126 | 0.070 | 0.091 |
| hetero | 20 | 0.069 | 0.076 | 0.098 | 0.062 | 0.105 |
| hetero | 50 | 0.066 | 0.066 | 0.064 | 0.058 | 0.111 |

TOST rejection rate at n = 50: bca/stud/basic/perc ≈ 0.045–0.066 in the normal and skewed cells and ≈ 0.017 in the hetero cell. Fisher z is 0.080 in the hetero cell.

### Takeaways for the `boot_cor_test` default

1. **The data don't support switching `"bca"` to `"stud"`.** For Pearson, the two are within MC error of each other in most cells. Stud is more conservative at n = 10 under normality (0.036 vs 0.056). It is a little more liberal than BCa at the bound with skewed data (0.083 vs 0.074 at n = 50). Neither is consistently better.
2. **Pearson "stud" isn't really studentized.** `.fisher_z_se()` returns the constant 1/√(n−3) for Pearson, so the pivot `(z* − z_obs)/se` is just the basic bootstrap on the Fisher-z scale. That explains why it tracks BCa (both respect the transformation) and is far better than the raw-scale `"basic"`. It also means it can't adapt to heteroscedasticity the way a real studentized bootstrap (with a per-resample SE) would. The hetero cells reflect this, at ~0.08–0.11.
3. **`"basic"` on the raw r scale is badly liberal** (0.22–0.33 at n = 10) and should never be the default.
4. **All bootstrap variants run ~0.06–0.08 at n = 20–50, even under bivariate normality**, where Fisher z is near nominal. This is a general bootstrap-for-r issue, not a BCa issue.
5. **Under heteroscedasticity, Fisher z fails and the failure doesn't shrink with n** (0.15 at n = 50) because its SE is wrong. The bootstrap is clearly better there (0.08 and decreasing). This is the main reason to use `boot_cor_test`.
6. **p/CI agreement:** p-values and CIs disagree in 0.2–0.7% of datasets for bca, stud and perc (basic and Fisher z always agree). These look like boundary/discreteness cases, but it's worth a look given this branch's issue.

## SMDs (two independent groups): `boot_smd_calc` vs `smd_calc`

nsim = 500 per cell, so the MC SE is ~0.0097 near 0.05; a cell is notable only when it falls outside about 0.03–0.07. These results are coarser than the correlation results.

Setup: equal n per group, function defaults (Hedges' g_av, var.equal = FALSE). Data are standardized draws scaled by group SDs, so the true d_av is exact.

- **normal**: N(·,1) in both groups
- **skewed**: standardized lognormal(0, 0.8) in both groups (same shape, shifted)
- **hetero**: normal, SDs 1 vs 3

Comparators: `smd_calc` with smd_ci/test_method set to z/z, t/t, and nct/z. For nct the p-value comes from the z test, so p and CI disagree in ~1% of datasets. That is expected.

**Failures at n = 10:** `boot_smd_calc` gave NA p-values or errors in 18–27 of 500 reps per cell (~4–5%, 136 per method in total). The rates below are among the reps that didn't fail. The cause is that the two-sample bootstrap resamples the pooled data (`R/boot_smd_calc.R:396`), so group sizes vary and sometimes drop to ≤1, which gives a NaN SMD. `boot_ses_calc` has the same pattern. The fix is within-group resampling, as in `boot_t_test`; this has been flagged as a separate task.

### Two-sided test at the true value (d = 0)

|  dist  |  n |  stud | basic |  perc |  bca  | smd_z | smd_t | smd_nct |
|--------|---:|------:|------:|------:|------:|------:|------:|--------:|
| normal | 10 | 0.050 | 0.017 | 0.066 | 0.039 | 0.038 | 0.020 | 0.038 |
| normal | 20 | 0.048 | 0.020 | 0.058 | 0.038 | 0.042 | 0.036 | 0.042 |
| normal | 50 | 0.040 | 0.038 | 0.046 | 0.042 | 0.042 | 0.036 | 0.042 |
| skewed | 10 | 0.072 | 0.027 | 0.066 | 0.040 | 0.016 | 0.008 | 0.016 |
| skewed | 20 | **0.124** | **0.112** | 0.078 | 0.060 | 0.048 | 0.034 | 0.048 |
| skewed | 50 | **0.084** | 0.076 | 0.052 | 0.036 | 0.030 | 0.026 | 0.030 |
| hetero | 10 | 0.042 | 0.011 | 0.074 | 0.032 | 0.032 | 0.026 | 0.032 |
| hetero | 20 | 0.036 | 0.022 | 0.062 | 0.028 | 0.034 | 0.024 | 0.034 |
| hetero | 50 | 0.038 | 0.034 | 0.050 | 0.036 | 0.040 | 0.040 | 0.040 |

### Equivalence at the bound (d = 0.5, bounds ±0.5)

One-sided rejection rate at the near bound (90% CI upper < 0.5):

|  dist  |  n |  stud | basic |  perc |  bca  | smd_z | smd_t | smd_nct |
|--------|---:|------:|------:|------:|------:|------:|------:|--------:|
| normal | 10 | 0.050 | 0.042 | 0.050 | 0.042 | 0.042 | 0.030 | 0.050 |
| normal | 20 | 0.030 | 0.030 | 0.036 | 0.032 | 0.030 | 0.026 | 0.032 |
| normal | 50 | 0.042 | 0.042 | 0.040 | 0.040 | 0.038 | 0.038 | 0.040 |
| skewed | 10 | 0.065 | 0.063 | 0.013 | 0.023 | 0.022 | 0.020 | 0.036 |
| skewed | 20 | **0.086** | **0.086** | 0.040 | 0.046 | 0.062 | 0.054 | 0.070 |
| skewed | 50 | 0.070 | 0.068 | 0.022 | 0.046 | 0.062 | 0.060 | 0.068 |
| hetero | 10 | 0.059 | 0.040 | 0.055 | 0.040 | 0.040 | 0.030 | 0.060 |
| hetero | 20 | 0.046 | 0.042 | 0.050 | 0.046 | 0.046 | 0.046 | 0.048 |
| hetero | 50 | 0.048 | 0.042 | 0.046 | 0.050 | 0.048 | 0.046 | 0.048 |

TOST rejection is 0 for n ≤ 20 (the bounds are too narrow relative to the SE). At n = 50 every method is 0.02–0.06.

### Takeaways for `boot_smd_calc`

1. **Under normality or unequal variances, everything is at or below nominal.** The analytic `smd_calc` is fine there and slightly conservative.
2. **With skewed data, the current `"stud"` default is the most liberal option** (0.124 at n = 20 and 0.084 at n = 50, two-sided; 0.086 at the bound). **BCa stays within MC error of nominal in every cell** (0.023–0.060). The likely reason is that the stud pivot divides by the normal-theory SE of d. That SE ignores skewness and kurtosis, so the pivot isn't pivotal under non-normal data, unlike the t-statistic pivot in `boot_t_test`.
3. So the t-test finding (prefer stud) doesn't carry over to SMDs. If anything, these results point toward BCa for `boot_smd_calc`. That should be confirmed with a larger nsim after the resampling fix, since the pooled resampling also changes the bootstrap distribution.

## Bottom line

- **`boot_cor_test`: keep `"bca"`.** For Pearson, `"stud"` is basic-on-Fisher-z and performs about the same as BCa, so switching would gain nothing. If a real studentized option is wanted, it would need a per-resample, non-normal-theory SE (e.g., a jackknife or delta-method SE of z within each resample).
- **`boot_smd_calc`: consider switching the default from `"stud"` to `"bca"`**, after fixing the two-sample resampling bug and rerunning with nsim ≥ 2000.

## Addendum: skewed SMD cells rerun after the resampling fix

Within-group resampling, nsim = 1000 (MC SE ~0.007), results in `junk/sim_boot_smd_calc_type1_skewed_fixed_results.rds`. 0 failed calls.
Two-sided rows show the rejection rate at d = 0. Equivalence rows show the one-sided rejection rate at the near bound (d = 0.5).

|  hyp        |  n |  stud | basic |  perc |  bca  | smd_z | smd_t | smd_nct |
|-------------|---:|------:|------:|------:|------:|------:|------:|--------:|
| two-sided   | 10 | 0.080 | 0.038 | 0.071 | 0.042 | 0.014 | 0.008 | 0.014 |
| two-sided   | 20 | **0.114** | 0.091 | 0.071 | 0.048 | 0.030 | 0.022 | 0.030 |
| two-sided   | 50 | **0.104** | 0.094 | 0.060 | 0.047 | 0.039 | 0.036 | 0.039 |
| equivalence | 10 | 0.090 | 0.084 | 0.022 | 0.029 | 0.026 | 0.014 | 0.052 |
| equivalence | 20 | 0.073 | 0.070 | 0.019 | 0.033 | 0.042 | 0.037 | 0.055 |
| equivalence | 50 | 0.060 | 0.060 | 0.018 | 0.036 | 0.051 | 0.050 | 0.055 |

The bug fix removed the failures but did **not** fix stud. Stud stays at ~0.10–0.11 and isn't improving with n. BCa is at or below nominal in every cell. The likely cause is the normal-theory SE inside the pivot (no skewness or kurtosis terms).
