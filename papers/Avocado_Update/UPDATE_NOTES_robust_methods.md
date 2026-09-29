# Avocado_Update.Rmd: suggested updates after the permutation/robust-methods changes

*Prepared 2026-09-29 against branch `120-a-potential-lack-of-alignment-between-the-ci-and-p-value-from-the-permutation-welch-t-test`.
Line numbers refer to `Avocado_Update.Rmd` as of this date.*

## What changed in the package (context for the edits)

1. **`perm_t_test` confidence intervals (issue #120).**
   - **Before:** the CI was a percentile interval of the raw permuted mean differences, while the p-value came from the studentized (Welch) test. The two could disagree under unequal variances.
   - **Now:** the CI inverts the same studentized test, `estimate − SE · q`. So the CI and the p-value always agree, for all alternatives including equivalence (the 90% CI) and minimal-effect tests.
   - **Scope:** the same fix was applied to `hodges_lehmann()` and `brunner_munzel(test_method = "perm")`.
2. **Floating point tie fix.** Permutation p-values could be too small with tied data. The exact Brunner–Munzel test with 5 vs 5 ordinal data rejected 8% of the time under the null.
3. **New `brunner_munzel(test_method = "perm_logit")`.** A studentized permutation test on the logit scale, with range-preserving CIs that always agree with its p-value.
4. **New warnings.** `brunner_munzel` now warns when `"t"` or `"perm"` is used for equivalence, minimal-effect, or any null other than 0.5.
5. **New docs on sharp vs. weak null hypotheses** in `?perm_t_test` and `?boot_t_test` (citing Wu & Ding, 2020). New "Choosing a test method" guidance in `?brunner_munzel`.

The simulation evidence behind these changes is in `junk/` (not part of the package). The main reports are:
- `junk/bm_t_vs_perm_sample_size.md` (two-sample Brunner–Munzel, `t` vs `perm`)
- `junk/bm_paired_t_vs_perm_sample_size.md` (paired Brunner–Munzel)
- `junk/sim_perm_t_test_ci.R` and `junk/sim_perm_t_test_type1.R` (`perm_t_test`)
- `junk/sim_brunner_munzel_tost.R` (Brunner–Munzel equivalence and minimal-effect tests)

---

## Priority 1: statements that are now wrong or would mislead

### 1.1 Choosing between bootstrap and permutation (Resampling Methods, line 787)

This is the most important paragraph to revise. Several of its contrasts no longer hold.

- **"[bootstrap] provides readily interpretable confidence intervals via CI inversion ... and ensures consistency between reported intervals and p-values."** This is no longer a bootstrap-only advantage. `perm_t_test` intervals are now obtained by inverting the permutation test and always agree with the p-value.
- **"In small samples, permutation approaches tend to have better finite-sample calibration because the null distribution is constructed exactly..."** This is true only under the **sharp null** (identical distributions). Under the **weak null** (equal means, different variances or shapes), no permutation test is exact. Simulation results (`junk/sim_perm_t_test_type1.R`, 1,000 reps, MC SE ≈ .007):

  | scenario | Welch | perm (R = 999) | perm (full enumeration) |
  |---|---|---|---|
  | identical normals, 7 v 7 | .046 | .047 | .048 |
  | n = 5 (SD 4) v n = 10 (SD 1) | .051 | **.086** | **.087** |
  | n = 5 (SD 1) v n = 10 (SD 4) | .052 | .041 | .038 |
  | equal n, SD 1 vs 4 (7 v 7) | .042 | .060 | .060 |
  | skewed, lognormal 5 v exponential 10 | .136 | .139 | .135 |

  The studentized permutation test can be liberal when the **smaller group is the more variable one** and n is very small, where Welch holds. No method fixes the small-sample problem of comparing means of differently skewed distributions.
- **"researchers comparing two independent, randomized groups may find the studentized permutation test most defensible"** This is right, and now has a sharper justification (Wu & Ding, 2020):
  - In a randomized experiment the permutation distribution *is* the randomization distribution. So the studentized test is exact under the sharp null and asymptotically valid under the weak null, whatever the population.
  - With random samples from populations (observational data), the justification shifts to an exchangeability assumption, and the bootstrap, which targets the weak null directly, is an equally natural choice.
  - Neither supports population inference from convenience samples.
- **Missing, and specific to equivalence testing:** exactness under the sharp null applies only to the test of **no effect** (`mu = 0`).
  - TOST tests at the bounds, and `perm_t_test` compares the observed statistic (centered at each bound) against the permutation distribution of the **unshifted** data.
  - So a permutation TOST is approximate (asymptotically valid) rather than exact, even in a randomized experiment.
  - This weakens the "exact" argument specifically for the use case this paper is about, and should be stated.

**Suggested replacement for line 787** (edit freely):

> The choice between the bootstrap and permutation t-tests involves several tradeoffs. Both target the same estimand (the difference in means), and in TOSTER both return confidence intervals that are consistent with their p-values (the permutation interval is obtained by inverting the studentized permutation test). The key difference is what justifies each method. The permutation test is exact under the *sharp* null hypothesis that the two groups have identical distributions, at any sample size. When studentized, it is also asymptotically valid under the *weak* null hypothesis that only the means are equal [@janssen1997; @chung2013]. In a randomized experiment the permutation distribution is the randomization distribution itself, so this justification comes from the design rather than from assumptions about the population [@Wu_Ding_2020]. This makes the studentized permutation test a natural choice for randomized experiments, particularly small ones. The bootstrap targets the weak null directly and is only asymptotically valid (never exact). It is an equally natural choice when the data are random samples from two populations, and it extends to statistics beyond the mean difference, such as standardized effect sizes. Two caveats apply. First, validity under the weak null is only asymptotic for both methods. In very small samples with unequal variances the studentized permutation test can be somewhat liberal, particularly when the smaller group has the larger variance. Neither method solves the small-sample problem of comparing means of differently skewed distributions. Second, for equivalence testing the permutation test is evaluated at the equivalence bounds rather than at zero, so it is approximate rather than exact even under the sharp null. When outliers are a concern, both methods accept the `tr` argument for trimmed means. Neither approach supports population inference from non-random (e.g., convenience) samples.

Optionally follow with a short table:

| Situation | Suggested default |
|---|---|
| Small randomized two-group experiment, test of no difference | `perm_t_test` (studentized) |
| Random samples from populations, or interest in SMDs or other statistics | `boot_t_test` / `boot_t_TOST` |
| Paired data with clearly asymmetric differences | `boot_t_test` (sign-flipping requires symmetry) |
| Very small n with the smaller group clearly more variable | Welch `t_TOST` or `boot_t_test`; treat permutation results with caution |

### 1.2 Resampling Methods intro (lines 781 and 783)

- **Line 781:** "the permutation approach constructs an exact or near-exact null distribution". Add "under the null hypothesis of identical distributions (the sharp null)", or point forward to the new 1.1 paragraph.
- **Line 783:** "@Arboretti_Pesarin_Salmaso_2020 extended this logic to the equivalence testing setting, showing that the studentized permutation of TOST maintains nominal Type I error rates when variances differ across groups."
  - **Check this claim against the paper.** Arboretti et al. develop a nonparametric-combination approach that permutes **margin-shifted** data and discuss **calibration**.
  - TOSTER's `perm_t_test` permutes the **unshifted** data and is **uncalibrated**; the function docs say it is asymptotically equivalent and can be conservative in the equivalence direction.
  - Suggested rewording: "@Arboretti_Pesarin_Salmaso_2020 developed permutation approaches to equivalence testing; TOSTER's implementation is a simpler, uncalibrated variant that is asymptotically equivalent under studentization (see `?perm_t_test`)."

### 1.3 Brunner–Munzel example and text (lines 755–772)

- **Line 755:** "three test methods are available: `"t"` (default), `"logit"` ..., and `"perm"`". There are now **four**. Add `"perm_logit"` (studentized permutation on the logit scale; range-preserving intervals). Suggested addition:

  > For equivalence and minimal effect tests (or any null value other than 0.5), `"perm_logit"` or `"logit"` is recommended; in simulations the `"t"` and `"perm"` methods produced overly conservative equivalence tests and inflated Type I error for minimal effect tests, and the function now warns when they are used for these hypotheses.

- **Code chunk at lines 757–768.** This runs an **equivalence** test with the default `test_method = "t"`.
  - It now triggers the new warning. Readers won't see it, because the setup chunk sets `warning = FALSE`, but the example contradicts the package's recommendation.
  - In the simulations (bounds 0.2–0.8 through 0.4–0.6, 7–30 per group), `t` equivalence tests held their level but were very conservative. Power at p = 0.5 with bounds (0.2, 0.8) and 10 per group was .32 for `t`, .53 for `logit`, and .64 for `perm_logit`.
  - Suggested chunk:

    ```r
    set.seed(2071)
    bm_test = brunner_munzel(
      formula = LDHF ~ Gender,
      data = bugs_g,
      # equivalence hypothesis
      alternative = "equivalence",
      # equivalence bounds
      # on the stochastic superiority scale
      mu = c(0.3, 0.7),
      # permutation test on the logit scale
      # (recommended for equivalence testing)
      test_method = "perm_logit",
      R = 9999
    )
    print(bm_test)
    ```
  - Any prose interpreting this output (p-value, CI) will need re-checking after re-knitting.
  - **Checked on this branch** (61 female, 29 male ratings):

    | test_method | p (equivalence) | 90% CI | runtime |
    |---|---|---|---|
    | `"t"` (current paper) | 0.0275 | (0.475, 0.682) | instant |
    | `"perm_logit"`, R = 9999, seed 2071 | 0.0177 | (0.476, 0.675) | ~1 s |

    The estimate is 0.579 either way, and the conclusion (equivalence within 0.3–0.7) doesn't change.
- **Line 772 (Checking Assumptions):**
  - **Bug in the text:** `(method = "perm")` should be `(test_method = "perm")`.
  - **Grammar:** "regardless of whether the shapes of the distributions" → "regardless of the shapes of the distributions".
  - **Update the small-sample recommendation to match the new docs:**
    - null of 0.5, two-sample: `"perm"` when the smaller group has fewer than 30
    - paired: `"perm"` with at least 6 pairs (with 5 or fewer pairs the permutation test cannot reject at α = .05)
    - equivalence and minimal effect: `"perm_logit"` or `"logit"`
  - A one-sentence pointer to "Choosing a test method" in `?brunner_munzel` may be enough.
- **Line 751:** "the Brunner–Munzel test remains valid in that setting". Consider adding "asymptotically".
  - In small samples where the smaller group is more variable, both the `t` and `perm` versions can be liberal (about 8% at 5 v 10 in the simulations).
  - Minor, but it keeps the paper consistent with the new docs.

### 1.4 Permutation t-test section (lines 833–850)

- **Add how the interval is obtained**, parallel to the bootstrap section's treatment. Suggested addition after the first paragraph (line 833):

  > The confidence interval is obtained by inverting the same studentized permutation test used for the p-value, so the interval and the p-value always agree: the interval excludes the equivalence bounds if and only if the equivalence test rejects.^[Previous versions of the package reported a percentile interval of the permuted mean differences, which ignored the studentization and could disagree with the p-value when variances differed between groups.]

- **Line 850:** "The two approaches can, however, diverge in small samples, under strong heteroscedasticity..." This is fine. Optionally cross-reference the new 1.1 text.
- **Checked on this branch:** the paper's example (`sleep`, paired, bounds ±0.5, `set.seed(8812)`, R = 999) gives p = 0.999 with a 90% CI of (−2.26, −0.89). Compare with the currently rendered PDF and update any prose that quotes the old interval.

### 1.5 Permutation "Checking Assumptions" (line 854)

- **"the studentized version is valid under exchangeability ... and provides asymptotic validity even when variances differ"** This is fine, but add the two caveats from 1.1:
  - the small-sample liberal case (smaller group more variable)
  - permutation TOST at nonzero bounds is approximate rather than exact
- **The paired sign-flip symmetry discussion is accurate.**

---

## Priority 2: accuracy and consistency improvements

- **Estimand table (line 106):** the pseudo-median row lists only `wilcox_TOST`. `hodges_lehmann()` also targets the pseudo-median, with asymptotic or permutation inference (and the fixed CI). Consider adding it.
- **Estimand table (line 107):** consider noting `test_method = "perm_logit"` for equivalence, or leave the recommendation to the text.
- **Conclusions (line 1096):** the advice to compare Student's t vs permutation vs bootstrap still stands. Optionally add: "for equivalence tests of stochastic superiority, compare `"perm_logit"` and `"logit"`".

## Priority 3: typos (line 779)

- "botha" → "both a"
- "evulate" → "evaluate"
- Line 785: "The rank-based tests are useful, but as answer different questions" → "... useful, but as they answer different questions".

---

## Citations to add

Add to `interactcadsample.bib` (the key used in the suggested text is `Wu_Ding_2020`):

```bibtex
@article{Wu_Ding_2020,
  title={Randomization Tests for Weak Null Hypotheses in Randomized Experiments},
  volume={116},
  DOI={10.1080/01621459.2020.1750415},
  number={536},
  journal={Journal of the American Statistical Association},
  publisher={Informa UK Limited},
  author={Wu, Jason and Ding, Peng},
  year={2020},
  pages={1898--1913}
}
```

`janssen1997`, `chung2013`, `neubert2007`, `Arboretti_Pesarin_Salmaso_2020`, and `Thulin_2024` are already in the bib.

---

## Online supplement (`toster_assumption_checks_supplement.qmd`)

No code in the supplement uses permutation CIs (the replication-stability example at lines 268–278 uses `boot_t_test`, which is unchanged), so its **numbers are unaffected**. Three text points are worth aligning:

- **Line 348 and the summary table (line 777):** "provides asymptotic validity under unequal variances". Add the small-sample caveat (the smaller, more variable group can make the test liberal at very small n).
- **Line 611:** "Brunner-Munzel gives up a location-difference interpretation in exchange for validity under unequal shapes and variances". Add "asymptotic" (same caveat as 1.3).
- **"When to Use What" (line 711):** optionally add a bullet on choosing the Brunner–Munzel `test_method` (see 1.3), especially `"perm_logit"`/`"logit"` for equivalence.

---

## After editing

1. **Re-knit** `Avocado_Update.Rmd`. The Brunner–Munzel output changes if you switch to `perm_logit`, and the `perm_t_test` equivalence example's CI changes because of the CI fix.
   - The example's p-value is unchanged unless ties are present. The sleep data do have ties, so check.
   - Re-read any prose that quotes numbers from these outputs.
2. **Mind the global `warning = FALSE`.** It hides the new `brunner_munzel` warnings in the rendered paper. That's fine for presentation, but check the chunk outputs interactively once so no example silently uses a combination the package warns about.
3. **Regenerate** `Avocado_Update.R` (purl), `.tex`, `.pdf` and `.docx` from the updated Rmd.
