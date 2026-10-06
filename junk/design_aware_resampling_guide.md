# Design-aware `perm_t_test()` / `boot_t_test()`: implementation guide

*Planning document. Not for the package. Written 2026-09-29, after the issue #120
(CI/p-value alignment) work in `perm_t_test()`.*

## 0. Goal and scope

Let `perm_t_test()` and `boot_t_test()` respect common randomization/sampling designs
without adding mixed-model machinery:

| Design | New argument(s) | Section |
|---|---|---|
| Blocked / stratified randomization | `strata` | 3 |
| Cluster randomization (assignment at cluster level) | `cluster` | 4 |
| Pre-post / one-sample data with clustering (students in schools) | `cluster` + `paired = TRUE` (or one-sample) | 5 |
| 2x2 AB/BA crossover | `crossover_diff()` helper, then the *existing* two-sample tests | 6 |
| Strata + cluster (clusters nested in blocks) | `strata` + `cluster` | 3.4 |

**Out of scope for now:** `tr > 0` combined with `cluster`/`strata`, `boot_t_TOST()`
(inherits later), designs with more than 2 periods/treatments, multi-level nesting deeper
than one cluster level (resample/permute at the *top* level and say so in the docs), and
user-supplied observation weights (see 8).

### Decisions already made

- **Estimand weighting:** individual-weighted by default. New argument
  `weights = NULL` (individual-weighted) or `weights = "cluster"` (each cluster counts once).
  Nothing else is accepted in v1. The argument is meaningful only when `cluster` is supplied.
- **Sequencing:** finish issue #120 first (done on this branch); implement this in a
  separate branch so the two changes stay reviewable.

## 1. What the current code does (what constrains the design)

Read from `R/perm_t_test.R` and `R/boot_t_test.R`:

- **`perm_t_test.default`** works on plain vectors `x`, `y`. `paired = TRUE` is converted to
  a one-sample problem on `x - y` early, so everything downstream is "one-sample" or
  "two-sample".
  - Permutation generators return *index/sign matrices*: `perm_signs(n, R)` (sign flips of
    centered data) and `perm_groups(nx, ny, R)` (`combn` for exact, unique random draws
    otherwise). Both dedupe rows; the `"plusone"` p-value's exactness argument (Phipson &
    Smyth) relies on sampling *without* replacement, so any new generator must too.
  - The statistic is a **loop** over permutations, recomputing means and (if `perm_se`) the SE.
  - CI = est -/+ stderr * `perm_crit(TSTAT, ...)`, i.e. inversion of the studentized test
    using the *unshifted* permutation distribution (see `junk/perm_t_test_shifted_null.md`).
    **Any new engine must produce a `TSTAT` vector that is independent of `mu`, or the
    #120 alignment fix is lost.**
- **`boot_t_test.default`** handles `tr == 0` with vectorized `matrix(sample(...))`, and
  `tr > 0` with a loop. It bootstraps x and y **separately** (already stratified by arm),
  uses centered data for `TSTAT` (null distribution) and uncentered data for `m_vec` /
  `m_se_vec` (CI). BCa uses a delete-one jackknife over all observations. For `tr == 0`
  the initial estimate comes from `simple_htest()` -> `t.test()`, which cannot know about
  clusters.
- **Formula methods** split the response by a 2-level factor into `x`/`y` and call the
  default method. Extra columns (cluster, strata) would be lost at that step.
- `boot_pvalue()` (`R/corr_calcs.R`), `stud()`, `basic()`, `perc()`, `bca_ci()`,
  `bca_params()` (`R/others.R`) are reusable as-is; they only see `m_vec`, `m_se_vec`,
  `TSTAT`.
- `sample_size` and `estimate` names are consumed elsewhere (6 files reference
  `sample_size`; formula methods rewrite names). Grep for consumers before adding entries.

**Consequence:** don't rewrite the existing paths. Add a *separate* design path
(new files) that is entered only when `strata` or `cluster` is non-NULL. The current
behavior stays byte-identical, and the new path is validated against it (section 9).

## 2. Architecture

```
R/design_helpers.R    # validation, design object, summaries, samplers, stat kernels
R/perm_t_test_design.R  # perm_t_test_design(): internal, called from perm_t_test.default
R/boot_t_test_design.R  # boot_t_test_design(): internal, called from boot_t_test.default
R/crossover_diff.R    # exported helper + print method
```

In both `.default` methods, after argument validation and before any data handling:

```r
if (!is.null(strata) || !is.null(cluster)) {
  return(perm_t_test_design(x = x, y = y, strata = strata, cluster = cluster,
                            weights = weights, paired = paired, alternative = alternative,
                            mu = mu, var.equal = var.equal, alpha = alpha, tr = tr, R = R,
                            symmetric = symmetric, perm_se = perm_se, p_method = p_method,
                            keep_perm = keep_perm, dname = dname, call = match.call()))
}
```

### 2.1 New arguments

Add to both `.default` signatures (after existing args, before `...`):

```r
strata  = NULL,   # vector aligned with the data (see 2.2)
cluster = NULL,   # vector aligned with the data (see 2.2)
weights = NULL    # NULL = individual-weighted; "cluster" = equal weight per cluster
```

Validation (`design_check_args()`):

- `weights` must be `NULL` or `"cluster"`. If `weights` is non-NULL and `cluster` is NULL: error.
- `tr > 0` with `cluster`/`strata`: error in v1 (except `weights = "cluster"` could allow
  trimming of cluster means via the existing code path; defer, but leave a hook in the error
  message).
- `var.equal = TRUE` with `cluster`: message that it is ignored (cluster-robust variance is
  Welch-type); with `strata` only: same.
- NA handling: **drop rows, not values**. The current code does `x[!is.na(x)]`, which
  misaligns any companion vector. In the design path, build one data frame
  (`value`, `group`, `cluster`, `strata`, `pair_id`), and `complete.cases()` on it.

### 2.2 Aligning `cluster`/`strata` with `x`, `y`

- **One-sample / paired:** vectors of length `length(x)` (number of subjects or pairs).
- **Two-sample:** vectors of length `length(x) + length(y)`, in `c(x, y)` order. Document
  this clearly. The formula interface avoids the awkwardness (below).
- **Formula method:** add `strata` and `cluster` as arguments evaluated in `data`, the same
  way `subset` is. Implement by *not* using the split-by-factor shortcut when either is
  non-NULL:

  ```r
  perm_t_test.formula <- function(formula, data, subset, na.action,
                                  strata = NULL, cluster = NULL, ...)
  ```

  Put `strata`/`cluster` into the `model.frame` call (capture them with `substitute`, add
  to `m` like `subset`; `model.frame` then handles `subset` and `na.action` consistently
  across all columns). Then build `x`, `y`, and the aligned `cluster`/`strata` vectors as
  `c(cluster[g == lvl1], cluster[g == lvl2])`. Keep the existing group-label renaming
  logic and extend the `sample_size` renaming (below).

  Paired + formula stays discouraged as it is today. For paired + cluster the recommended
  entry point is the default method with `x`, `y` vectors, or wide-format data.

### 2.3 The design object

`design_prep()` returns a list, the single source of truth for both engines:

```r
list(
  type      = "two_sample" | "one_sample",   # paired -> one_sample on d
  value     = numeric,          # observation-level y (or d = x - y)
  group     = factor(2 levels) or NULL,
  cluster   = integer index or NULL (NULL => every row is its own cluster),
  strata    = integer index or NULL (NULL => one stratum),
  n_obs, n_clusters (per group, per stratum),
  weights   = NULL | "cluster",
  assign_level = "cluster" | "unit"   # what gets permuted/resampled
)
```

Checks inside `design_prep()`:

1. **Treatment constant within cluster** (two-sample with cluster): each cluster must belong
   to a single group. If not, error: cluster-level permutation would be the wrong scheme;
   tell the user this is not a cluster-randomized design (within-cluster assignment would
   need a paired/stratified analysis with `strata = cluster` instead).
2. **Clusters nested in strata:** each cluster in exactly one stratum.
3. **Minimum sizes:** at least 2 clusters per arm (per stratum, if strata; see 3.3 for the
   stratum rule). Warn below ~10 clusters per arm (p-values coarse or bootstrap unreliable).
4. **Singleton reduction:** `cluster = NULL` is treated as one cluster per row, which makes
   the individual-level code a special case of the cluster code (used for testing, see 9).

## 3. Blocked / stratified randomization (`strata`)

### 3.1 Estimator

For strata b = 1..B, arms x, y, with stratum arm sizes `n_xb`, `n_yb`, `n_b = n_xb + n_yb`,
`N = sum n_b`:

- estimate: `d = sum_b W_b (xbar_b - ybar_b)`
- default weights `W_b = n_b / N` (individual-weighted, matches the decision above);
  `weights = "cluster"` does **not** apply to strata-only designs (error: use the default).
- variance (Welch-type): `V = sum_b W_b^2 (s2_xb / n_xb + s2_yb / n_yb)`
- df: Satterthwaite over the 2B variance components.
- stat: `t = (d - mu) / sqrt(V)`.

### 3.2 Permutation and bootstrap schemes

- **Permutation:** permute group labels *within* each stratum, preserving `n_xb`. Number of
  distinct assignments = `prod_b choose(n_b, n_xb)`.
  - Exact enumeration: `expand.grid` (or iterative Kronecker construction) over per-stratum
    `combn` results. **Add a guard** (e.g. `max_exact = 1e6`): if the count exceeds it and
    `R = NULL`, error asking for a finite `R`. The current unstratified code would silently
    try to enumerate.
  - Sampling: per stratum, draw `R` random subsets of size `n_xb` with
    `t(apply(matrix(runif(R * n_b), R), 1, rank)) <= n_xb` (fast enough for typical
    strata counts), `cbind` the strata, then dedupe rows with the same retry loop as
    `perm_groups()`.
  - Output: an `R_used x n_units` logical matrix `in_x` plus the same `exact`,
    `R.used`, `max_perms` fields `perm_groups()` returns, so the downstream p-value and
    CI code (`compute_perm_pval`, `perm_count`, `perm_crit`) is unchanged.
- **Bootstrap:** resample observations within (stratum x arm). Resample counts equal the
  original cell counts, so arm sizes per stratum are preserved.

### 3.3 Variance in small strata

Welch variances need `n_xb >= 2` and `n_yb >= 2`. **v1: error** if any stratum has fewer,
with a message pointing to (a) merging strata or (b) a paired analysis for pairs (1 vs 1).
Possible v2: a pooled-residual variance
`V = sum_b W_b^2 s2_pooled (1/n_xb + 1/n_yb)` with `s2_pooled` across strata.

### 3.4 Strata + cluster

Clusters nested in strata, assignment at cluster level within strata. Same estimator as 3.1
but with the cluster-ratio means and variances of section 4 computed per stratum. The
sampler operates on *clusters* within strata (3.2 applied to the cluster table).

## 4. Cluster randomization (`cluster`)

### 4.1 Cluster summaries (the key simplification)

Collapse to one row per cluster: `k` clusters, with size `n_i`, sum `S_i`, group `g_i`,
stratum `b_i`. Define per-cluster weight and weighted total:

| `weights` | `w_i` | `T_i` |
|---|---|---|
| `NULL` (individual-weighted) | `n_i` | `S_i` (cluster sum) |
| `"cluster"` | 1 | `S_i / n_i` (cluster mean) |

Arm mean (ratio estimator): `m_g = sum_{i in g} T_i / sum_{i in g} w_i`.

Cluster-robust variance of `m_g` (ratio-estimator / linearization):

```
V_g = [k_g / (k_g - 1)] * sum_{i in g} (T_i - m_g * w_i)^2 / (sum_{i in g} w_i)^2
```

Welch-type combination: `se = sqrt(V_x + V_y)`,
`df = se^4 / (V_x^2 / (k_x - 1) + V_y^2 / (k_y - 1))`.

**Sanity property (used as a test):** with `w_i = 1` and `T_i = x_i` (every observation its
own cluster, or `weights = "cluster"` with singleton clusters), `V_g = s_g^2 / k_g`, so the
result reduces exactly to the existing Welch t-test. With `weights = "cluster"` the whole
thing equals the existing code applied to cluster means.

The estimand and label should say so, e.g. "difference in means (individual-weighted)".

### 4.2 Permutation engine (vectorized, no per-permutation loop)

The statistic depends only on `(w_i, T_i)` and group membership, so the whole permutation
distribution comes from matrix products. With `M` = `R_used x k` 0/1 matrix (row = cluster
assignment to arm x; the complement is arm y):

```r
Wx  <- M %*% w                     # sum w in arm x
Tx  <- M %*% T
mx  <- Tx / Wx
SSx <- M %*% T^2 - 2 * mx * (M %*% (T * w)) + mx^2 * (M %*% w^2)   # sum (T - m w)^2
Vx  <- (kx / (kx - 1)) * SSx / Wx^2          # kx is constant across perms
# arm y: same with (1 - M)
TSTAT <- ((mx - my)) / sqrt(Vx + Vy)        # perm_se = TRUE
```

- `perm_se = FALSE`: divide by the observed `se` (as the current code does).
- Sampler: `perm_groups(kx, ky, R)` applied to the cluster table (no strata) or the
  strata-wise sampler of 3.2 (strata). Reuse `perm_groups` for the unstratified cluster
  case unchanged: it already returns `idx_x` rows of length `kx`; convert to `M`.
- The permuted statistic is centered implicitly: the current code compares
  `t_obs = (d - mu)/se` with the permutation distribution of the *uncentered permuted*
  `(mx* - my*)/se*`, and `EFF <- (mx_perm - my_perm) + diff_means`. Mirror this exactly:
  `TSTAT[i] = (mx* - my*) / se*`, `EFF[i] = (mx* - my*) + d`. Do not center the cluster
  totals separately; the permutation distribution is symmetric about the pooled value.
- **Design-consistency of `Wx`:** note `w_i` under individual weighting makes the arm
  sample sizes vary across permutations (`Wx` changes). That is correct for this
  estimand. It also means permuting *cluster labels* with unequal cluster sizes gives
  a wider reference distribution than permuting individuals: which is the whole point.
- Few clusters: `choose(k, kx)` small means coarse p-values. Emit the existing
  "insufficient R for alpha" message logic; add: if `max_perms < 1/alpha`, warn that the
  test cannot reject at the chosen alpha (the CI code already returns `+/-Inf` in that case).

### 4.3 Bootstrap engine

Resample **clusters** within arm (within stratum x arm if strata). With `idx` an `R x k_g`
matrix of resampled cluster indices for arm g:

```r
Wg   <- rowSums(matrix(w[idx], nrow = R))
Tg   <- rowSums(matrix(T[idx], nrow = R))
mg   <- Tg / Wg
SSg  <- rowSums((matrix(T[idx], nrow = R) - mg * matrix(w[idx], nrow = R))^2)
Vg   <- (kg / (kg - 1)) * SSg / Wg^2
```

Mirror the existing structure:

- **CI vectors (`m_vec`, `m_se_vec`):** resample the *uncentered* cluster table:
  `m_vec = mx* - my*`, `m_se_vec = sqrt(Vx* + Vy*)`.
- **Null distribution (`TSTAT`):** the existing code centers x and y to a common mean
  `mz` (`x - mx + mz`) before resampling. For clusters, centering shifts every
  observation, so `T_i^c = T_i - w_i * (m_g - m_z)` for individual weighting (each
  observation shifted; the sum shifts by `n_i * shift`); for `weights = "cluster"`
  (`w_i = 1`, `T_i` = mean) the shift is subtracted directly. `m_z` = pooled ratio mean
  over both arms. Then `TSTAT = (mx*_c - my*_c) / se*_c`. (The current code resamples
  centered data for TSTAT and uncentered for the CI: keep both sets.)
- **BCa:** jackknife = delete-one-*cluster* over all clusters (both arms), computing
  `m_x - m_y` each time. This parallels the current pooled delete-one-observation
  jackknife. Do not delete individual observations.
- Reuse `stud()`, `basic()`, `perc()`, `bca_ci()`, `bca_params()`, `boot_pvalue()` unchanged.
- Return `boot = m_vec`, `stderr = sd(m_vec)`, as now.

## 5. Pre-post (paired) or one-sample with clusters

Convert to differences `d_i = x_i - y_i` per pair (as the code already does), then treat as
one-sample data with `cluster` (students in schools). Cluster summaries as in 4.1, but only
one arm: `m = sum T_i / sum w_i`, `V = [k/(k-1)] sum (T_i - m w_i)^2 / (sum w_i)^2`,
`df = k - 1`, `t = (m - mu) / sqrt(V)`.

- **Bootstrap:** resample clusters (with all their pairs); `TSTAT` from centered totals
  `T_i^c = T_i - w_i * m`. Same as 4.3 with one arm.
- **Permutation:** sign flips. Two levels are conceivable; choose one default:
  - `cluster` (recommended default): flip the sign of *all* `d_i` in a cluster together;
    i.e. flip `T_i` (and keep `w_i`). `2^k` sign patterns; reuse `perm_signs(k, R)`. This
    respects within-cluster dependence of the changes and matches the
    "unit of inference = cluster" logic.
  - unit-level flips: valid only if the pre/post labels are exchangeable within each
    student *and* the statistic's variance accounts for clustering. Not offered in v1.
  Vectorized: `S %*% ...` with `S` the `R x k` sign matrix: `m* = (S %*% (T^c)) / sum(w)`,
  `V* = [k/(k-1)] * (S^2-free)`: note `(T_i^c s_i - m* w_i)^2` requires `m*`, so compute
  `SS* = (S^2) %*% (T^c)^2 - 2 m* (S %*% (T^c w)) + m*^2 sum(w^2)` where `S^2 = 1`.
- **Caveat to document:** paired sign-flipping tests symmetry about zero of the (cluster
  total) differences, not only equality of means. Skewed differences make all methods
  over-reject (see `junk/perm_t_test_shifted_null.md`); clustering with few schools makes
  this worse. Advise `boot_t_test()` and a look at the differences when `k < ~15`.
- `weights = "cluster"` reduces to the existing paired test on school-mean differences.

## 6. `crossover_diff()`: 2x2 AB/BA crossover helper

### 6.1 The math (why a helper, and why not a new engine)

Model for subject `i` in period `j` on treatment `k`:
`Y_ij = mu + pi_j + tau_k + s_i + e_ij`. Let `P_i = Y_i1 - Y_i2` (period-1 minus period-2).

- Sequence AB: `E[P] = (pi_1 - pi_2) + (tau_A - tau_B)`
- Sequence BA: `E[P] = (pi_1 - pi_2) - (tau_A - tau_B)`

So `E[P_AB] - E[P_BA] = 2 (tau_A - tau_B)` and the period effect cancels in the
*between-sequence* difference. Define `u_i = P_i / 2`. Then

> **`mean(u | AB) - mean(u | BA) = tau_A - tau_B`**

which is exactly the mean-difference estimand of an ordinary two-sample t-test with
sequence as the grouping factor. The subject random effect `s_i` drops out because it is
common to both periods.

Consequences:

- **No new engine.** `perm_t_test(diff ~ sequence, data = out, ...)` and
  `boot_t_test(...)` already do the right thing, including TOST and MET via
  `alternative`. Bounds `mu` are in original outcome units for `tau_A - tau_B`.
- **Why the equal-weighted between-sequence difference matters.** The naive alternative,
  pooling `d_i = (Y_A - Y_B)` over all subjects, equals
  `[n_AB (tau + delta) + n_BA (tau - delta)] / N` where `delta = pi_1 - pi_2`, so it is
  **biased whenever `n_AB != n_BA`**. The two-sample formulation has no such problem.
  This is a reason to prefer the helper over "let users compute A - B themselves".
- **Randomization justification:** subjects were randomized to sequences. Under
  `H0: tau_A = tau_B` with a period effect, `P_i` is exchangeable across sequences,
  so permuting sequence labels across subjects is a valid randomization test (up to the
  studentization, as elsewhere). The bootstrap resamples subjects within sequence, which
  is what the two-sample bootstrap already does.
- **Site/center blocking:** if randomization was stratified by site, pass
  `strata = site` (section 3) to the tests. This is the main reason to build the
  strata feature first.
- **Assumptions to state in the docs:** no differential carryover, no
  treatment-by-period interaction, both periods observed. The carryover test (comparing
  subject totals `Y_i1 + Y_i2` between sequences) is low-powered and using it as a gate
  distorts error rates (Freeman, 1989, *Statistics in Medicine*). Provide the totals as a
  diagnostic column but do not automate a two-stage procedure.

### 6.2 API

```r
crossover_diff(data, outcome, subject, period, treatment, sequence = NULL,
               log = FALSE, ref = NULL, na.action = c("drop", "error"))
```

- `data`: long format, one row per subject-period.
- `outcome`, `subject`, `period`, `treatment`, `sequence`: column names (character) or bare
  names via tidy-eval-free `substitute()` + `eval(..., data)`; keep it base-R (the package
  avoids dplyr). Character strings is simpler and consistent with the rest of the package.
- `sequence = NULL`: derive from the data (`treatment` in period-order per subject, e.g.
  `"AB"`); if supplied, verify it is consistent with `treatment`/`period`.
- `log = TRUE`: log-transform the outcome first (bioequivalence convention); contrast becomes
  log(A) - log(B); tell users to set bounds as `log(c(0.8, 1.25))` (works because
  `perm_t_test` accepts asymmetric two-value `mu` for equivalence).
- `ref`: which treatment is the reference "B" (contrast is `test - ref`). Default: the
  second alphabetical/factor level of `treatment`; state it in the printed output.

**Returns** a data.frame with class `c("crossover_diff", "data.frame")`, one row per
subject:

| column | meaning |
|---|---|
| `subject` | id |
| `sequence` | factor, levels ordered `c("<test>-<ref>", "<ref>-<test>")` so that level 1 = test-first sequence (becomes `x` in the formula method) |
| `diff` | `(Y_period1 - Y_period2) / 2` (**the value to analyze**) |
| `total` | `Y_period1 + Y_period2` (carryover diagnostic only) |
| `strata` | (optional) carried through if the input has a `site`-like column, see below |

Attributes: `contrast = "test - ref"` labels, `log`, `n_dropped`, `treatments`.

Sign convention check: the formula method takes `x` = first level, `y` = second, and
`estimate = mean(x) - mean(y)`. With level 1 = test-first sequence, `mean(diff|test-first) -
mean(diff|ref-first) = tau_test - tau_ref`. Verify this in a unit test with simulated data
where the truth is known.

### 6.3 Implementation steps

1. **Validate inputs:** columns exist; exactly 2 unique periods, 2 unique treatments, 2
   unique sequences (else `stop()` with a message that only 2x2 designs are supported and
   pointing to mixed models for anything larger).
2. **Reshape** without dependencies: `split()` by subject, or `reshape()`; sort by
   `period`. For each subject require exactly one row per period. Subjects with a missing
   period (row missing or `NA` outcome): dropped with a message
   (`na.action = "drop"`) and counted in `n_dropped`; error if `"error"`.
3. **Derive sequence** = `paste(treatment_in_period1, treatment_in_period2, sep = "-")`;
   if a `sequence` column was supplied, check agreement. Require both sequences to occur
   and each to have at least 2 subjects.
4. **Check design assumptions cheaply:** each subject received both treatments (else error).
5. **Compute** `diff = (y1 - y2) / 2` (after `log()` if requested; error on non-positive
   outcomes) and `total = y1 + y2`.
6. **Order factor levels** so the test-first sequence is level 1; record `ref`.
7. **Return** the classed data.frame, plus:
   - `print.crossover_diff()`: header (design, n per sequence, dropped subjects, contrast,
     log status) then the head of the data. Short and informative.
   - Optionally `summary`/`carryover_check()`: not in v1. Mention in docs how to run
     `t.test(total ~ sequence)` as a diagnostic, with the Freeman caveat.
8. **Passing strata:** allow an optional `strata = "site"` column argument that is carried
   into the result (subject-constant; validate constancy). The user then calls
   `perm_t_test(diff ~ sequence, data = out, strata = strata, ...)`.

### 6.4 Example (for the roxygen and vignette)

```r
cx <- crossover_diff(dat, outcome = "auc", subject = "id", period = "period",
                     treatment = "trt", log = TRUE)
perm_t_test(diff ~ sequence, data = cx, alternative = "equivalence",
            mu = log(c(0.8, 1.25)), R = 9999)
boot_t_test(diff ~ sequence, data = cx, alternative = "equivalence",
            mu = log(c(0.8, 1.25)))
```

## 7. Return values and printing

- Keep class `"htest"`. Extend `sample_size` only in the design path:
  - individual: `c(nx = , ny = )` as now (observations),
  - add `clusters_x`, `clusters_y` (or `clusters` for one-sample) when `cluster` supplied,
  - the formula methods currently rename `sample_size` only when its length is 2: update
    that rule (or rename by name) so the extended vector is not mangled. Check the other
    consumers (`htest_helpers.R`, `simple_htest.R`, `compare_smd.R`, `hodges_lehmann.R`)
    before adding names.
- `method` strings: "Exact Permutation Cluster-Robust Welch Two Sample t-test",
  "Randomization Permutation Stratified Welch Two Sample t-test", and the bootstrap
  equivalents; mention weighting ("cluster-weighted") when `weights = "cluster"`.
- `estimate`: arm means under the chosen weighting; label the weighting explicitly.
- `null.value` names unchanged.
- Add `design` element to the returned list (strata/cluster counts, weighting, permutation
  unit) so `describe()`-style methods can report it later. Not required by `as_htest()`.

## 8. Explicit non-goals / future

- `weights` as a numeric vector of observation weights. The ratio-estimator machinery
  (`w_i`, `T_i`) already generalizes to it; defer until there is a use case.
- `tr > 0` with clusters (Yuen-type) is unresolved conceptually (trim what: individuals
  or cluster means?). Error for now.
- `boot_t_TOST()` duplicates a lot of `boot_t_test()`. Once the design engine exists it
  should be a thin wrapper; consider that refactor after this work.
- Generalizing to `wilcox_TOST()`, `brunner_munzel()`: the sampler layer (`in_x`/sign
  matrices) is reusable; the statistic kernels are not.

## 9. Testing plan

New file `tests/testthat/test-design_perm_boot.R` (use `hush()`; seeds fixed):

1. **Reduction tests (most valuable):**
   - `cluster = seq_along(x)` (singletons), `weights = NULL` reproduces the existing
     `perm_t_test` (exact permutation, small n, `R = NULL`): same `estimate`, `stderr`,
     `parameter` (df), `p.value`, and `conf.int` to 1e-8. Note the `M` matrix here permutes
     the same labels as `perm_groups`.
   - `weights = "cluster"` equals the existing function applied to hand-computed cluster
     means (exact permutation).
   - `strata` with a single stratum equals the unstratified call.
   - Paired + singleton cluster equals the existing paired test.
2. **Hand-computed values:** a tiny dataset (e.g. 3 clusters per arm, unequal sizes) where
   `m_g`, `V_g`, `df` are computed in the test by explicit formulas; check the cluster
   estimator and SE. Test against published/known values where available (e.g. compare the
   ratio-estimator SE with a `sandwich`-style computation done by hand; don't add a
   dependency).
3. **Crossover:** simulate AB/BA with known `tau`, unequal sequence sizes and a period
   effect; verify `perm_t_test(diff ~ sequence)` estimates `tau` and that the naive pooled
   `mean(d_i)` does not. Round trip: `log = TRUE` ratio, sequence derivation, dropped
   subjects, error on non-2x2 input.
4. **Input validation:** treatment varying within a cluster (error), cluster spanning strata
   (error), stratum with < 2 per arm (error), `weights` without `cluster` (error),
   `tr > 0` with cluster (error), misaligned lengths (error), NA row-dropping keeps
   alignment (test that a missing `x` value drops its cluster/strata entry too).
5. **CI/p-value alignment (issue #120 regression):** for `perm_t_test` in the design path,
   assert `p.value <= alpha` iff the null value lies outside the CI (as in the existing
   #120 tests), across `alternative` values including equivalence/minimal.effect.
6. **Simulation study** (`junk/sim_design_aware_type1.R`, in the style of the existing
   `junk/sim_perm_t_test_*` scripts): type I error (and CI coverage) at the null and at
   the equivalence bounds for
   - cluster-randomized data with ICC in {0, 0.05, 0.2}, cluster counts {6, 10, 20},
     equal and unequal sizes, **comparing the naive `cluster = NULL` analysis against the
     design-aware one**. This shows the anti-conservatism the feature fixes and supports
     the vignette section.
   - blocked randomization with block effects, unequal allocation across blocks;
   - clustered pre-post with school effects, normal and skewed differences;
   - AB/BA crossover with unequal sequences and a period effect.

   Expected shape of results: naive analysis inflated (TOST: false equivalence claims)
   when ICC > 0; design-aware near nominal with >= ~10 clusters/arm; coarse/conservative
   permutation p-values with < 8 clusters/arm.
7. Coverage target stays >= 90% (project goal 100%); the new files are small, so cover
   every validation branch.

## 10. Documentation

- Roxygen for `perm_t_test`/`boot_t_test`: document `strata`, `cluster`, `weights`, the
  alignment rule for two-sample data, the estimand statement, the variance estimator
  (equations in `\eqn{}`), and the few-clusters caveat. `@family Robust tests`/`TOST`
  unchanged.
- New exported `crossover_diff()` with `@family htest` (or a new `@family design`),
  examples, and the carryover caveat.
- Vignette: new section (or a short new vignette, e.g. `design_aware_tests.Rmd`) with:
  a cluster-randomized example, a blocked example, a pre-post-in-schools example, an
  AB/BA bioequivalence example, and the simulation figure.
- `NEWS.md` entry; bump nothing until release.

## 11. Suggested build order

1. `design_prep()` + validation + NA/alignment handling (with tests).
2. Cluster summary kernels (4.1) and the reduction tests. These are the bedrock.
3. Permutation: unstratified cluster (`perm_groups` on clusters), then paired/one-sample
   cluster sign flips, then the stratified sampler and strata (+cluster).
4. Bootstrap: cluster resampling (two-sample, one-sample/paired), strata, BCa jackknife.
5. Formula-method plumbing for `strata`/`cluster`, label and `sample_size` handling.
6. `crossover_diff()` and its tests.
7. Simulation study, then vignette/docs.

Each step is independently testable, and steps 1-2 can ship without touching the public API.

## 12. Open questions for review

1. **Few-cluster policy:** warn only, or refuse below some minimum (e.g. < 4 clusters/arm)?
   Suggest error for < 2, warn for < 10.
2. **Strata-only estimator weights:** `n_b / N` by default (as proposed); should equal
   stratum weights be reachable? `weights = "cluster"` is not applicable there and errors
   in the plan above.
3. **Paired + cluster permutation unit:** confirm cluster-level sign flips as the only v1
   option.
4. **`var.equal` with cluster/strata:** silently ignored with a message, or error?
5. **Argument names:** `strata`/`cluster` vs. `block`/`cluster`. `strata` mirrors
   survey/stratified terminology (and `survival::strata`); `block` is more natural to
   trialists. Pick before the docs are written.
