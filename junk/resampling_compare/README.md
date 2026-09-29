# Resampling comparison: perm_t_test vs boot_t_test (stud, BCa)

A simulation study of how `perm_t_test()` (studentized permutation) and
`boot_t_test()` (`boot_ci = "stud"` and `"bca"`) perform across two-sample and paired
designs, with an analytic reference: Welch's t for `tr = 0`, Yuen's trimmed t for
`tr > 0`. The goal is to identify the scenarios in which each method holds its
nominal error rates, for both standard tests and TOST.

## Files

| file | purpose |
|---|---|
| `setup.R` | distributions, true values, reference tests, one replication, cell definitions |
| `validate.R` | pre-flight checks and runtime estimate (run first) |
| `run.R` | runs one block, saves after every cell, resumes if interrupted |
| `summarize.R` | combines finished blocks into `summary_long.csv` and `summary_tables.md` |
| `results/` | created by `run.R` |

## Requirements

- R with `devtools`, `parallel` (base), and optionally `WRS2` (used only to check
  the Yuen implementation in `validate.R`).
- The TOSTER source from the branch with the `perm_t_test` CI fix (issue #120),
  **not** the CRAN release. The scripts stop with an error if an older TOSTER is
  loaded, because its permutation intervals are the old percentile intervals.

## How to run

All commands run from the **package root** (the folder containing `DESCRIPTION`).

```bash
Rscript junk/resampling_compare/validate.R
```

```bash
SIM_BLOCK=A Rscript junk/resampling_compare/run.R
```

```bash
SIM_BLOCK=B Rscript junk/resampling_compare/run.R
```

```bash
SIM_BLOCK=C Rscript junk/resampling_compare/run.R
```

```bash
Rscript junk/resampling_compare/summarize.R
```

- **Parallel runs:** blocks are independent, so they can run at the same time, e.g.
  on different machines or in separate sessions with `SIM_CORES` split between them.
  Each block writes its own file.
- **Interruptions:** re-running the same block resumes where it stopped. Completed
  cells are skipped, and the settings must match the original run.
- **Collecting results:** copy the `results/results_block_*.rds` files into one
  `results/` folder and run `summarize.R`. It also works with a subset of blocks.
- **Windows:** set variables with `set SIM_BLOCK=A` (cmd) or `$env:SIM_BLOCK='A'`
  (PowerShell) before calling `Rscript`.

### Settings (environment variables)

| variable | default | meaning |
|---|---|---|
| `SIM_BLOCK` | (required) | `A`, `B`, or `C` |
| `SIM_NSIM` | 2000 | replications per cell (MC SE near .05 is .0049; 4000 gives .0034) |
| `SIM_R` | 1999 | permutations and bootstrap resamples (same for both methods) |
| `SIM_CORES` | all cores minus one | worker processes |
| `SIM_SD_RATIO` | 2 | SD ratio in the "unequal SD" cells |
| `SIM_PKG` | `.` | TOSTER source folder (used with `devtools::load_all`) |
| `SIM_DIR` | `junk/resampling_compare` | location of these scripts |
| `SIM_OUT_DIR` | `<SIM_DIR>/results` | where results are written |

Results are reproducible for a given `SIM_NSIM`, `SIM_R`, and `SIM_SD_RATIO`
regardless of the number of cores, because each replication has its own seed.

### Expected runtime

Measured with `validate.R` on a 13-core Windows machine, R = 1999, nsim = 2000:

| block | cells | ~time per replication | ~hours on 13 cores |
|---|---|---|---|
| A | 130 | 1.1 s | 6.3 |
| B | 105 | 2.3 s | 10.4 |
| C | 62 | 2.3 s | 6.2 |

Time scales linearly with `SIM_NSIM` and roughly inversely with `SIM_CORES`.
`validate.R` prints the estimate for your machine.

## Design

### Methods (all evaluated on the same data sets)

- `perm_t_test`: package defaults (studentized, Welch SE), `R = SIM_R`. The 95%
  interval uses the default symmetric two-sided rule. The 90% interval uses
  `symmetric = FALSE`, which is the equal-tailed interval `perm_t_test` uses for TOST.
- `boot_t_test`, `boot_ci = "stud"` and `"bca"`, `R = SIM_R`.
- Reference: Welch's t (`tr = 0`), Yuen's trimmed t with Yuen–Welch df (`tr > 0`).
  For paired designs, the one-sample (trimmed) t on the differences.

`validate.R` checks that the 90% intervals are exactly the intervals each function
uses for TOST.

### What is measured

Everything is computed from each method's confidence intervals, which agree with the
method's own p-values by construction. Data are generated with a known true value θ
of the estimand. θ is the mean difference, or the difference in population trimmed
means when `tr > 0`, computed by integrating the quantile function (checked by Monte
Carlo in `validate.R`).

| metric | definition | nominal |
|---|---|---|
| `type1` | 95% CI misses θ (two-sided Type I error of the standard test at θ) | .05 |
| `err_below`, `err_above` | 90% CI entirely below / above θ | .05 each |
| `one_sided` | max(`err_below`, `err_above`) | .05 |
| `width95` | mean 95% CI width | lower is better |
| `power0` | 95% CI excludes 0 (power when θ ≠ 0) | higher is better |
| `tost05`, `tost1` | 90% CI inside ±0.5 / ±1 (TOST power when θ = 0) | higher is better |
| `fail` | method errored or returned a non-finite interval | 0 |

**Why `one_sided` answers the TOST question:** if θ were an equivalence bound, the
TOST Type I error at that bound is the probability that the 90% CI lies entirely on
the wrong side of it.
- At an upper bound: TOST error = `err_below`, minimal-effect error = `err_above`.
- At a lower bound: the roles swap.

So `one_sided` is the worst-case TOST (and minimal-effect) Type I error.
- **Exactness:** this holds exactly for the bootstrap and the analytic tests, whose
  intervals shift with the data. For `perm_t_test` it is a close approximation. The
  permutation reference distribution is built from unshifted data, and a separate
  study (`junk/perm_t_test_shifted_null.md`) found this makes no practical
  difference under equal shapes.

### Scenarios

All continuous distributions are standardized to mean 0 and SD 1. `x` is always the
distinctive group: the more variable one in "unequal SD" and "unequal spread" cells,
and the skewed one in "different shape" cells.

**Distributions:**
- continuous: normal; lognormal (σ = 0.5, moderate skew); exponential (strong skew);
  t₃ (heavy tails); contaminated normal (10% from N(0, 5²), outliers)
- discrete: 5-point Likert (mid, wide, skewed); paired differences on −2..2
  (symmetric, skewed)

**Two-sample layouts:**
- balanced n = 8, 15, 30, 50
- 1:2 with n = 8, 15, 30 in the smaller group
- for non-identical structures, both "x smaller (1:2)" and "x larger (2:1)"

| block | contents |
|---|---|
| A (`tr = 0`) | **Two-sample:** identical distributions (sharp null) for all continuous distributions and Likert-mid; unequal SD (ratio `SIM_SD_RATIO`) for all continuous distributions; Likert unequal spread (wide vs mid). **Paired:** all continuous distributions plus discrete symmetric and skewed differences, n = 8, 15, 30, 50. |
| B (`tr = 0.2`) | As block A for the continuous distributions (Likert and discrete cells excluded). |
| C | Different shapes: exponential vs normal at `tr` = 0 and 0.2, and skewed vs mid Likert at `tr` = 0 (the weak null for means). Power: true shift 0.5 for normal, exponential, t₃ and contaminated data, n = 15 and 30, two-sample and paired, `tr` = 0 and 0.2. |

TOST power (θ = 0, bounds ±0.5 and ±1) comes from the identical-distribution and
symmetric paired cells in blocks A and B at no extra cost.

## Output

`summarize.R` writes to `results/`:

- `summary_long.csv`: one row per cell × method with all metrics. Use this for
  plots or custom tables.
- `summary_tables.md`: markdown tables by block and design.
  - two-sided Type I error, and one-sided (TOST) error
  - rates more than 2 MC SE above .05 in **bold** (liberal), more than 2 MC SE
    below in _italics_ (conservative)
  - power tables, TOST power tables, and a table of method failures

## Notes

- **Unequal SD with trimmed skewed data:** with `tr > 0` and skewed data, the true
  value in "unequal SD" cells is not 0, because scaling a skewed distribution changes
  its trimmed mean. The true value is computed exactly, and all metrics are relative
  to it.
- **BCa failures:** BCa can fail (undefined acceleration) with heavily tied or
  degenerate small samples, e.g. Likert data at n = 8. These cases are counted in
  `fail` rather than silently dropped, and the other metrics are averaged over
  successful data sets. `summary_long.csv` includes `n_valid`, the smallest number
  of successful data sets for any method in that cell.
