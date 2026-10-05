# boot_cor_test() studentized bootstrap simulations (2026-10-05)

Scratch scripts behind the `boot_cor_test()` changes on this branch: real
studentization (influence-function SEs), the HC4 correction for Pearson,
`boot_scale`, and the `boot_ci = "auto"` default.

| File | What it does |
|---|---|
| `stud_check.R`, `stud_check2.R` | First coverage checks: constant-SE "stud" vs ADF / jackknife SEs, z vs r scale (Pearson only, self-contained) |
| `se_check.R` | Accuracy of `.cor_se()`: mean SE vs Monte Carlo SD for each method and DGP (n = 50, 300) |
| `kendall_se30.R` | Kendall SE bias and variability at n = 20, 30 |
| `cov_sim_par.R` | Main coverage study: 5 DGPs x 3 methods x n = 30, 80; 1000 reps, R = 999; 12 parallel workers |
| `cov_sim_results_2026-10-05.txt` | Results of the main study (29 of 30 cells; column definitions in its header) |

Notes:

- All scripts call `devtools::load_all("C:/GitHub/TOSTER")`, so they use the
  package code as it is when you run them.
- In `cov_sim_par.R`, `stud_IF` uses the package's `.cor_se()`. When the saved
  results were produced, that was the plain ADF SE for Pearson. The package now
  applies the HC4 correction inside `.cor_se()`, so a rerun gives
  `stud_IF == stud_HC4` for Pearson.
- `SMOKE=1 NSIM=3` runs a quick serial check of the t3 Pearson cells (output
  goes to a temp file). The full run (~70 min on 12 cores) writes
  `cov_sim_par.out` in this folder.
