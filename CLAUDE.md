# CLAUDE.md - TOSTER Package Development Guide

## Package Overview

**TOSTER** (Two One-Sided Tests) - Version 0.9.0

Implements the two one-sided tests (TOST) procedure to test equivalence for t-tests, correlations, differences between proportions, and meta-analyses. Includes power analysis for t-tests and correlations. Allows specifying equivalence bounds in raw scale units or effect sizes.

Primary citation: Lakens (2017) <doi:10.1177/1948550617697177>

## Development Commands

```r
devtools::load_all()     # Load package for interactive development
devtools::test()         # Run testthat tests
devtools::check()        # Full R CMD check
devtools::document()     # Generate documentation from roxygen2
```

## Code Style

### Naming Conventions
- **Main functions**: `snake_case` (e.g., `t_TOST`, `boot_t_TOST`, `brunner_munzel`, `wilcox_TOST`)
- **Legacy functions**: `camelCase` or mixed (e.g., `TOSTone`, `TOSTpaired`, `powerTOSTtwo`)
- **Internal helpers**: `snake_case` with descriptive names (e.g., `d_est_pair`, `perm_loop`)
- **S3 methods**: `generic.class` pattern (e.g., `print.TOSTt`, `plot.TOSTt`)

### Function Patterns
- Use `match.arg()` for argument validation with predefined choices
- Provide both `default` and `formula` methods for data analysis functions
- Include `...` for extensibility in generic functions
- Use lifecycle badges for stability indicators

### Roxygen2 Documentation
- Markdown enabled (`Roxygen: list(markdown = TRUE)` in DESCRIPTION)
- Use `@family` tags to group related functions:
  - `@family TOST` - main equivalence test functions
  - `@family power` - power analysis functions
  - `@family Robust tests` - non-parametric alternatives
  - `@family htest` - hypothesis test utilities
- Include comprehensive `@details` sections with mathematical notation using `\eqn{}`
- Provide multiple `@examples` demonstrating different use cases

## S3 Methods

### Output Classes

**`TOSTt`** - Primary class for t-test based TOST results
- Components: `TOST` (data.frame), `eqb`, `alpha`, `method`, `hypothesis`, `effsize`, `smd`, `decision`, `data.name`, `call`
- Methods: `print.TOSTt`, `plot.TOSTt`, `describe.TOSTt`

**`TOSTnp`** - Non-parametric TOST results (Wilcoxon-based)
- Methods: `print.TOSTnp`, `describe.TOSTnp`

**`htest`** - Standard R hypothesis test class (used by `brunner_munzel`, correlation tests)
- Use `as_htest()` to convert `TOSTt`/`TOSTnp` to `htest`

**`power.htest`** - Power analysis results (standard R class)
- Returned by `power_t_TOST` and similar power functions
- Contains: `power`, `beta`, `alpha`, `n`, `delta`, `sd`, `bounds`, `note`, `method`

### Method Implementation Pattern
```r
# Generic function
describe <- function(x, ...) {
  UseMethod("describe")
}

# Method implementation
describe.TOSTt <- function(x, digits = 3, ...) {
  # Implementation
}
```

## Package Structure

### R/
- **Core TOST functions**: `t_TOST.R`, `wilcox_TOST.R`, `log_TOST.R`, `brunner_munzel.R`
- **Legacy functions**: `TOSTone.R`, `TOSTtwo.R`, `TOSTpaired.R`, `TOSTr.R`, `TOSTmeta.R`
- **Effect size calculators**: `smd_calc.R`, `ses_calc.R`, `cohend_calcs.R`, `rbs_calcs.R`
- **Power functions**: `power_t_TOST.R`, `powerTOSTone.R`, `powerTOSTtwo.R`, `power_correlations.R`
- **Bootstrap functions**: `boot_t_TOST.R`, `boot_smd_calc.R`, `boot_cor_test.R`
- **S3 methods**: `methods.TOSTt.R`, `methods.TOSTnp.R`
- **Plotting**: `plot.R`, `plot_smd.R`, `plot_corr.R`, `gg_curv_t.R`
- **Utilities**: `htest.R`, `htest_helpers.R`, `others.R`
- **Jamovi interface**: `00jmv.R`, `datatostone.*.R`, `datatostpaired.*.R`

### tests/
- Framework: testthat (edition 3)
- Main test file: `tests/testthat.R`
- Test files follow `test-*.R` naming convention
- Key test patterns:
  - `test-known_results.R` - regression tests against published values
  - `test-tTOST.R` - comprehensive t_TOST function tests
  - `test-power_t_TOST.R` - power calculation validation

### vignettes/
- `IntroTOSTt.Rmd` - Introduction to t-test based TOST
- `IntroductionToTOSTER.Rmd` - Package overview
- `SMD_calcs.Rmd` - Standardized mean difference calculations
- `robustTOST.Rmd` - Robust/non-parametric methods
- `correlations.Rmd` - Correlation equivalence tests
- `the_ftestTOSTER.Rmd` - F-test equivalence (ANOVA)

## Testing

### Running Tests
```r
devtools::test()                    # Run all tests
testthat::test_file("tests/testthat/test-tTOST.R")  # Single file
```

### Test Coverage
```r
covr::package_coverage()            # Check coverage
covr::report()                      # HTML coverage report
```

Target: 100% code coverage (minimum 90%)

### Writing Tests
- Test against known published results with minimal tolerance:
```r
expect_equal(result$p.value, 0.0006040021, tolerance = 1e-5)
```
- Use `hush()` helper to suppress console output during tests
- Test both error conditions (`expect_error()`) and warnings (`expect_warning()`)
- Verify consistency between related functions (e.g., `t_TOST` vs `tsum_TOST`)

## Documentation

### Roxygen2 Patterns
```r
#' @title Function Title
#' @description
#' `r lifecycle::badge('stable')`
#'
#' Description text...
#'
#' @param x Description of x
#' @param ... further arguments passed to methods
#'
#' @details
#' Detailed explanation...
#'
#' @return Description of return value
#'
#' @examples
#' # Example code
#'
#' @family TOST
#' @export
```

### Lifecycle Badges
- `stable` - API is stable
- `maturing` - API may change
- `deprecated` - Use alternative

## Dependencies

### Imports (required)
- **stats**: Core statistical functions
- **graphics**: Base plotting
- **ggplot2**, **ggdist**, **distributional**, **cowplot**: Visualization
- **tidyr**: Data manipulation
- **jmvcore**: Jamovi integration
- **R6**: OOP for Jamovi classes
- **lifecycle**: Deprecation management

### Suggests (optional)
- **testthat** (>= 3.0.0): Testing
- **knitr**, **rmarkdown**: Vignettes
- **broom**: Tidy output
- **car**, **afex**: ANOVA extensions

## References

- `/references/general/` - Persistent package documentation and style references (tracked in git)
- `/references/working/` - Branch-specific methodology papers and notes (gitignored)

When implementing new equivalence testing methods, check `/references/working/` first for source material. New statistical functions should follow established approaches in the equivalence testing literature.

## Equivalence Testing Context

### Core Concepts
- **TOST (Two One-Sided Tests)**: Procedure to demonstrate that an effect falls within specified bounds
- **Equivalence bounds**: The smallest effect size considered practically meaningful
- **Minimal effects testing (MET)**: Testing whether an effect exceeds a meaningful threshold

### Hypothesis Types
- `hypothesis = "EQU"`: Equivalence testing (null: effect outside bounds)
- `hypothesis = "MET"`: Minimal effects testing (null: effect inside bounds)

### Key Parameters
- `eqb`: Equivalence bounds (single value for symmetric, two values for asymmetric)
- `eqbound_type`: "raw" (recommended) or "SMD"
- `alpha`: Significance level (default 0.05)

### Function Requirements
New statistical test functions typically need:
1. Effect size specification (raw and/or standardized)
2. Equivalence bounds handling
3. Confidence intervals (typically 1-2*alpha for TOST)
4. Both `default` and `formula` S3 methods
5. Appropriate return class (`TOSTt`, `htest`, or custom)
