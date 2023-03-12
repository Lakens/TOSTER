NEWS
================

**TOSTER R package and jamovi module**

# TOSTER v0.7.1

- Fixing the ggplot2 error related to `after_stat` update from that package.

# TOSTER v0.7.0

- Add correlation functions `z_cor_test`, `boot_cor_test`, `corsum_test`, `boot_compare_cor`, and `simple_htest`
- Add `describe` method to provide verbose output for analyses
  - Alternative `describe_htest` function for htest objects

# TOSTER v0.6.0

- Changed Glass's delta SE for paired samples (minor).
- Added `smd_calc` and `ses_calc` for just calculating the standardized effect sizes (no tests).
- Default CIs for for SMDs are now NCT rather than the Goulet method.
- `compare_smd` can be supplied with user provided standard errors.
- Add `log_TOST` and `boot_log_TOST` function for comparing ratios of means.
- Reduce the amount of text in the `print` methods

# TOSTER v0.5.1

- Formatting requirements for jamovi (all superficial changes)

# TOSTER v0.5.0

- Added "compare" functions. 
  - `compare_smd`: Compare 2 SMDs from summary statistics
  - `boot_compare_smd`: Compare 2 SMDs from raw data
  - `compare_cor`: Compare 2 independent correlations
- Added additional SMD options
  - Confidence intervals can now be estimated using other methods
  - smd_ci can be used to set the confidence interval method
  - Glass's delta can now be calculated using the `glass` argument
- Added additional standardized effect sizes for `wilcox_TOST`
  - `ses` argument can be set to "r", "odds", or "cstat" 
  - Respectively, these will provide the rank-biserial correlation, odds, or concordance probability
  
# TOSTER v0.4.2

- Fix to Cohen's drm calculation


# TOSTER v0.4.1

- jmv removed as dependency


# TOSTER v0.4.0

## “Avocado TOST” (Release date:  2022-02-05)

Changes:

-   *New* t\_TOST function
    -   Generalized function for TOST for *any* type of t-test
    -   S3 generic methods to print and plot results
-   *New* power\_t\_TOST function
    -   Generalized function for t_TOST power analysis
    -   Outputs the an object of `power.htest` class
-   Updated all jamovi functions to allow minimal effect tests
    -   Direction of one-sided tests now allows
-   Added equ\_anova and equ\_ftest
    -   Now allows equivalence (also called non-inferiority)
-   jamovi functions using t-test have more plotting options
-   Error in powerTOSTtwo fixed when determining N
-   All old t-test based TOST functions now have message telling users they are defunct
-   All TOST procedures have text results changed to provide more appropriate feedback.

# TOSTER v0.3.4

## (Release date: 2018-08-05)

Changes:

-   Added a verbose = FALSE option to all functions.
-   Added Fisher’s z transformed CI to output of TOSTr.
-   Cleaned up the text and numeric output of the functions.
-   Renamed some variable names in the TOSTtwo.prop function for
    similarity with TOSTmeta.
-   Added warnings and error messages for possible incorrect input when
    using the functions.

# TOSTER v0.3.3

## (Release date: 2018-05-08)

Changes:

-   Error in order in which p-value 1 and p-value 2 were reported for
    TOSTr due to incorrect use of abs(r) in function. (thanks to Nils
    Kroemer and Dan Quintana)
-   Fixed error in powerTOSTr function - instead of based on conversion
    to d, new function calculates directly from r and is more accurate.
-   Added testthat folder and initial unit tests.

# TOSTER v0.3.2

## (Release date: 2018-04-14)

Changes:

-   Error in sensitivity analyses of powerTOSTpaired.raw,
    powerTOSTone.raw, powerTOSTtwo.raw where output was in Cohen’s d -
    now changed to output in raw scores, added sd to examples. (thanks
    to Lisa DeBruine)
-   Raw power functions still output ceiling N in message, but exact N
    in output value

# TOSTER v0.3.1

## (Release date: 2018-04-06)

Changes:

-   Error in TOSTpaired.raw where t-test for TOST multiplied by sdiff
    (copied from TOSTpaired function) removed. CI were correct, but
    p-values did not match. Now tests are correct. (thanks to
    ontogenerator)
-   Changed text in TOSTr function (text copied from t-test script now
    changed to correlation)
