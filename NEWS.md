NEWS
================

**TOSTER R package and jamovi module**

# TOSTER v0.9.0
- Added `perm_t_test` function to allow for permutation tests

# TOSTER v0.8.7
- Update documentation to make it clear what the "eqb" argument does within the `wilcox_TOST` function.

# TOSTER v0.8.6
- Add warning message about error control with MET 
- Fix unit tests for boot_ses_calc to catch errors when estimates contain infinite values or when ses is not "rb"

# TOSTER v0.8.5
- Big update to package documentation to make things more detailed.
- Added extra message when permutation tests are used for the Brunner-Munzel test.
- Added more alternative hypotheses to `boot_compare_smd`.
- Expanded functionality of `power_eq_f` function.

# TOSTER v0.8.4
- Added simple plot for `TOSTt` methods
- Small fix to output for printed method for `TOSTt`
- Added fix to `power_z_cor` to provide "zero" power for scenarios where the effect is undetectable

# TOSTER v0.8.3
- Change in the standard error formulation to Glass delta for independent samples
  - Hat tip to Paul Dudgeon for catching an error in the code that led to this development
  - Other small changes to standard error calcs (see vignettes for details)
- Small modification to `plot_smd` to catch errors when attempting to plot
- Brunner-Munzel updates
  - Added more warning messages
  - Changed default for `simple_htest` to 0.5 rather than 0

# TOSTER v0.8.2

- Fixed error with `describe` method for minimal effects test for `TOSTt` objects.

# TOSTER v0.8.1

- Small correction to the displayed equation for Cohen's ds standard error. Thank you to Matthew B Jané for finding this error.
- Added bootstrap options such as `boot_smd_calc` and `boot_ses_calc`.
  - Many functions also now allow for different CI methods for bootstrapped results.

# TOSTER v0.8.0

- Added Brunner-Munzel test
- Updated documentation to include lifecycle labels
- Created new function for two proportions tests (`twoprop_test`)
  - And power `power_twoprop`
- Created new function for power for correlations (`power_z_cor`)
- Deprecated old functions

# TOSTER v0.7.1

- Fixing the ggplot2 error related to `after_stat` update from that package.
- Update documentation.

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
