NEWS
================

**TOSTER R package and jamovi module**

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
