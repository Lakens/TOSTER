# Test tiers --------
# CRAN:     always runs (analytic tests + fast resampling smoke tests)
# Full:     skip_on_cran(); runs when NOT_CRAN=true (devtools::test(), GitHub CI)
# Extended: skip_if_not_extended(); runs only when TOSTER_EXTENDED_TESTS=true
#           (always set on GitHub CI; opt-in locally)

skip_if_not_extended <- function() {
  testthat::skip_on_cran()
  if (!identical(tolower(Sys.getenv("TOSTER_EXTENDED_TESTS")), "true")) {
    testthat::skip("Extended tests: set TOSTER_EXTENDED_TESTS=true")
  }
}
