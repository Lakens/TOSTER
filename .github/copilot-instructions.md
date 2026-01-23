# GitHub Copilot Instructions - TOSTER Package

## Package Context

TOSTER implements equivalence testing (TOST procedure) for t-tests, correlations, proportions, and meta-analyses in R. Uses roxygen2 with markdown enabled.

## Code Style Requirements

- Follow existing patterns in the repository
- Function names: `snake_case` for new functions (e.g., `t_TOST`, `boot_t_TOST`)
- Arguments: lowercase with dots or underscores (e.g., `low_eqbound`, `var.equal`)
- Use `match.arg()` for parameter validation with predefined choices
- Maintain consistency with surrounding code in all suggestions

## For Issue Fixes

When tagged to fix an issue:
- Make minimal, targeted changes addressing only the specific problem
- Do not refactor unrelated code
- Preserve existing function signatures unless the issue requires changing them
- If the fix needs a new dependency, flag for human review instead of adding it
- Check if similar patterns exist elsewhere in the codebase and follow them

## For PR Reviews

Check for:
- **Style consistency**: Matches existing code patterns
- **Return classes**: Test functions should return `htest` or `TOSTt`/`TOSTnp`; power functions should return `power.htest`
- **Documentation**: Complete roxygen2 with `@param`, `@return`, `@examples`, and `@family` tags
- **Test coverage**: New code paths need corresponding tests
- **Hardcoded values**: Should be parameterized where appropriate
- **Edge cases**: NA handling, zero-length inputs, type validation
- **Input validation**: Numeric bounds checking, argument type checks

## Limitations - Do Not Attempt

Without human review, do not:
- Modify C++ code in `src/`
- Change GitHub Actions workflows
- Add or remove package dependencies
- Refactor S3 method dispatch patterns
- Change exported function APIs (breaking changes)
- Modify jamovi interface files (`R/*jmv*.R`, `R/data*.R`)

If a fix requires any of the above, describe what's needed and recommend human review.

## Testing Expectations

- Every code fix should include a suggested test case
- Use testthat (edition 3)
- Numerical comparisons: `expect_equal(result, expected, tolerance = 1e-5)`
- Reference known published results when available for validation
- Test both success paths and error conditions (`expect_error()`)

## When Uncertain

- If multiple approaches exist, describe options rather than picking one
- Flag anything affecting backward compatibility
- Recommend human review for changes touching more than 2-3 files
- Note if a fix might interact with the jamovi interface
