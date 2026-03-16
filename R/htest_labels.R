# Internal label helpers for htest-returning functions --------
#
# Shared utilities that build names(estimate) labels and SMD notation
# strings across TOSTER's t-test, Wilcoxon, and SMD functions.
# The core principle: default methods use generic group names ("x"/"y");
# formula methods resolve actual factor level names; numeric-looking
# names get quoted via quote_if_numeric() (defined in trans_rank_prob.R).

# Internal: Build estimate labels for t-test style htest objects --------
#
# Constructs names(estimate) labels for t-test, Wilcoxon, bootstrap,
# and permutation test functions. Follows the same XNAME/YNAME resolution
# pattern as brunner_munzel().
#
# @param type character; one of "t" (t-test), "wilcoxon" (Wilcoxon/Mann-Whitney)
# @param xname character; label for the first group (default: "x")
# @param yname character or NULL; label for the second group (default: "y", NULL for one-sample)
# @param paired logical; whether the test is paired
# @return character vector of label(s) for names(estimate)
# @noRd
ttest_estimate_label <- function(type = c("t", "wilcoxon"),
                                 xname = "x",
                                 yname = "y",
                                 paired = FALSE) {
  type <- match.arg(type)

  xq <- quote_if_numeric(xname)

  if (is.null(yname)) {
    # One-sample
    qty <- switch(type,
      "t"        = "mean",
      "wilcoxon" = "(pseudo)median"
    )
    return(paste(qty, "of", xq))
  }

  yq <- quote_if_numeric(yname)

  if (type == "t") {
    if (paired) {
      # Paired t-test: single estimate (test on difference scores z = x - y)
      paste0("mean of the differences (z = ", xq, " - ", yq, ")")
    } else {
      # Two-sample t-test: two estimates (group means)
      c(paste("mean of group", xq),
        paste("mean of group", yq))
    }
  } else {
    # Wilcoxon: single Hodges-Lehmann estimate of location shift
    if (paired) {
      # Paired: test on difference scores z = x - y
      paste0("Hodges-Lehmann estimate (z = ", xq, " - ", yq, ")")
    } else {
      # Two-sample
      paste0("Hodges-Lehmann estimate (", xq, " - ", yq, ")")
    }
  }
}

# Internal: Build SMD notation label --------
#
# Constructs the ((XNAME-YNAME)/SD_*) notation string for SMD method descriptions.
# Used by smd_calc() and boot_smd_calc().
#
# @param xname character; label for the first group (default: "x")
# @param yname character or NULL; label for the second group (NULL for one-sample)
# @param denom_label character; the SD subscript label (e.g., "pooled", "z", or a group name for Glass)
# @return character string like "((x-y)/SD_pooled)"
# @noRd
smd_notation_label <- function(xname = "x",
                               yname = NULL,
                               denom_label = "pooled") {
  xq <- quote_if_numeric(xname)

  if (is.null(yname)) {
    # One-sample: numerator is just the variable
    numerator <- xq
  } else {
    yq <- quote_if_numeric(yname)
    numerator <- paste0(xq, "-", yq)
  }

  paste0("(", numerator, ")/SD_", denom_label)
}

# Internal: Resolve SD subscript label for SMD notation --------
#
# Maps the resolved denominator state to a human-readable subscript.
#
# @param denom character; the user-facing denom argument value
# @param int_denom character; the internally resolved denominator code ("z", "rm", "d", "glass1", "glass2")
# @param xname character; first group label (used for glass1)
# @param yname character; second group label (used for glass2)
# @return character string for the SD subscript
# @noRd
resolve_sd_label <- function(denom, int_denom, xname = "x", yname = "y",
                             var.equal = TRUE) {
  if (denom != "auto") {
    # User explicitly chose denom - use it directly (except glass -> group name)
    switch(denom,
      "z"      = "z",
      "rm"     = "rm",
      "pooled" = "pooled",
      "avg"    = "avg",
      "glass1" = quote_if_numeric(xname),
      "glass2" = quote_if_numeric(yname)
    )
  } else {
    # Auto-resolved - map from int_denom
    switch(int_denom,
      "z"      = "z",
      "rm"     = "rm",
      "d"      = if (var.equal) "pooled" else "avg",
      "glass1" = quote_if_numeric(xname),
      "glass2" = quote_if_numeric(yname)
    )
  }
}
