#' @title Comparing Standardized Mean Differences (SMDs) Between Independent Studies
#' @description
#' `r lifecycle::badge('stable')`
#'
#' A function to compare standardized mean differences (SMDs) between independent studies.
#' This function is intended to be used to compare the compatibility of original studies
#' with replication studies (lower p-values indicating lower compatibility).
#'
#' @param smd1,smd2 SMDs from study 1 & 2, respectively.
#' @param n1,n2 Sample size(s) from study 1 & 2, respectively. Can be a single number (total sample size)
#'   or a vector of 2 numbers (group sizes) for independent samples designs.
#' @param se1,se2 User supplied standard errors (SEs). This will override the internal calculations for the standard error.
#' @param paired A logical indicating whether the SMD is from a paired or independent samples design.
#'   If a one-sample design, then paired must be set to TRUE.
#' @param alternative A character string specifying the alternative hypothesis:
#'   * "two.sided": difference is not equal to null (default)
#'   * "greater": difference is greater than null
#'   * "less": difference is less than null
#'   * "equivalence": difference is within the equivalence bounds (TOST)
#'   * "minimal.effect": difference is outside the equivalence bounds (TOST)
#'
#'   You can specify just the initial letter.
#' @param null A number or vector indicating the null hypothesis value(s):
#'   * For standard tests: a single value representing the null difference (default = 0)
#'   * For equivalence/minimal effect tests: either a single value (symmetric bounds ±value will be created)
#'     or a vector of two values representing the lower and upper bounds
#' @param TOST Defunct: use alternative argument. Logical indicator (default = FALSE) to perform
#'   two one-sided tests of equivalence (TOST).
#'
#' @details
#' This function tests for differences between SMDs from independent studies (e.g., original vs replication).
#' It is particularly useful for:
#'
#' * Comparing effect sizes between an original study and its replication
#' * Meta-analytic comparisons between studies
#' * Testing if effect sizes from different samples are equivalent
#'
#' The function handles both paired/one-sample designs and independent samples designs:
#'
#' * For paired/one-sample designs (`paired = TRUE`), standard errors are calculated
#'   for Cohen's dz, and n1 and n2 must be single values.
#'
#' * For independent samples designs (`paired = FALSE`), standard errors are calculated
#'   for Cohen's ds, and n1 and n2 can be either single values (total sample size) or
#'   vectors of length 2 (group sizes).
#'
#' * For all other SMDs, you should supply your own standard errors using the `se1` and `se2` arguments.
#'
#' The function supports both standard hypothesis testing and equivalence/minimal effect testing:
#'
#' * For standard tests (two.sided, less, greater), the function tests whether the difference
#'   between SMDs differs from the null value (typically 0).
#'
#' * For equivalence testing ("equivalence"), it determines whether the difference falls within
#'   the specified bounds, which can be set asymmetrically.
#'
#' * For minimal effect testing ("minimal.effect"), it determines whether the difference falls
#'   outside the specified bounds.
#'
#' When performing equivalence or minimal effect testing:
#' * If a single value is provided for `null`, symmetric bounds ±value will be used
#' * If two values are provided for `null`, they will be used as the lower and upper bounds
#'
#' @return A list with class "htest" containing the following components:
#'
#' * **statistic**: z-score with name "z"
#' * **p.value**: numeric scalar containing the p-value for the test under the null hypothesis
#' * **estimate**: difference in SMD between studies
#' * **null.value**: the specified hypothesized value(s) for the null hypothesis
#' * **alternative**: character string indicating the alternative hypothesis
#' * **method**: description of the method used for comparison
#' * **data.name**: "Summary Statistics" to denote summary statistics were utilized
#' * **smd**: list containing the SMDs used in the comparison
#' * **sample_sizes**: list containing the sample sizes used in the comparison
#' * **call**: the matched call
#'
#' @examples
#' # Example 1: Comparing two independent samples SMDs (standard test)
#' compare_smd(smd1 = 0.5, n1 = c(30, 30),
#'             smd2 = 0.3, n2 = c(25, 25),
#'             paired = FALSE, alternative = "two.sided")
#'
#' # Example 2: Comparing two paired samples SMDs
#' compare_smd(smd1 = 0.6, n1 = 40,
#'             smd2 = 0.4, n2 = 45,
#'             paired = TRUE, alternative = "two.sided")
#'
#' # Example 3: Testing for equivalence between SMDs
#' # Testing if the difference between SMDs is within ±0.2
#' compare_smd(smd1 = 0.45, n1 = c(25, 25),
#'             smd2 = 0.35, n2 = c(30, 30),
#'             paired = FALSE, alternative = "equivalence", null = 0.2)
#'
#' # Example 4: Testing for minimal effects
#' # Testing if the difference between SMDs is outside ±0.3
#' compare_smd(smd1 = 0.7, n1 = 30,
#'             smd2 = 0.3, n2 = 35,
#'             paired = TRUE, alternative = "minimal.effect", null = 0.3)
#'
#' # Example 5: Using asymmetric bounds for equivalence testing
#' compare_smd(smd1 = 0.45, n1 = c(30, 30),
#'             smd2 = 0.35, n2 = c(25, 25),
#'             paired = FALSE, alternative = "equivalence", null = c(-0.1, 0.3))
#'
#' # Example 6: Using user-supplied standard errors
#' compare_smd(smd1 = 0.5, n1 = 50, se1 = 0.15,
#'             smd2 = 0.7, n2 = 45, se2 = 0.16,
#'             paired = TRUE, alternative = "two.sided")
#'
#' @name compare_smd
#' @family compare studies
#' @export compare_smd

compare_smd = function(smd1,
                       n1,
                       se1 = NULL,
                       smd2,
                       n2,
                       se2 = NULL,
                       paired = FALSE,
                       alternative = c("two.sided", "less", "greater",
                                       "equivalence", "minimal.effect"),
                       null = 0,
                       TOST = FALSE){
  alternative <- match.arg(alternative)
  if(TOST && alternative %in% c("two.sided","greater","less")){
    alternative = "equivalence"
  }
  if(alternative %in% c("equivalence", "minimal.effect")){
    if(length(null) == 1){
      null = c(null, -1*null)
    }
    TOST = TRUE
  } else {
    if(length(null) > 1){
      stop("null can only have 1 value for non-TOST procedures")
    }
    TOST = FALSE
  }
  if(length(n1) > 2 || length(n2) > 2 || !is.numeric(n1) || !is.numeric(n1)){
    stop("n1 and n2 must be a numeric vector of a length of 1 or 2.")
  }
  if(paired && (length(n1) > 1 || length(n2) > 1)){
    stop("n1 and n2 must be a length of 1 if paired is TRUE.")
  }



  if(!is.null(se1) || !is.null(se2)){
    message("User supplied standard errors. Proceed with caution.")
  }

  # calculate standard errors
  if(paired){
    if(!is.null(se1) || !is.null(se2)){
      meth = "Difference in SMDs (user-supplied SE)"
    } else{
      meth = "Difference in Cohen's dz (paired)"
    }
    if(is.null(se1)){
      se1 = se_dz(smd1, n1)
    } else{
      se1 = se1
    }
    if(is.null(se2)){
      se2 = se_dz(smd2, n2)
    } else{
      se2 = se2
    }

  } else{
    if(!is.null(se1) || !is.null(se2)){
      meth = "Difference in SMDs (user-supplied SE)"
    } else{
      meth = "Difference in Cohen's ds (two-sample)"
    }
    if(is.null(se1)){
      se1 = se_ds(smd1, n1)
    } else{
      se1 = se1
    }
    if(is.null(se2)){
      se2 = se_ds(smd2, n2)
    } else{
      se2 = se2
    }

  }

  # z-score and p-value
  # difference in SMD minus null hypothesis
  if(TOST){
    se_diff = sqrt(se1^2 + se2^2)
    nullhi = max(null)
    nulllo = min(null)
    d_difflo = smd1 - smd2 - nulllo
    zlo = d_difflo/se_diff
    d_diffhi = smd1 - smd2 - nullhi
    zhi = d_diffhi/se_diff
    if(alternative == "equivalence"){
      plo = p_from_z(zlo, alternative = 'greater')
      phi = p_from_z(zhi, alternative = 'less')
      if(phi >= plo){
        pvalue = phi
        z_test = zhi
      } else {
        pvalue = plo
        z_test = zlo
      }
    }
    if(alternative == "minimal.effect"){
      plo = p_from_z(zlo, alternative = 'less')
      phi = p_from_z(zhi, alternative = 'greater')
      if(phi <= plo){
        pvalue = phi
        z_test = zhi
      } else {
        pvalue = plo
        z_test = zlo
      }
    }

    z = z_test
    names(z) = "z"
    pval = pvalue

  }else {
    d_diff = smd1 - smd2 - null
    se_diff = sqrt(se1^2 + se2^2)
    z = d_diff/se_diff
    names(z) = "z"
    pval = p_from_z(z, alternative = alternative)
  }


  # Equivalence Testing
  #if(TOST){
  #  d_diff2 = abs(smd1 - smd2) - null
  #  z2 = d_diff2/se_diff
  #  if(abs(z2) > abs(z)){
  #    z = z2
  #    names(z) = "z"
  #  }
  #  pval = p_from_z(z2, alternative = "less")
  #  alternative = "less"
  #}

  est2 = smd1 - smd2
  names(est2) = "difference in SMDs"
  null2 = null
  names(null2) = rep("difference in SMDs",length(null2))
  # Store as htest
  rval <- list(statistic = z, p.value = pval,
               #conf.int = cint,
               estimate = est2,
               null.value = null2,
               alternative = alternative,
               method = meth,
               smd = list(smd1 = smd1,
                          smd2 = smd2),
               sample_sizes = list(n1 = n1,
                                   n2 = n2),
               data.name = "Summary Statistics",
               call = match.call())
  class(rval) <- "htest"
  return(rval)
}
