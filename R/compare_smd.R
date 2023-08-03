#' @title Comparing SMDs between independent studies
#' @description
#' `r lifecycle::badge('stable')`
#'
#' A function to compare standardized mean differences (SMDs) between studies. This function is intended to be used to compare the compatibility of original studies with replication studies (lower p-values indicating lower compatibility).
#'
#'
#' @param smd1,smd2 SMDs from study 1 & 2, respectively.
#' @param n1,n2 sample size(s) from study 1 & 2, respectively (can be 1 number or vector of 2 numbers).
#' @param se1,se2 User supplied standard errors (SEs). This will override the internal calculations.
#' @param paired a logical indicating whether the SMD is from a paired or independent samples design. If a one-sample design, then paired must be set to TRUE.
#' @param null a number indicating the null hypothesis. For TOST, this would be equivalence bound.
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater", "less", "equivalence" (TOST), or "minimal.effect" (TOST). You can specify just the initial letter.
#' @param TOST Defunct: use alternative argument. Logical indicator (default = FALSE) to perform two one-sided tests of equivalence (TOST).
#'
#' @details This function tests for differences between SMDs from independent studies (e.g., original vs replication).
#'
#' @return A list with class "htest" containing the following components:
#'
#'   - "statistic": z-score.
#'   - "p.value": numeric scalar containing the p-value for the test under the null hypothesis.
#'   - "estimate": difference in SMD between studies.
#'   - "null.value": the specified hypothesized value for the null hypothesis.
#'   - "alternative": character string indicating the alternative hypothesis (the value of the input argument alternative). Possible values are "greater", "less", or "two-sided".
#'   - "method": Type of SMD.
#'   - "data.name": "Summary Statistics" to denote summary statistics were utilized to obtain results.
#'   - "smd": SMDs input for the function.
#'   - "sample_sizes": Sample sizes input for the function.
#'   - "call": the matched call.
#'
#' @name compare_smd
#' @family compare studies
#' @export compare_smd
#'

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
