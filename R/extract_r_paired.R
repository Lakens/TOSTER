#' @title Extract Paired Correlation
#' @description A function for estimating the correlation from a paired samples t-test.
#' Useful for tsum_TOST when the correlation is not available.
#' @param tstat The t-value from a paired samples t-test
#' @param pvalue The two-tailed p-value from a paired samples t-test
#' @param n Sample size (number of pairs)
#' @inheritParams tsum_TOST
#' @return An estimate of the correlation.
#' @references
#' Lajeunesse, M. J. (2011). On the meta‐analysis of response ratios for studies with correlated and multi‐group designs. Ecology, 92(11), 2049-2055
#' @importFrom stats na.omit setNames terms
#' @export
extract_r_paired = function(m1,
                            sd1,
                            m2,
                            sd2 = NULL,
                            n,
                            tstat = NULL,
                            pvalue = NULL){

  if(is.null(sd2)){
    sd2 = sd1
  }
  if(is.null(tstat)){
    if(is.null(pvalue)){
      stop("tstat or pvalue must be provided")
    }

    tstat = qt(pvalue/2, n-1, lower.tail=FALSE)

  }

  corr = (sd2^2 + sd1^2 - tstat^(-2)*n*(m1-m2)^2) / (2*sd2*sd1)

  return(corr)
}
