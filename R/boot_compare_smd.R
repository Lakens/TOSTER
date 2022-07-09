#' @title Comparing SMDs between independent studies with Bootstrapping
#' @description A function to compare standardized mean differences (SMDs) between studies. This function is intended to be used to compare the compatibility of original studies with replication studies (lower p-values indicating lower compatibility)
#' @param x1 	a (non-empty) numeric vector of data values from study 1.
#' @param y1 an optional (non-empty) numeric vector of data values from study 1.
#' @param x2 a (non-empty) numeric vector of data values from study 2.
#' @param y2 an optional (non-empty) numeric vector of data values from study 2.
#' @param null a number indicating the null hypothesis. For TOST, this would be equivalence bound.
#' @param paired a logical indicating whether the SMD is from a paired or independent samples design.
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less". You can specify just the initial letter.
#' @param R number of bootstrap replicates
#' @param alpha alpha level (default = 0.05)
#' @return A list with class "htest" containing the following components:
#' \describe{
#'   \item{\code{"statistic"}}{z-score}
#'   \item{\code{"p.value"}}{numeric scalar containing the p-value for the test under the null hypothesis.}
#'   \item{\code{"estimate"}}{difference in SMD between studies}
#'   \item{\code{"conf.int}}{percentile (bootstrap) confidence interval for difference in SMDs}
#'   \item{\code{"null.value"}}{the specified hypothesized value for the null hypothesis.}
#'   \item{\code{"alternative"}}{character string indicating the alternative hypothesis (the value of the input argument alternative). Possible values are "greater", "less", or "two-sided".}
#'   \item{\code{"method"}}{Type of SMD}
#'   \item{\code{"data.name"}}{"Boostrapped" to denote summary statistics were utilized to obtain results.}
#'   \item{\code{"smd"}}{SMDs input for the function.}
#'   \item{\code{"df_ci"}}{Data frame of confidence intervals.}
#'   \item{\code{"boot_res"}}{List of bootstrapped results.}
#'   \item{\code{"call"}}{the matched call.}
#' }
#' @name boot_compare_smd
#' @export boot_compare_smd
#'


boot_compare_smd = function(x1,
                            y1 = NULL,
                            x2,
                            y2 = NULL,
                            null = 0,
                            paired = FALSE,
                            alternative = c("two.sided", "less", "greater"),
                            R = 1999,
                            alpha = 0.05){
  alternative <- match.arg(alternative)
  if(!missing(null) && (length(null) != 1 || is.na(null))) {
    stop("'null' must be a single number")
  }


  if(!missing(alpha) && (length(alpha) != 1 || !is.finite(alpha) ||
                         alpha < 0 || alpha > 1)) {
    stop("'alpha' must be a single number between 0 and 1")
  }
  if(alternative == "two.sided"){
    conf_level = 1-alpha
  } else {
    conf_level = 1-alpha*2
  }



  #if(is.null(y1) && is.null(y2)){
  #  paired = TRUE
  #}

  if(paired){
    if(is.null(y1) && is.null(y2)){
      df1 = na.omit(data.frame(z = x1))
      df2 = na.omit(data.frame(z = x2))
    } else {
      df1 = na.omit(data.frame(z = x1 - y1))
      df2 = na.omit(data.frame(z = x2 - y2))
    }
  } else if(is.null(y1) && is.null(y2)){

      df1 = na.omit(data.frame(z = x1))
      df2 = na.omit(data.frame(z = x2))

  } else {
    x1 = na.omit(x1)
    y1 = na.omit(y1)
    x2 = na.omit(x2)
    y2 = na.omit(y2)
    df1 = data.frame(y = c(x1,y1),
                     group = c(rep("x",length(x1)),
                               rep("y",length(y1))))
    df2 = data.frame(y = c(x2,y2),
                     group = c(rep("x",length(x2)),
                               rep("y",length(y2))))
  }


  smd1_vec = rep(NA, times=length(R))
  smd2_vec = rep(NA, times=length(R))
  d_diff_vec = rep(NA, times=length(R))
  z_stat_vec = rep(NA, times=length(R))
  zdiff_stat_vec = rep(NA, times=length(R))
  #m_vec <- rep(NA, times=length(R)) # mean difference vector
  if(ncol(df1) == 1){
    if(paired){
      meth = "Bootstrapped Differences in SMDs (paired)"
    } else {
      meth = "Bootstrapped Differences in SMDs (one-sample)"
    }

    md1 = mean(df1$z)
    sd1 = sd(df1$z)
    md2 = mean(df2$z)
    sd2 = sd(df2$z)

    smd1 = md1/sd1
    smd2 = md2/sd2
    se1 = se_dz(smd1, length(df1$z))
    se2 = se_dz(smd2, length(df2$z))
    d_diff = smd1 - smd2
    z_se = sqrt(se1^2+se2^2)
    z_stat = d_diff/z_se
    for(i in 1:R){
      df1_boot = df1[sample(row.names(df1), nrow(df1), replace=TRUE), ]
      df2_boot = df2[sample(row.names(df2), nrow(df2), replace=TRUE), ]

      md1_boot = mean(df1_boot)
      sd1_boot = sd(df1_boot)
      md2_boot = mean(df2_boot)
      sd2_boot = sd(df2_boot)

      smd1_boot = md1_boot/sd1_boot
      smd2_boot = md2_boot/sd2_boot
      se1_boot = se_dz(smd1_boot, length(df1_boot))
      se2_boot = se_dz(smd2_boot, length(df2_boot))
      d_diff_boot = smd1_boot - smd2_boot
      z_se_boot = sqrt(se1_boot^2 + se2_boot^2)
      z_stat_boot = d_diff_boot / z_se_boot
      zdiff_stat_boot = (d_diff_boot - d_diff) / z_se_boot

      smd1_vec[i] = smd1_boot
      smd2_vec[i] = smd2_boot
      d_diff_vec[i] = d_diff_boot
      z_stat_vec[i] = z_stat_boot
      zdiff_stat_vec[i] = zdiff_stat_boot
    }
  } else{
    meth = "Bootstrapped Differences in SMDs (two-sample)"
    md1 = mean(subset(df1,
                      group == "x")$y) -
      mean(subset(df1,
                  group == "y")$y)
    sd1 = poolSD(subset(df1,
                        group == "x")$y,
                 subset(df1,
                        group == "y")$y)
    md2 = mean(subset(df2,
                      group == "x")$y) -
      mean(subset(df2,
                  group == "y")$y)
    sd2 = poolSD(subset(df2,
                        group == "x")$y,
                 subset(df2,
                        group == "y")$y)

    smd1 = md1/sd1
    smd2 = md2/sd2
    se1 = se_ds(smd1, nrow(df1))
    se2 = se_ds(smd2, nrow(df2))
    d_diff = smd1 - smd2
    z_se = sqrt(se1^2+se2^2)
    z_stat = d_diff/z_se
    for(i in 1:R){
      df1_boot = df1[sample(row.names(df1), nrow(df1), replace=TRUE), ]
      df2_boot = df2[sample(row.names(df2), nrow(df2), replace=TRUE), ]
      md1_boot = mean(subset(df1_boot,
                        group == "x")$y) -
        mean(subset(df1_boot,
                    group == "y")$y)
      sd1_boot = poolSD(subset(df1_boot,
                          group == "x")$y,
                   subset(df1_boot,
                          group == "y")$y)
      md2_boot = mean(subset(df2_boot,
                        group == "x")$y) -
        mean(subset(df2_boot,
                    group == "y")$y)
      sd2_boot = poolSD(subset(df2_boot,
                          group == "x")$y,
                   subset(df2_boot,
                          group == "y")$y)

      smd1_boot = md1_boot/sd1_boot
      smd2_boot = md2_boot/sd2_boot
      se1_boot = se_ds(smd1_boot, nrow(df1_boot))
      se2_boot = se_ds(smd2_boot, nrow(df2_boot))
      d_diff_boot = smd1_boot - smd2_boot
      z_se_boot = sqrt(se1_boot^2+se2_boot^2)
      z_stat_boot = d_diff_boot/z_se_boot
      zdiff_stat_boot = (d_diff_boot - d_diff) / z_se_boot

      smd1_vec[i] = smd1_boot
      smd2_vec[i] = smd2_boot
      d_diff_vec[i] = d_diff_boot
      z_stat_vec[i] = z_stat_boot
      zdiff_stat_vec[i] = zdiff_stat_boot
    }
  }

  smd1_ci = ci_perc(smd1_vec,
                    alternative = alternative,
                    alpha = alpha)
  smd2_ci = ci_perc(smd2_vec,
                    alternative = alternative,
                    alpha = alpha)
  d_diff_ci = ci_perc(d_diff_vec,
                    alternative = alternative,
                    alpha = alpha)
  df_ci = data.frame(estimate = c(d_diff, smd1, smd2),
                     lower.ci = c(d_diff_ci[1],smd1_ci[1],smd2_ci[1]),
                     upper.ci = c(d_diff_ci[2],smd1_ci[2],smd2_ci[2]),
                     row.names = c("Difference in SMD", "SMD1", "SMD2"))
  # Calculate p-value
  if(alternative == "greater"){
    pval = sum(zdiff_stat_vec >= z_stat)/length(zdiff_stat_vec)
  } else if(alternative == "less"){
    pval = sum(zdiff_stat_vec <= z_stat)/length(zdiff_stat_vec)
  } else {
    pval1 = sum(zdiff_stat_vec >= z_stat)/length(zdiff_stat_vec)
    pval2 = sum(zdiff_stat_vec <= z_stat)/length(zdiff_stat_vec)
    pval = 2*min(pval1,pval2)
    if(pval > 1){
      pval = 1
    }
  }
  par = R
  names(par) = "R"
  names(z_stat) = "z (observed)"
  names(d_diff) = "difference in SMDs"
  names(null) = "difference in SMDs"
  attr(d_diff_ci,"conf.level") <- conf_level
  # Store as htest
  rval <- list(statistic = z_stat, p.value = pval,
               conf.int = d_diff_ci,
               estimate = d_diff,
               null.value = null,
               alternative = alternative,
               method = meth,
               df_ci = df_ci,
               boot_res = list(
                 smd1 = smd1_vec,
                 d_diff = d_diff_vec,
                 z_stat = z_stat_vec,
                 zdiff_stat = zdiff_stat_vec
               ),
               data.name = "Bootstrapped",
               call = match.call())
  class(rval) <- "htest"
  return(rval)
}



