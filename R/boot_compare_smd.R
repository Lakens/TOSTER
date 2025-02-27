#' @title Comparing Standardized Mean Differences (SMDs) Between Independent Studies with Bootstrapping
#' @description
#' `r lifecycle::badge('stable')`
#'
#' A function to compare standardized mean differences (SMDs) between independent studies using
#' bootstrap methods. This function is intended to be used to compare the compatibility of
#' original studies with replication studies (lower p-values indicating lower compatibility).
#'
#' @param x1 A numeric vector of data values from study 1 (first group for two-sample designs,
#'   or the only group for one-sample/paired designs).
#' @param y1 An optional numeric vector of data values from study 1 (second group for two-sample designs,
#'   or second measurement for paired designs). Set to NULL for one-sample designs.
#' @param x2 A numeric vector of data values from study 2 (first group for two-sample designs,
#'   or the only group for one-sample/paired designs).
#' @param y2 An optional numeric vector of data values from study 2 (second group for two-sample designs,
#'   or second measurement for paired designs). Set to NULL for one-sample designs.
#' @param null A number or vector indicating the null hypothesis value(s):
#'   * For standard tests: a single value representing the null difference (default = 0)
#'   * For equivalence/minimal effect tests: either a single value (symmetric bounds ±value will be created)
#'     or a vector of two values representing the lower and upper bounds
#' @param paired A logical indicating whether the SMD is from a paired or independent samples design.
#'   If a one-sample design, then paired should be set to TRUE.
#' @param alternative A character string specifying the alternative hypothesis:
#'   * "two.sided": difference is not equal to null (default)
#'   * "greater": difference is greater than null
#'   * "less": difference is less than null
#'   * "equivalence": difference is within the equivalence bounds (TOST)
#'   * "minimal.effect": difference is outside the equivalence bounds (TOST)
#'
#'   You can specify just the initial letter.
#' @param R Number of bootstrap replications (default = 1999).
#' @param alpha Alpha level (default = 0.05).
#'
#' @details
#' This function tests for differences between standardized mean differences (SMDs) from
#' independent studies using bootstrap resampling methods. Unlike the `compare_smd` function,
#' which works with summary statistics, this function works with raw data and uses
#' bootstrapping to estimate confidence intervals and p-values.
#'
#' The function supports both paired/one-sample designs and independent samples designs:
#'
#' * For **paired/one-sample designs** (`paired = TRUE`):
#'   * If `y1` and `y2` are provided, the function calculates differences between paired measures
#'   * If `y1` and `y2` are NULL, the function treats `x1` and `x2` as one-sample data
#'   * SMDs are calculated as Cohen's dz (mean divided by standard deviation of differences)
#'
#' * For **independent samples designs** (`paired = FALSE`):
#'   * Requires `x1`, `y1`, `x2`, and `y2` (first and second groups for both studies)
#'   * If `y1` and `y2` are NULL, the function treats `x1` and `x2` as one-sample data with paired = TRUE
#'   * SMDs are calculated as Cohen's ds (mean difference divided by pooled standard deviation)
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
#' The bootstrap procedure follows these steps:
#' 1. Calculate SMDs for both studies using the original data
#' 2. Calculate the difference between SMDs and its standard error
#' 3. Generate R bootstrap samples by resampling with replacement
#' 4. Calculate SMDs and their difference for each bootstrap sample
#' 5. Calculate test statistics for each bootstrap sample
#' 6. Calculate confidence intervals using the percentile method
#' 7. Compute p-values by comparing the observed test statistics to their bootstrap distributions
#'
#' **Note on p-value calculation**: The function uses the bootstrap distribution of test statistics
#' (z-scores) rather than the raw differences to calculate p-values. This approach is analogous to
#' traditional hypothesis testing and estimates the probability of obtaining test statistics as
#' extreme as those observed in the original data under repeated sampling.
#'
#' @return A list with class "htest" containing the following components:
#'
#' * **statistic**: z-score (observed) with name "z (observed)"
#' * **p.value**: The p-value for the test under the null hypothesis
#' * **conf.int**: Bootstrap confidence interval for the difference in SMDs
#' * **estimate**: Difference in SMD between studies
#' * **null.value**: The specified hypothesized value(s) for the null hypothesis
#' * **alternative**: Character string indicating the alternative hypothesis
#' * **method**: Description of the SMD type and design used
#' * **df_ci**: Data frame containing confidence intervals for the difference and individual SMDs
#' * **boot_res**: List containing the bootstrap samples for SMDs, their difference, and test statistics
#' * **data.name**: "Bootstrapped" to indicate bootstrap methods were used
#' * **call**: The matched call
#'
#' @examples
#' # Example 1: Comparing two independent samples SMDs (standard test)
#' set.seed(123)
#' # Study 1 data
#' x1 <- rnorm(30, mean = 0)
#' y1 <- rnorm(30, mean = 0.5, sd = 1)
#' # Study 2 data
#' x2 <- rnorm(25, mean = 0)
#' y2 <- rnorm(25, mean = 0.3, sd = 1)
#'
#' # Two-sided test for independent samples (use fewer bootstraps for example)
#' boot_compare_smd(x1, y1, x2, y2, paired = FALSE,
#'                 alternative = "two.sided", R = 99)
#'
#' # Example 2: Testing for equivalence between SMDs
#' # Testing if the difference between SMDs is within ±0.2
#' boot_compare_smd(x1, y1, x2, y2, paired = FALSE,
#'                 alternative = "equivalence", null = 0.2, R = 99)
#'
#' # Example 3: Testing for minimal effects
#' # Testing if the difference between SMDs is outside ±0.3
#' boot_compare_smd(x1, y1, x2, y2, paired = FALSE,
#'                 alternative = "minimal.effect", null = 0.3, R = 99)
#'
#' # Example 4: Comparing paired samples SMDs
#' # Study 1 data (pre-post measurements)
#' pre1 <- rnorm(20, mean = 10, sd = 2)
#' post1 <- rnorm(20, mean = 12, sd = 2)
#' # Study 2 data (pre-post measurements)
#' pre2 <- rnorm(25, mean = 10, sd = 2)
#' post2 <- rnorm(25, mean = 11, sd = 2)
#'
#' # Comparing paired designs
#' boot_compare_smd(x1 = pre1, y1 = post1, x2 = pre2, y2 = post2,
#'                 paired = TRUE, alternative = "greater", R = 99)
#'
#' # Example 5: Using asymmetric bounds for equivalence testing
#' boot_compare_smd(x1, y1, x2, y2, paired = FALSE,
#'                 alternative = "equivalence", null = c(-0.1, 0.3), R = 99)
#'
#' @family compare studies
#' @name boot_compare_smd
#' @export boot_compare_smd
#'


boot_compare_smd = function(x1,
                            y1 = NULL,
                            x2,
                            y2 = NULL,
                            null = 0,
                            paired = FALSE,
                            alternative = c("two.sided", "less", "greater",
                                            "equivalence", "minimal.effect"),
                            R = 1999,
                            alpha = 0.05){
  alternative <- match.arg(alternative)

  # Parameter validation
  if(!missing(alpha) && (length(alpha) != 1 || !is.finite(alpha) ||
                         alpha < 0 || alpha > 1)) {
    stop("'alpha' must be a single number between 0 and 1")
  }

  # Handle TOST/minimal effects test and validate null parameter
  if(alternative %in% c("equivalence", "minimal.effect")){
    if(length(null) == 1){
      null = c(null, -1*null)
    } else if(length(null) != 2) {
      stop("For equivalence/minimal effect testing, 'null' must be either a single value or a vector of length 2")
    }
    # Ensure lower bound is less than upper bound
    null = c(min(null), max(null))
    TOST = TRUE
  } else {
    if(length(null) > 1){
      stop("null can only have 1 value for non-TOST procedures")
    }
    TOST = FALSE
  }

  # Set confidence level based on alternative
  if(alternative != "two.sided"){
    conf_level = 1 - alpha*2
  } else {
    conf_level = 1 - alpha
  }

  # Prepare data frames for analysis
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

  # Initialize vectors for bootstrap results
  smd1_vec = rep(NA, times=R)
  smd2_vec = rep(NA, times=R)
  d_diff_vec = rep(NA, times=R)
  z_stat_vec = rep(NA, times=R)
  zdiff_stat_vec = rep(NA, times=R)
  z_se_vec = rep(NA, times=R)

  # One-sample or paired design
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

    # Calculate difference and its standard error
    d_diff = smd1 - smd2
    z_se = sqrt(se1^2+se2^2)

    # Bootstrap loop
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

      # Calculate difference and its standard error for bootstrap sample
      d_diff_boot = smd1_boot - smd2_boot
      z_se_boot = sqrt(se1_boot^2 + se2_boot^2)

      # Store bootstrap results
      smd1_vec[i] = smd1_boot
      smd2_vec[i] = smd2_boot
      d_diff_vec[i] = d_diff_boot
      z_se_vec[i] = z_se_boot

      # For comparability of bootstrap distribution
      zdiff_stat_vec[i] = (d_diff_boot - d_diff) / z_se_boot
    }
  } else {  # Two-sample design
    meth = "Bootstrapped Differences in SMDs (two-sample)"
    md1 = mean(subset(df1, group == "x")$y) - mean(subset(df1, group == "y")$y)
    sd1 = poolSD(subset(df1, group == "x")$y, subset(df1, group == "y")$y)
    md2 = mean(subset(df2, group == "x")$y) - mean(subset(df2, group == "y")$y)
    sd2 = poolSD(subset(df2, group == "x")$y, subset(df2, group == "y")$y)

    smd1 = md1/sd1
    smd2 = md2/sd2
    n1_1 = nrow(subset(df1, group == "x"))
    n2_1 = nrow(subset(df1, group == "y"))
    n1_2 = nrow(subset(df2, group == "x"))
    n2_2 = nrow(subset(df2, group == "y"))
    se1 = se_ds(smd1, c(n1_1, n2_1))
    se2 = se_ds(smd2, c(n1_2, n2_2))

    # Calculate difference and its standard error
    d_diff = smd1 - smd2
    z_se = sqrt(se1^2+se2^2)

    # Bootstrap loop
    for(i in 1:R){
      df1_boot = df1[sample(row.names(df1), nrow(df1), replace=TRUE), ]
      df2_boot = df2[sample(row.names(df2), nrow(df2), replace=TRUE), ]

      md1_boot = mean(subset(df1_boot, group == "x")$y) - mean(subset(df1_boot, group == "y")$y)
      sd1_boot = poolSD(subset(df1_boot, group == "x")$y, subset(df1_boot, group == "y")$y)
      md2_boot = mean(subset(df2_boot, group == "x")$y) - mean(subset(df2_boot, group == "y")$y)
      sd2_boot = poolSD(subset(df2_boot, group == "x")$y, subset(df2_boot, group == "y")$y)

      smd1_boot = md1_boot/sd1_boot
      smd2_boot = md2_boot/sd2_boot

      n1_1_boot = nrow(subset(df1_boot, group == "x"))
      n2_1_boot = nrow(subset(df1_boot, group == "y"))
      n1_2_boot = nrow(subset(df2_boot, group == "x"))
      n2_2_boot = nrow(subset(df2_boot, group == "y"))

      se1_boot = se_ds(smd1_boot, c(n1_1_boot, n2_1_boot))
      se2_boot = se_ds(smd2_boot, c(n1_2_boot, n2_2_boot))

      # Calculate difference and its standard error for bootstrap sample
      d_diff_boot = smd1_boot - smd2_boot
      z_se_boot = sqrt(se1_boot^2 + se2_boot^2)

      # Store bootstrap results
      smd1_vec[i] = smd1_boot
      smd2_vec[i] = smd2_boot
      d_diff_vec[i] = d_diff_boot
      z_se_vec[i] = z_se_boot

      # For comparability of bootstrap distribution
      zdiff_stat_vec[i] = (d_diff_boot - d_diff) / z_se_boot
    }
  }

  # Calculate confidence intervals
  smd1_ci = ci_perc(smd1_vec,
                    alternative = alternative,
                    alpha = alpha)
  smd2_ci = ci_perc(smd2_vec,
                    alternative = alternative,
                    alpha = alpha)
  d_diff_ci = ci_perc(d_diff_vec,
                      alternative = alternative,
                      alpha = alpha)

  # Create data frame of confidence intervals
  df_ci = data.frame(estimate = c(d_diff, smd1, smd2),
                     lower.ci = c(d_diff_ci[1], smd1_ci[1], smd2_ci[1]),
                     upper.ci = c(d_diff_ci[2], smd1_ci[2], smd2_ci[2]),
                     row.names = c("Difference in SMD", "SMD1", "SMD2"))


  if(alternative == "equivalence"){
    # For equivalence, test if difference is within bounds
    z_stat_l = (d_diff - min(null)) / z_se  # Test statistic for lower bound
    z_stat_u = (d_diff - max(null)) / z_se  # Test statistic for upper bound


    # Compute p-values comparing observed statistics to bootstrap distribution
    p_l = mean(zdiff_stat_vec <= z_stat_l)  # Proportion less than or equal to lower bound statistic
    p_u = mean(zdiff_stat_vec >= z_stat_u)  # Proportion greater than or equal to upper bound statistic

    pval = max(p_l, p_u)  # Take the maximum (most conservative)
    z_stat = ifelse(p_l > p_u, z_stat_l, z_stat_u)  # Report the less significant statistic

  } else if(alternative == "minimal.effect"){
    # For minimal effects, test if difference is outside bounds
    z_stat_l = (d_diff - min(null)) / z_se  # Test statistic for lower bound
    z_stat_u = (d_diff - max(null)) / z_se  # Test statistic for upper bound



    # Compute p-values comparing observed statistics to bootstrap distribution
    p_l = mean(z_stat_vec >= z_stat_l)  # Proportion greater than or equal to lower bound statistic
    p_u = mean(z_stat_vec <= z_stat_u)  # Proportion less than or equal to upper bound statistic

    pval = min(p_l, p_u)  # Take the minimum (most significant)
    z_stat = ifelse(p_l < p_u, z_stat_l, z_stat_u)  # Report the more significant statistic

  } else {
    # For standard hypothesis tests
    z_stat = (d_diff - null) / z_se  # Test statistic for null hypothesis



    if(alternative == "greater"){
      # Test if difference > null
      pval = mean(zdiff_stat_vec <= z_stat)
    } else if(alternative == "less"){
      # Test if difference < null
      pval = mean(zdiff_stat_vec >= z_stat)
    } else {  # "two.sided"
      # Test if difference ≠ null
      pval = 2 * min(mean(zdiff_stat_vec <= z_stat), mean(zdiff_stat_vec >= z_stat))
      if(pval > 1){
        pval = 1
      }
    }
  }

  # Format output
  names(z_stat) = "z (observed)"
  names(d_diff) = "difference in SMDs"

  # For TOST, report the null range; for standard tests, report single null value
  if(TOST){
    names(null) = rep("difference in SMDs", length(null))
  } else {
    names(null) = "difference in SMDs"
  }

  attr(d_diff_ci, "conf.level") <- conf_level

  # Store as htest object
  rval <- list(statistic = z_stat,
               p.value = pval,
               conf.int = d_diff_ci,
               estimate = d_diff,
               null.value = null,
               alternative = alternative,
               method = meth,
               df_ci = df_ci,
               boot_res = list(
                 smd1 = smd1_vec,
                 smd2 = smd2_vec,
                 d_diff = d_diff_vec,
                 z_stat = z_stat_vec,
                 zdiff_stat = zdiff_stat_vec,
                 z_se = z_se_vec
               ),
               data.name = "Bootstrapped",
               call = match.call())
  class(rval) <- "htest"
  return(rval)
}



