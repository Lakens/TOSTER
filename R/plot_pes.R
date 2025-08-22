#' @title Plot Partial Eta-Squared
#' @description
#' `r lifecycle::badge('stable')`
#'
#' Creates consonance plots (confidence curves and/or consonance density functions) for
#' partial eta-squared values from ANOVA models, allowing visualization of uncertainty
#' around effect size estimates.
#'
#' @param Fstat The F-statistic from the F-test.
#' @param df1 Degrees of freedom for the numerator (effect degrees of freedom).
#' @param df2 Degrees of freedom for the denominator (error degrees of freedom).
#' @param type Choose which plot(s) to create:
#'   * "c": consonance function only (p-values across potential parameter values)
#'   * "cd": consonance density function only (distribution of plausible parameter values)
#'   * c("c", "cd"): both plots together (default)
#' @param levels Numeric vector of confidence levels to display (default: c(.68, .9, .95, .999)).
#'   These correspond to the confidence intervals shown on the plot.
#'
#' @details
#' Consonance plots provide a graphical representation of the full range of confidence intervals
#' for partial eta-squared values at different confidence levels. These plots help visualize
#' the uncertainty around effect size estimates and go beyond the traditional approach of
#' reporting only a single confidence interval (typically 95%).
#'
#' Partial eta-squared (\eqn{\eta^2}) is a measure of effect size commonly used in ANOVA, representing
#' the proportion of variance in the dependent variable attributed to a specific factor,
#' while controlling for other factors in the model. Values range from 0 to 1, with larger
#' values indicating stronger effects.
#'
#' The function creates two types of visualizations:
#'
#' 1. **Consonance function** ("c"): Shows how p-values change across different possible values
#'    of partial eta-squared. The x-axis represents possible parameter values, and the
#'    y-axis represents the corresponding p-values from two-sided hypothesis tests.
#'
#' 2. **Consonance density** ("cd"): Shows the distribution of plausible values for the
#'    partial eta-squared. This can be interpreted as showing where the "weight of evidence"
#'    is concentrated.
#'
#' These plots are particularly useful for:
#' * Visualizing uncertainty around effect size estimates
#' * Understanding the precision of effect size estimates
#' * Comparing the relative plausibility of different effect sizes
#' * Going beyond the binary "significant vs. non-significant" interpretation
#'
#' The required inputs (F-statistic, df1, df2) can typically be extracted from standard
#' ANOVA output in R, such as from `aov()`, `Anova()`, or `afex_aov()` functions.
#'
#' These types of plots are discussed by Schweder & Hjort (2016) and Rafi & Greenland (2020).
#'
#' @return A `ggplot2` object or plot grid from `cowplot`.
#'
#' @examples
#' \dontrun{
#' # Example 1: Basic consonance plot for partial eta-squared
#' # For an F-statistic of 4.5 with df1 = 2, df2 = 60
#' plot_pes(Fstat = 4.5, df1 = 2, df2 = 60)
#'
#' # Example 2: Consonance function only (p-value curve)
#' plot_pes(Fstat = 3.2, df1 = 1, df2 = 45, type = "c")
#'
#' # Example 3: Consonance density only
#' plot_pes(Fstat = 6.8, df1 = 3, df2 = 80, type = "cd")
#'
#' # Example 4: Custom confidence levels
#' plot_pes(Fstat = 5.1, df1 = 2, df2 = 50, levels = c(0.5, 0.8, 0.95))
#'
#' # Example 5: Using with actual ANOVA results
#' # aov_result <- aov(DV ~ IV, data = your_data)
#' # aov_summary <- summary(aov_result)[[1]]
#' # F_value <- aov_summary$"F value"[1]
#' # df1 <- aov_summary$Df[1]
#' # df2 <- aov_summary$Df[2]
#' # plot_pes(Fstat = F_value, df1 = df1, df2 = df2)
#'
#' # Example 6: Saving and further customizing the plot
#' library(ggplot2)
#' p <- plot_pes(Fstat = 4.5, df1 = 2, df2 = 60)
#' p + theme_minimal() +
#'   labs(title = "Consonance Plot for Partial Eta-Squared",
#'        subtitle = "F(2, 60) = 4.5")
#'}
#' @references
#' Schweder, T., & Hjort, N. L. (2016). Confidence, likelihood, probability:
#' Statistical inference with confidence distributions. Cambridge University Press.
#' ISBN: 9781316445051
#'
#' Rafi, Z., & Greenland, S. (2020). Semantic and cognitive tools to aid statistical science:
#' Replace confidence and significance by compatibility and surprise. BMC Medical Research
#' Methodology, 20, 244. doi:10.1186/s12874-020-01105-9
#'
#' @family plotting functions
#' @export

plot_pes <- function(Fstat,
                     df1,
                     df2,
                     type = c("c","cd"),
                     levels = c(.68,.9,.95,.999)){

  dat = pes_curv(Fstat,
                 df1,
                 df2,
                 steps = 5000)

  resplot = gg_curv_t(
    dat,
    type = type,
    levels = levels,
    position = "pyramid",
    xaxis = expression(paste("partial ",eta^2)),
    yaxis1 = expression(paste(italic(p),
                              "-value")),
    yaxis2 = "Confidence Interval (%)",
    color = "black",
    fill = "skyblue",
    alpha_shade = .5
  )



  return(resplot)

}
