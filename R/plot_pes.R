#' Function to produce plots of the distribution of the standardized mean difference
#' @param Fstat The F-statistic from the F-test.
#' @param df1 Degrees of freedom for the numerator.
#' @param df2 Degrees of freedom for the denominator.
#' @param type Choose whether to plot a "consonance" function ("c"), consonance density ("cd"), or both (c("c","cd"); defualt option).
#' @param levels Numeric vector of confidence levels to display
#' @details
#' This function was created so that users could create consonance plots of partial eta-squared from ANOVA-level effects.
#' These types of plots are discussed by Schweder T, Hjort NL. (2016, ISBN:9781316445051) and Rafi Z, Greenland S. (2020) <doi:10.1186/s12874-020-01105-9>.
#' @return Returns plot of the distribution of partial eta-squared
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
