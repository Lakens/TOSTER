#' @title Plot correlation coefficients
#' @description
#' `r lifecycle::badge('stable')`
#'
#' Function to produce plots of the distribution of standard correlation coefficients
#' @param r The observed correlation coefficient.
#' @param n Total number of observations (sample size).
#' @param method The method by which the coefficient was calculated: pearson, spearman, or kendall (default is "pearson")
#' @param type Choose whether to plot a "consonance" function ("c"), consonance density ("cd"), or both (c("c","cd"); defualt option).
#' @param levels Numeric vector of confidence levels to display
#' @details
#' This function was created so that users could create consonance plots of Pearson's correlation coefficient.
#' These types of plots are discussed by Schweder T, Hjort NL. (2016, ISBN:9781316445051) and Rafi Z, Greenland S. (2020) <doi:10.1186/s12874-020-01105-9>.
#' @return Returns plot of the distribution of the correlation coefficient.
#' @family plotting functions
#' @export

plot_cor <- function(r,
                     n,
                     method = c("pearson","spearman","kendall"),
                     type = c("c", "cd"),
                     levels = c(.68, .9, .95, .999)){
  method = match.arg(method)
  dat = corr_curv(r = r,
                  n = n,
                  type = method,
                  steps = 5000)

  resplot = gg_curv_t(
    dat,
    type = type,
    levels = levels,
    position = "pyramid",
    xaxis = "Correlation Coefficent",
    yaxis1 = expression(paste(italic(p),
                              "-value")),
    yaxis2 = "Confidence Interval (%)",
    color = "black",
    fill = "skyblue",
    alpha_shade = .5
  )
}
