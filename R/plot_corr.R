#' Function to produce plots of the distribution of Pearson's correlation coefficient
#' @param r The observed correlation coefficient.
#' @param n Total number of observations (sample size).
#' @param type Choose whether to plot a "consonance" function ("c"), consonance density ("cd"), or both (c("c","cd"); defualt option).
#' @param levels Numeric vector of confidence levels to display
#' @details
#' This function was created so that users could create consonance plots of Pearson's correlation coefficient.
#' These types of plots are discussed by Schweder T, Hjort NL. (2016, ISBN:9781316445051) and Rafi Z, Greenland S. (2020) <doi:10.1186/s12874-020-01105-9>.
#' @return Returns plot of the distribution of the correlation coefficient.
#' @export

plot_corr <- function(r,
                      n,
                      type = c("c","cd"),
                      levels = c(.68,.9,.95,.999)){
  dat = corr_curv(r = r,
                 n = n,
                 steps = 5000)

  resplot = gg_curv_t(
    dat,
    type = type,
    levels = levels,
    position = "pyramid",
    xaxis = "Pearson's Correlation Coefficent (r)",
    yaxis1 = expression(paste(italic(p),
                              "-value")),
    yaxis2 = "Confidence Interval (%)",
    color = "black",
    fill = "skyblue",
    alpha_shade = .5
  )
}
