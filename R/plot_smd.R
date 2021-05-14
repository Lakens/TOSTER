#' plot_smd
#' @description Function to produce plots of the distribution of the standardized mean difference
#' @param d Estimate of the standardized mean difference
#' @param df degrees of freedom for the standardized mean difference
#' @param lambda The non-centrality parameter for the standardized mean difference
#' @param smd_label Label for the x-axis indicating the SMD measure
#' @param type Choose whether to plot a "consonance" function ("c"), consonance density ("cd"), or both (c("c","cd"); defualt option).
#' @param levels Numeric vector of confidence levels to display
#' @details
#' This function was created so that users could create plots from their own SMD calculations and were inspired by the concurve R package.
#' The difficulty is that specific information must be past onto this function.
#' The calculations for the standardized mean difference can be found in the vignettes of this package.
#' These types of plots are discussed by Schweder T, Hjort NL. (2016, ISBN:9781316445051) and Rafi Z, Greenland S. (2020) <doi:10.1186/s12874-020-01105-9>.
#' @return Returns plot of the distribution of the standardized mean difference.
#' @export

plot_smd <- function(d,
                     df,
                     lambda,
                     smd_label,
                     type = c("c","cd"),
                     levels = c(.5,.9,.95,.999)){

  if (missing(smd_label)){

    dat = d_curv_raw(d = d,
                      df = df,
                      lambda=lambda)

    resplot = gg_curv_t(dat,
                        type = type,
                        levels = levels,
                        position = "pyramid",
                        xaxis = "SMD",
                        yaxis1 = expression(paste("two-tailed ", italic(p),
                                                  "-value")),
                        yaxis2 = "Confidence Interval (%)",
                        color = "black",
                        fill = "skyblue",
                        alpha_shade = .5)

  } else {
    dat = d_curv_raw(d = d,
                      df = df,
                      lambda=lambda)

    resplot = gg_curv_t(dat,
                        type = type,
                        levels = levels,
                        position = "pyramid",
                        xaxis = as.character(smd_label),
                        yaxis1 = expression(paste("two-tailed ", italic(p),
                                                  "-value")),
                        yaxis2 = "Confidence Interval (%)",
                        color = "black",
                        fill = "skyblue",
                        alpha_shade = .5)

  }

  return(resplot)

}

