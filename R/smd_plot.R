#' Function to produce plots of the distribution of standardized mean difference
#' @param lambda The non-centrality parameter for the standardized mean difference
#' @param df Degrees of freedom for the standardized mean difference
#' @param SE The standard error of the for the standardized mean difference (square root of the variance of the SMD)
#' @param smd_label The label that is desired for the plot
#' @param ci_shade At which confidence intervals to change the color shade of the distribution
#' @param ci_lines The confidence intervals to display on the line at the bottom of the plots
#' @details
#' This function was created so that users could create plots from their own SMD calculations.
#' The difficulty is that specific information must be past onto this function.
#' The non-centrality parameter (lambda), degrees of freedom (df), and standard error (SE) calculations can be found in the vignettes of this package.
#' @return Returns plot of the distribution of the standardized mean difference.
#' @import ggdist
#' @import distributional
#' @export

plot_smd <- function(lambda,
                     df,
                     SE,
                     smd_label,
                     ci_shade = c(0.5,0.90,0.95,0.999),
                     ci_lines = 0.90){

  if (!missing(smd_label)){
    points = data.frame(mu = 0,
                        param = df,
                        sigma = SE,
                        ncp = lambda,
                        smd_label = smd_label)

    sets <- ci_shade

    p1 = ggplot(data = points,
                aes_string(y = 0)) +
      stat_dist_halfeye(aes(
        dist = dist_student_t(
          mu = mu,
          df = param,
          sigma = sigma,
          ncp = lambda
        ),
        fill = stat(cut_cdf_qi(p=cdf,
                               .width = sets))
      ),
      .width = c(ci_lines)) +
      scale_fill_brewer(direction = -1,
                        na.translate = FALSE) +
      labs(x = '', y = '',
           fill = "Confidence Interval") +
      theme_tidybayes() +
      theme(legend.position="top",
            strip.text = element_text(face = "bold", size = 10),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank()) +
      facet_wrap(~smd_label,
                 ncol = 2,
                 scales = "free")


  } else {
    points = data.frame(mu = 0,
                        param = df,
                        sigma = SE,
                        ncp = lambda)

    sets <- ci_shade

    p1 = ggplot(data = points,
                aes_string(y = 0)) +
      stat_dist_halfeye(aes(
        dist = dist_student_t(
          mu = mu,
          df = param,
          sigma = sigma,
          ncp = lambda
        ),
        fill = stat(cut_cdf_qi(p=cdf,
                               .width = sets))
      ),
      .width = c(ci_lines)) +
      scale_fill_brewer(direction = -1,
                        na.translate = FALSE) +
      labs(x = '', y = '',
           fill = "Confidence Interval") +
      theme_tidybayes() +
      theme(legend.position="top",
            strip.text = element_text(face = "bold", size = 10),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
  }

  return(p1)

}


