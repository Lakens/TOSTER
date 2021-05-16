#' Methods for TOSTt objects
#'
#' Methods defined for objects returned from the agree functions.
#'
#' @param x object of class \code{TOSTt} as returned from the reli_stats function
#' @param digits Number of digits to print for p-values
#' @param type Type of plot to produce. Default is a consonance plot "c" but consonance distribution plot can be produced with "cd".
#' @param ci_lines Confidence interval lines for plots. Default is 1-alpha*2 (e.g., alpha = 0.05 is 90\%)
#' @param ci_shades Confidence interval shades when plot type is "cd".
#' @param ... further arguments passed through, see description of return value
#'   for details.
#'   \code{\link{TOSTt-methods}}.
#' @return
#' \describe{
#'   \item{\code{print}}{Prints short summary of the Limits of Agreement}
#'   \item{\code{plot}}{Returns a plot of the data points used in the reliability analysis}
#' }
#'
#' @name TOSTt-methods


### methods for TOSTt objects

#' @rdname TOSTt-methods
#' @method print TOSTt
#' @export

print.TOSTt <- function(x,
                        digits = getOption("digits"),
                        ...){
  cat("\n")
  cat(strwrap(x$method), sep = "\n")
  cat(x$hypothesis, "\n", sep = "")
  cat("Equivalence Bounds (raw):",format(x$eqb$low_eq[1], digits = 3, nsmall = 3, scientific = FALSE)," & ",format(x$eqb$high_eq[1], digits = 3, nsmall = 3, scientific = FALSE), sep="")
  cat("\n")
  cat("Alpha Level:", x$alpha, sep="")
  cat("\n")
  cat(x$decision$TOST)
  cat("\n")
  cat(x$decision$ttest)
  cat("\n")
  cat("Conclusion: The effect is ",x$decision$combined,".",sep="")
  cat("\n")
  cat("\n")
  cat("TOST Results \n")
  print(x$TOST)
  cat("\n")
  cat("Effect Sizes \n")
  print(x$effsize)
  cat("\n")

}

#' @rdname TOSTt-methods
#' @method plot TOSTt
#' @import ggplot2
#' @import ggdist
#' @import distributional
#' @importFrom cowplot plot_grid get_legend
#' @importFrom stats density dt
#' @importFrom utils head
#' @export

plot.TOSTt <- function(x,
                       type = "cd",
                       ci_lines,
                       ci_shades,
                       ...){

  low_eqd = x$eqb$low_eq[2]
  high_eqd = x$eqb$high_eq[2]

  low_eqt = x$eqb$low_eq[1]
  high_eqt = x$eqb$high_eq[1]

  lenst = c(length(low_eqt),
            length(high_eqt))
  round_t = max(lenst)
  smd_type = x$smd$smd_label


  if(missing(ci_shades)){
    c1 = 1-x$alpha
    c2 = 1-x$alpha*2
    if(c1 > .68 && c2 > .68){
      sets = c(.68,c2,c1,.999)
    } else if(c2 <=.68 && c1 < .999) {
      sets = c(c2,c1,.999)
    } else {
      sets = c(.68,c2,c1)
    }
  }

  if(missing(ci_lines)){
    ci_levs = 1-x$alpha*2
  } else{
    ci_levs = ci_lines
  }
  # Get x-axis label
  if(grepl("one",x$method, ignore.case=TRUE)){
    x_label = "Mean"
  } else {
    x_label = "Mean Difference"
  }


  if(type == "c"){

    d_res = d_curv(x)

    d_plot <-
      gg_curv_t(
        data_list = d_res,
        type = "c",
        levels = ci_levs
      ) +
      geom_vline(xintercept = low_eqd,linetype = "dashed")+
      geom_vline(xintercept = high_eqd, linetype ="dashed")+
      scale_x_continuous(sec.axis = dup_axis(breaks=c(round(low_eqd,2),
                                                      round(high_eqd,2)),
                                             name = "")) +
      facet_grid(~as.character(x$smd$smd_label)) +
      theme_tidybayes() +
      theme(strip.text = element_text(face = "bold",
                                      size = 10),
            axis.title.x = element_blank())

    t_res = t_curv(x)

    t_plot <-
      gg_curv_t(
        data_list = t_res,
        type = "c",
        levels = ci_levs
      ) +
      geom_vline(xintercept = low_eqt,linetype = "dashed")+
      geom_vline(xintercept = high_eqt, linetype ="dashed")+
      scale_x_continuous(sec.axis = dup_axis(breaks=c(round(low_eqt,round_t),
                                                      round(high_eqt,round_t)),
                                             name = "")) +
      facet_grid(~as.character(x_label)) +
      theme_tidybayes() +
      theme(strip.text = element_text(face = "bold",
                                      size = 10),
            axis.title.x = element_blank())

    plts = plot_grid(d_plot,
                     t_plot,
                     ncol = 1)
    return(plts)
  }

  if(type == "cd"){

    if(!missing(ci_lines) && length(ci_lines)>1){
      warning("Multiple CI lines provided only first element will be used.")
    }


    d_res = d_curv(x)
    d_plot <- gg_curv_t(
      data_list = d_res,
      type = "cd",
      levels = sets,
      xaxis = ""
    ) +
      geom_point(data = data.frame(y = 0,
                                   x = x$effsize$estimate[2]),
                 aes(x = x, y = y),
                 size = 3) +
      annotate("segment",
               x = x$effsize$lower.ci[2],
               xend = x$effsize$upper.ci[2],
               y = 0, yend = 0,
               size = 1.5,
               colour = "black")+
      geom_vline(aes(xintercept = low_eqd),
                 linetype="dashed") +
      geom_vline(aes(xintercept = high_eqd),
                 linetype="dashed") +
      scale_x_continuous(sec.axis = dup_axis(breaks=c(round(low_eqd,2),
                                                      round(high_eqd,2)))) +
      facet_grid(~as.character(smd_type)) +
      labs(y = "")+
      theme_tidybayes() +
      theme(
        legend.position = "top",
        strip.text = element_text(face = "bold", size = 11),
        legend.text = element_text(face = "bold", size = 11),
        legend.title = element_text(face = "bold", size = 11),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(face = "bold", size = 11),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        legend.background = element_rect(fill = "transparent",colour = NA)
      )


    points = data.frame(
      type = x_label,
      mu = c(x$effsize$estimate[1]),
      param = c(round(unname(x$TOST$df[1]), 0)),
      sigma = c(x$TOST$SE[1]),
      lambda = c(0),
      est = c(x$effsize$estimate[1]),
      low = c(x$eqb$low_eq[1]),
      high = c(x$eqb$high_eq[1]),
      alpha = c(x$alpha),
      stringsAsFactors = FALSE
    )

    #sets = ci.cuts
    t_plot = ggplot(data = points,
                    aes_string(y = 0)) +
      stat_dist_halfeye(aes(
        dist = dist_student_t(
          mu = mu,
          df = param,
          sigma = sigma,
          ncp = lambda
        ),
        fill = stat(cut_cdf_qi(p = cdf,
                               .width = sets))
      ),
      .width = c2,
      slab_color = "black",
      slab_size = .5) +
      #scale_fill_brewer(direction = -1,
      #                  na.translate = FALSE) +
      scale_fill_viridis_d(option = "D",
                           direction = -1,
                           na.translate = FALSE) +
      labs(x = '', y = '',
           fill = "Confidence Interval") +
      geom_vline(aes(xintercept = low),
                 linetype = "dashed") +
      geom_vline(aes(xintercept = high),
                 linetype = "dashed") +
      #geom_text(aes(y=1.5, x=low,
      #              vjust=-.9, hjust=1),
      #          angle = 90,
      #          label='Lower Bound') +
      #geom_text(aes(y=1.5, x=high, vjust=1.5, hjust=1),
      #          angle = 90,
      #         label='Upper Bound') +
      facet_wrap( ~ type) +
      theme_tidybayes() +
      theme(
        legend.position = "top",
        strip.text = element_text(face = "bold", size = 11),
        legend.text = element_text(face = "bold", size = 11),
        legend.title = element_text(face = "bold", size = 11),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(face = "bold", size = 11),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        legend.background = element_rect(fill = "transparent",colour = NA)
      ) +
      scale_x_continuous(sec.axis = dup_axis(breaks = c(
        round(low_eqt, round_t),
        round(high_eqt, round_t)
      )))

    # extract the legend from one of the plots
    legend <- get_legend(t_plot)

    prow <- plot_grid(
      d_plot + theme(legend.position = "none"),
      t_plot + theme(legend.position = "none"),
      ncol = 1
    )

    # add the legend to the row we made earlier. Give it one-third of
    # the width of one plot (via rel_widths).


    plts = plot_grid(legend, prow, ncol = 1,
                     rel_heights = c(.1, 1))

    return(plts)
  }

}
