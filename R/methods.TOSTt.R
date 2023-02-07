#' Methods for TOSTt objects
#'
#' Methods defined for objects returned from the t_TOST and boot_t_TOST functions.
#'
#' @param x object of class \code{TOSTt}
#' @param digits Number of digits to print for p-values
#' @param type Type of plot to produce. Default is a consonance density plot "cd". Consonance plots (type = "cd") and null distribution plots (type = "tnull") can also be produced. Note: null distribution plots only available for estimates = "raw".
#' @param ci_lines Confidence interval lines for plots. Default is 1-alpha*2 (e.g., alpha = 0.05 is 90\%)
#' @param ci_shades Confidence interval shades when plot type is "cd".
#' @param estimates indicator of what estimates to plot; options include "raw" or "SMD". Default is is both: c("raw","SMD").
#' @param ... further arguments passed through, see description of return value
#'   for details.
#'   \code{\link{TOSTt-methods}}.
#' @return
#' \describe{
#'   \item{\code{print}}{Prints short summary of the tests.}
#'   \item{\code{plot}}{Returns a plot of the effects.}
#'   \item{\code{describe}}{Verbose description of results.}
#' }
#'
#' @name TOSTt-methods


### methods for TOSTt objects

#' @rdname TOSTt-methods
#' @method print TOSTt
#' @export

print.TOSTt <- function(x,
                        digits = 4,
                        ...){
  effsize = x$effsize
  TOST = x$TOST
  TOST$p.value = ifelse(TOST$p.value < 0.001,
                        "< 0.001",
                        round(TOST$p.value, 3))
  effsize$CI = paste0("[",
                      round(effsize$lower.ci,digits),
                      ", ",
                      round(effsize$upper.ci,digits),
                      "]")
  effsize = effsize[c("estimate", "SE", "CI","conf.level")]
  TOST = TOST[c("t","df","p.value")]
  colnames(effsize) = c("Estimate", "SE", "C.I.", "Conf. Level")
  cat("\n")
  cat(strwrap(x$method), sep = "\n")
  #cat(x$hypothesis, "\n", sep = "")
  #cat("Equivalence Bounds (raw):",format(x$eqb$low_eq[1], digits = 3, nsmall = 3, scientific = FALSE)," & ",format(x$eqb$high_eq[1], digits = 3, nsmall = 3, scientific = FALSE), sep="")
  #cat("\n")
  #cat("Alpha Level:", x$alpha, sep="")
  cat("\n")
  cat(x$decision$TOST)
  cat("\n")
  cat(x$decision$ttest)
  cat("\n")
  cat(x$decision$combined)
  cat("\n")
  cat("\n")
  cat("TOST Results \n")
  print(TOST, digits = digits)
  cat("\n")
  cat("Effect Sizes \n")
  print(effsize, digits = digits)

  if(grepl("Log-transformed",x$method, ignore.case=TRUE)){
    cat("\n")
  } else {
  if("boot" %in% names(x)){
    cat("Note: percentile bootstrap method utilized.")
  }else{
    cat("Note: SMD confidence intervals are an approximation. See vignette(\"SMD_calcs\").")
  }
  }
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
                       type = c("cd","c","tnull"),
                       estimates = c("raw","SMD"),
                       ci_lines,
                       ci_shades,
                       ...){
  type = match.arg(type)

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
  }else{
    sets = ci_shades
    if(length(sets) > 4){
      stop("ci_shades cannot have a length greater than 4")
    }
  }

  if(!missing(ci_lines)&& length(ci_lines >1) && type == "cd"){
    stop("ci_lines cannot be greater than 1 if type = \"cd\" ")
  }

  if(missing(ci_lines)){
    ci_levs = 1-x$alpha*2
    c2 = ci_levs
  } else{
    ci_levs = ci_lines
    c2 = ci_lines
  }

  # Get x-axis label
  if(grepl("one",x$method, ignore.case=TRUE)){
    x_label = "Mean"
  } else {
    x_label = "Mean Difference"
  }
  if(grepl("Log-transformed",x$method, ignore.case=TRUE)){
    x_label = "log(Means Ratio)"
  }


  if(type == "c"){
    # type c --------

    if("boot" %in% names(x)){
      warning("Consonance plots from bootstrapped result based on estimates not bootstrap samples")
    }

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

    if(grepl("Log-transformed",x$method, ignore.case=TRUE)){

      d_res = t_res
      d_res[[1]]$lower.limit = exp(d_res[[1]]$lower.limit)
      d_res[[1]]$upper.limit = exp(d_res[[1]]$upper.limit)
      d_res[[1]]$intrvl.width = d_res[[1]]$upper.limit - d_res[[1]]$lower.limit
      d_res[[2]] = exp(d_res[[2]])
      d_plot <-
        gg_curv_t(
          data_list = d_res,
          type = "c",
          levels = ci_levs
        ) +
        geom_vline(xintercept = low_eqt,linetype = "dashed")+
        geom_vline(xintercept = high_eqt, linetype ="dashed")+
        scale_x_continuous(sec.axis = dup_axis(breaks=c(round(low_eqd,2),
                                                        round(high_eqd,2)),
                                               name = "")) +
        facet_grid(~as.character(x_label)) +
        theme_tidybayes() +
        theme(strip.text = element_text(face = "bold",
                                        size = 10),
              axis.title.x = element_blank())

    } else {

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
    }

    if("SMD" %in% estimates && "raw" %in% estimates){
      plts = plot_grid(d_plot,
                       t_plot,
                       ncol = 1)
    }

    if("SMD" %in% estimates && !("raw" %in% estimates)){
      plts = d_plot
    }

    if(!("SMD" %in% estimates) && "raw" %in% estimates){
      plts = t_plot
    }

    return(plts)
  }

  if(type == "cd"){
    # type cd --------

    if(!missing(ci_lines) && length(ci_lines)>1){
      warning("Multiple CI lines provided only first element will be used.")
    }

    if(!("boot" %in% names(x))){

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


    if(grepl("Log-transformed",x$method, ignore.case=TRUE)){
      t_res = t_curv(x)
      d_res = t_res
      d_res[[1]]$lower.limit = exp(d_res[[1]]$lower.limit)
      d_res[[1]]$upper.limit = exp(d_res[[1]]$upper.limit)
      d_res[[1]]$intrvl.width = d_res[[1]]$upper.limit - d_res[[1]]$lower.limit
      d_res[[2]] = exp(d_res[[2]])
    } else {
      d_res = d_curv(x)
    }
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




    }

    if("boot" %in% names(x)){
      df_t = data.frame(val = x$boot$raw,
                        type = x_label,
                        est = c(x$effsize$estimate[1]),
                        low = c(x$eqb$low_eq[1]),
                        high = c(x$eqb$high_eq[1]),
                        alpha = c(x$alpha),
                        stringsAsFactors = FALSE)
      df_d = data.frame(val = x$boot$SMD,
                        type = x$smd$smd_label,
                        est = c(x$effsize$estimate[2]),
                        low = c(x$eqb$low_eq[2]),
                        high = c(x$eqb$high_eq[2]),
                        alpha = c(x$alpha),
                        stringsAsFactors = FALSE)
      t_plot = ggplot(data = df_t,
                      aes(y = 0, x = val)) +
        stat_halfeye(aes(fill = stat(cut_cdf_qi(
          cdf,
          .width = sets
        ))),
        .width = c2,
        slab_color = "black",
        slab_size = .5) +
        scale_fill_viridis_d(option = "D",
                             direction = -1,
                             na.translate = FALSE) +
        labs(x = '', y = '',
             fill = "Confidence Interval") +
        geom_vline(aes(xintercept = low),
                   linetype = "dashed") +
        geom_vline(aes(xintercept = high),
                   linetype = "dashed") +
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

      d_plot = ggplot(data = df_d,
                      aes(y = 0, x = val)) +
        stat_halfeye(aes(fill = stat(cut_cdf_qi(
          cdf,
          .width = sets
        ))),
        .width = c2,
        slab_color = "black",
        slab_size = .5) +
        scale_fill_viridis_d(option = "D",
                             direction = -1,
                             na.translate = FALSE) +
        labs(x = '', y = '',
             fill = "Confidence Interval") +
        geom_vline(aes(xintercept = low),
                   linetype = "dashed") +
        geom_vline(aes(xintercept = high),
                   linetype = "dashed") +
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
          round(low_eqd, 2),
          round(high_eqd, 2)
        )))
    }

    # extract the legend from one of the plots
    legend <- get_legend(t_plot)

    prow <- plot_grid(
      d_plot + theme(legend.position = "none"),
      t_plot + theme(legend.position = "none"),
      ncol = 1
    )

    # add the legend to the row we made earlier. Give it one-third of
    # the width of one plot (via rel_widths).

    if("SMD" %in% estimates && "raw" %in% estimates){
      plts = plot_grid(legend, prow, ncol = 1,
                       rel_heights = c(.1, 1))
    }

    if("SMD" %in% estimates && !("raw" %in% estimates)){
      plts = d_plot
    }

    if(!("SMD" %in% estimates) && "raw" %in% estimates){
      plts = t_plot
    }

    return(plts)
  }

  if(type == "tnull"){
    if("SMD" %in% estimates){
      message("SMD cannot be plotted if type = \"tnull\" ")
    }
    if(!missing(ci_lines) && length(ci_lines)>1){
      warning("Multiple CI lines provided; only first element will be used.")
    }

    if("Equilvalence" %in% x$hypothesis){
      METhyp = TRUE
    } else {
      METhyp = FALSE
    }
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
    points = data.frame(
      x_label = x_label,
      point = x$effsize$estimate[1],
      ci_high = x$effsize$lower.ci[1],
      ci_low = x$effsize$upper.ci[1],
      stringsAsFactors = FALSE
    )
    points_l = data.frame(
      mu = c(x$eqb$low_eq[1]),
      param = c(round(unname(x$TOST$df[1]), 0)),
      sigma = c(x$TOST$SE[1]),
      lambda = c(0),
      stringsAsFactors = FALSE
    )
    points_u = data.frame(
      mu = c(x$eqb$high_eq[1]),
      param = c(round(unname(x$TOST$df[1]), 0)),
      sigma = c(x$TOST$SE[1]),
      lambda = c(0),
      stringsAsFactors = FALSE
    )

    x_l = c(low_eqt - qnorm(1-x$alpha)*points_l$sigma,
            low_eqt + qnorm(1-x$alpha)*points_l$sigma)
    x_u = c(high_eqt - qnorm(1-x$alpha)*points_l$sigma,
                 high_eqt + qnorm(1-x$alpha)*points_l$sigma)

    t_plot = ggplot(data = points,
                    aes_string(y = 0)) +
      stat_dist_slab(data = points_l,
                     aes(fill = stat(x < x_l[1] | x > x_l[2]),
                         dist = dist_student_t(
                           mu = mu,
                           df = param,
                           sigma = sigma,
                           ncp = lambda
                         )),
                     alpha = .5,
                     #  fill = NA,
                     slab_color = "black",
                     slab_size = .5) +
      stat_dist_slab(data = points_u,
                     aes(fill = stat(x < x_u[1] | x > x_u[2]),
                         dist = dist_student_t(
                           mu = mu,
                           df = param,
                           sigma = sigma,
                           ncp = lambda
                         )),

                     alpha = .5,
                     #  fill = NA,
                     slab_color = "black",
                     slab_size = .5) +
      geom_point(data = data.frame(y = -.1,
                                   x = points$point),
                 aes(x = x, y = y),
                 size = 3) +
      annotate("segment",
               x = points$ci_low,
               xend = points$ci_high,
               y = -.1, yend = -.1,
               size = 1.5,
               colour = "black")+
      # set palettes  need true false
      scale_fill_manual(values = c("gray85", "green")) +
      geom_vline(aes(xintercept = low_eqt),
                 linetype = "dashed") +
      geom_vline(aes(xintercept = high_eqt),
                 linetype = "dashed") +
      facet_wrap( ~ x_label) +
      labs(caption = "Note: green indicates rejection region for null equivalence and MET hypotheses")+
      theme_tidybayes() +
      theme(
        legend.position = "none",
        strip.text = element_text(face = "bold", size = 11),
        legend.text = element_text(face = "bold", size = 11),
        legend.title = element_text(face = "bold", size = 11),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
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
    return(t_plot)
  }

}

#' @rdname TOSTt-methods
#' @export

describe <- function(x, ...) {
  UseMethod("describe")
}

#' @rdname TOSTt-methods
#' @method describe TOSTt
#' @export

describe.TOSTt <- function(x,
                           digits = 3,
                           ...){

  text2 = describe_TOST(x = x,
                        digits = digits,
                        ...)

  return(text2)
}


describe_TOST = function(x,
                         digits = 3,
                         ...){
  tosty = x
  htest = as_htest(x)

  type_tost = ifelse(htest$alternative == "equivalence",
                     "equivalence",
                     "minimal effect")
  nhst_null = ifelse(is.null(tosty$call$mu),
                     0,
                     tosty$call$mu)
  alt_nhst = paste0("true ",
                    names(htest$null.value[1]),
                    " is ",
                    "not equal to",
                    " ",
                    nhst_null)

  null_nhst = paste0("true ",
                     names(htest$null.value[1]),
                     " is ",
                     "equal to",
                     " ",
                     nhst_null)
  if(htest$alternative == "equivalence"){

    alt_tost = paste0("true ",
                      names(htest$null.value)[1],
                      " is ",
                      "between",
                      " ",
                      htest$null.value[1], " and ",
                      htest$null.value[2])

    null_tost = paste0("true ",
                       names(htest$null.value)[1],
                       " is ",
                       "more extreme than",
                       " ",
                       htest$null.value[1], " and ",
                       htest$null.value[2])

  }
  if(htest$alternative == "minimal.effect"){

    alt_tost = paste0("true ",
                      names(htest$null.value)[1],
                      " is ",
                      "less than ",
                      htest$null.value[1], " or ",
                      "greater than ",
                      htest$null.value[2])

    null_tost = paste0("true ",
                       names(htest$null.value)[1],
                       " is ",
                       "greater than ",
                       htest$null.value[1], " or ",
                       "less than ",
                       htest$null.value[2])

  }

  method_state = paste0("Using the ", htest$method,
                        ", a null hypothesis significance test (NHST)",
                        ", and a ",type_tost," test, via two one-sided tests (TOST), were performed",
                        " with an alpha-level of ", x$alpha, ".",
                        " These tested the null hypotheses that ", null_nhst,
                        " (NHST), and ", null_tost," (TOST).")



  pTOST = htest$p.value
  sigTOST = ifelse(pTOST < tosty$alpha, TRUE,FALSE)
  pNHST = tosty$TOST$p.value[1]
  sigNHST = ifelse(pNHST < tosty$alpha, TRUE,FALSE)
  stat_name = names(htest$statistic)
  smd_name = row.names(tosty$effsize)[2]


  if(sigTOST){
    if(!is.null(htest$parameter)){
      stat_print = paste0(stat_name, "(",
                          rounder_stat(htest$parameter, digits = digits),
                          ") = ",
                          rounder_stat(htest$statistic, digits = digits))
    } else{
      stat_print = paste0(stat_name, " = ",
                          rounder_stat(htest$statistic, digits = digits))
    }

    sig_text = paste0("The ", type_tost, " test",
                      " was significant, ",
                      stat_print,
                      ", ",
                      printable_pval(pTOST, digits = digits))

    claim_text = paste0("At the desired error rate, it can be stated that the ",
                        alt_tost, ".")
  } else if(sigNHST){
    if(!is.null(htest$parameter)){
      stat_print = paste0(stat_name, "(",
                          rounder_stat(htest$parameter, digits = digits),
                          ") = ",
                          rounder_stat(tosty$TOST[1,1], digits = digits))
    } else{
      stat_print = paste0(stat_name, " = ",
                          rounder_stat(tosty$TOST[1,1], digits = digits))
    }
    sig_text = paste0("The ", type_tost, " test was not significant (",
                      printable_pval(pTOST, digits = digits),").",
                      " The NHST",
                      " was significant, ",
                      stat_print,
                      ", ",
                      printable_pval(pNHST, digits = digits))

    claim_text = paste0("At the desired error rate, it can be stated that the ",
                        alt_nhst,
                        " (i.e., no ",type_tost,
                        ").")
  } else {
    sig_text = paste0("Both the ", type_tost, " test (",
                      printable_pval(pTOST, digits = digits),
                      "), and the NHST (",
                      printable_pval(pNHST, digits = digits),
                      ")",
                      " were not significant")

    claim_text = paste0("Therefore, the results are inconclusive: neither null hypothesis can be rejected.")
  }

  stat_text = paste0(sig_text,
                     " (",
                     names(htest$null.value[1]),
                     " = ",
                     rounder_stat(htest$estimate,
                                  digits = digits),
                     " ",
                     100 * attr(htest$conf.int, "conf.level"),
                     "% C.I.[",
                     rounder_stat(min(htest$conf.int),
                                  digits = digits),
                     ", ",
                     rounder_stat(max(htest$conf.int),
                                  digits = digits),
                     "]",
                     "; ",
                     smd_name,
                     " = ",
                     rounder_stat(tosty$effsize$estimate[2],
                                  digits = digits),
                     " ",
                     100 * attr(htest$conf.int, "conf.level"),
                     "% C.I.[",
                     rounder_stat(tosty$effsize$lower.ci[2],
                                  digits = digits),
                     ", ",
                     rounder_stat(tosty$effsize$upper.ci[2],
                                  digits = digits),
                     "]",
                     ").")

  text2 = paste(method_state, stat_text,
                claim_text,
                sep = " ")

  return(text2)
}
