#' Methods for TOSTt objects
#'
#' Methods defined for objects returned from the t_TOST and boot_t_TOST functions.
#'
#' @param x object of class `TOSTt`.
#' @param digits Number of digits to print for p-values
#' @param type Type of plot to produce. Default is "simple" which shows point estimates with confidence intervals. Other options include consonance plots ("c"), consonance density plots ("cd"), and null distribution plots ("tnull"). Note: null distribution plots only available for estimates = "raw".
#' @param ci_lines Confidence interval lines for plots. Default is 1-alpha*2 (e.g., alpha = 0.05 is 90%)
#' @param ci_shades Confidence interval shades when plot type is "cd".
#' @param estimates indicator of what estimates to plot; options include "raw" or "SMD". Default is is both: c("raw","SMD").
#' @param layout Layout for displaying multiple estimates. Options are "stacked" (default, separate plots stacked vertically) or "combined" (single faceted plot with shared legend). Only applies when both "raw" and "SMD" are in estimates.
#' @param ... further arguments passed through, see description of return value for details..
#'
#' @return
#'   - print: Prints short summary of the tests.
#'   - plot: Returns a plot of the effects.
#'   - describe: Verbose description of results.
#'
#'
#' @examples
#' # example code
#' # Print
#'
#' res1 = t_TOST(mpg ~ am, data = mtcars, eqb = 3)
#'
#' res1
#' # Print with more digits
#' print(res1, digits = 6)
#'
#' # Plot with density plot - only raw values (SLOW)
#' #plot(res1, type = "cd", estimates = "raw")
#' # Plot with consonance - only raw values (SLOW)
#' #plot(res1, type = "c", estimates = "raw")
#' # Plot null distribution - only raw values
#' #plot(res1, type = "tnull", estimates = "raw")
#'
#' # Get description of the results
#' describe(res1)
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
    if(x$call$boot_ci == "perc"){
      cat("Note: percentile bootstrap ci method utilized.")
    }
    if(x$call$boot_ci == "basic"){
      cat("Note: basic bootstrap ci method utilized.")
    }

    if(x$call$boot_ci == "stud"){
      cat("Note: studentized bootstrap ci method utilized.")
    }

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
                       type = c("simple","cd","c","tnull"),
                       estimates = c("raw","SMD"),
                       ci_lines,
                       ci_shades,
                       layout = c("stacked", "combined"),
                       ...){
  type = match.arg(type)
  layout = match.arg(layout)

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

  if(type == "simple"){
    # type simple --------

    ci_print = x$effsize$conf.level[1]

    # Build subtitle with decision text and equivalence bounds
    eqb_text = paste0("Equivalence bounds: [",
                      round(low_eqt, 3), ", ",
                      round(high_eqt, 3), "] (raw)")
    subtitle_text = paste0(x$decision$TOST, "\n",
                           x$decision$ttest, "\n",
                           eqb_text)

    # Get estimates for mean ----
    df_t = x$effsize[1,]

    # Get estimates for SMD ----
    df_d = x$effsize[2,]

    # Check if we should use combined layout
    both_estimates = "SMD" %in% estimates && "raw" %in% estimates

    if(both_estimates && layout == "combined"){
      # Combined faceted plot ----
      # Create combined data frame with scale indicator
      df_combined = rbind(
        data.frame(
          estimate = df_t$estimate,
          lower.ci = df_t$lower.ci,
          upper.ci = df_t$upper.ci,
          type = paste0(x_label, " (raw)"),
          low_eq = low_eqt,
          high_eq = high_eqt,
          stringsAsFactors = FALSE
        ),
        data.frame(
          estimate = df_d$estimate,
          lower.ci = df_d$lower.ci,
          upper.ci = df_d$upper.ci,
          type = paste0(x$smd$smd_label, " (standardized)"),
          low_eq = low_eqd,
          high_eq = high_eqd,
          stringsAsFactors = FALSE
        )
      )
      # Preserve order: raw first, then SMD
      df_combined$type = factor(df_combined$type,
                                levels = c(paste0(x_label, " (raw)"),
                                           paste0(x$smd$smd_label, " (standardized)")))

      plts <- ggplot(df_combined,
                     aes(x = estimate,
                         y = 1,
                         xmin = lower.ci,
                         xmax = upper.ci)) +
        geom_pointrange() +
        geom_vline(aes(xintercept = low_eq), linetype = "dashed") +
        geom_vline(aes(xintercept = high_eq), linetype = "dashed") +
        facet_wrap(~type, scales = "free_x") +
        theme_tidybayes() +
        labs(title = subtitle_text,
             caption = paste0(ci_print*100, "% Confidence Interval")) +
        theme(strip.text = element_text(face = "bold", size = 10),
              plot.title = element_text(size = 10),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank())

      return(plts)
    }

    # Stacked layout (default) ----
    # Build facet labels with scale indicator
    raw_facet_label = paste0(x_label, " (raw)")
    smd_facet_label = paste0(x$smd$smd_label, " (standardized)")

    # Raw plot (now shown on top with subtitle)
    t_plot <-
      ggplot(df_t,
             aes(x=estimate,
                 y = 1,
                 xmin=lower.ci,
                 xmax=upper.ci)) +
      geom_pointrange() +
      geom_vline(xintercept = low_eqt,linetype = "dashed")+
      geom_vline(xintercept = high_eqt, linetype ="dashed")+
      scale_x_continuous(sec.axis = dup_axis(breaks=c(round(low_eqt,round_t),
                                                      round(high_eqt,round_t)),
                                             name = "")) +
      facet_grid(~as.character(raw_facet_label)) +
      theme_tidybayes() +
      labs(caption = paste0(ci_print*100,"% Confidence Interval"),
           title = subtitle_text) +
      theme(strip.text = element_text(face = "bold",
                                      size = 10),
            plot.title = element_text(size = 10),
            axis.title.x = element_blank(),
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank())

    # SMD plot (now shown on bottom)
    d_plot <-
      ggplot(df_d,
             aes(x=estimate,
                 y = 1,
                 xmin=lower.ci,
                 xmax=upper.ci)) +
      geom_pointrange() +
      geom_vline(xintercept = low_eqd,linetype = "dashed")+
      geom_vline(xintercept = high_eqd, linetype ="dashed")+
      scale_x_continuous(sec.axis = dup_axis(breaks=c(round(low_eqd,round_t),
                                                      round(high_eqd,round_t)),
                                             name = "")) +
      facet_grid(~as.character(smd_facet_label)) +
      labs(caption = paste0(ci_print*100,"% Confidence Interval")) +
      theme_tidybayes() +
      theme(strip.text = element_text(face = "bold",
                                      size = 10),
            axis.title.x = element_blank(),
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank())

    # Stack plots: raw on top, SMD on bottom
    if(both_estimates){
      plts = plot_grid(t_plot,
                       d_plot,
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
      #legboth = get_legend(d_plot)
      plts = plot_grid(d_plot,
                       t_plot,
                       #rel_heights = c(.1,.4,.4),
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
        fill = after_stat(cut_cdf_qi(p = cdf,
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
        legend.position = "bottom",
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
          legend.position = "bottom",
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
        stat_halfeye(aes(fill = after_stat(cut_cdf_qi(
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
          legend.position = "bottom",
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
        stat_halfeye(aes(fill =after_stat(cut_cdf_qi(
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
          legend.position = "bottom",
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

  if(type == "tnull"){
    # tnull ------
    if("SMD" %in% estimates){
      message("SMD cannot be plotted if type = \"tnull\" ")
    }
    if(!missing(ci_lines) && length(ci_lines)>1){
      warning("Multiple CI lines provided; only first element will be used.")
    }

    # Determine if this is equivalence or MET hypothesis
    is_equivalence = grepl("Equivalence", x$hypothesis, ignore.case = TRUE)

    # Common parameters
    se = x$TOST$SE[1]
    df_t = round(unname(x$TOST$df[1]), 0)

    # Data for point estimate and CI
    points = data.frame(
      x_label = x_label,
      point = x$effsize$estimate[1],
      ci_high = x$effsize$lower.ci[1],
      ci_low = x$effsize$upper.ci[1],
      stringsAsFactors = FALSE
    )

    points_l = data.frame(
      mu = c(low_eqt),
      param = c(df_t),
      sigma = c(se),
      lambda = c(0),
      stringsAsFactors = FALSE
    )
    points_u = data.frame(
      mu = c(high_eqt),
      param = c(df_t),
      sigma = c(se),
      lambda = c(0),
      stringsAsFactors = FALSE
    )

    # Calculate one-sided critical values
    # For equivalence: lower bound uses right tail (greater), upper bound uses left tail (less)
    # For MET: lower bound uses left tail (less), upper bound uses right tail (greater)
    if (is_equivalence) {
      crit_l_right = low_eqt + qnorm(1 - x$alpha) * se
      crit_u_left = high_eqt - qnorm(1 - x$alpha) * se
    } else {
      crit_l_left = low_eqt - qnorm(1 - x$alpha) * se
      crit_u_right = high_eqt + qnorm(1 - x$alpha) * se
    }

    # Build plot with one-sided rejection regions
    t_plot = ggplot(data = points,
                    aes_string(y = 0))

    if (is_equivalence) {
      t_plot = t_plot +
        stat_dist_slab(data = points_l,
                       aes(fill = after_stat(x > crit_l_right),
                           dist = dist_student_t(
                             mu = mu,
                             df = param,
                             sigma = sigma,
                             ncp = lambda
                           )),
                       alpha = .5,
                       slab_color = "black",
                       slab_size = .5) +
        stat_dist_slab(data = points_u,
                       aes(fill = after_stat(x < crit_u_left),
                           dist = dist_student_t(
                             mu = mu,
                             df = param,
                             sigma = sigma,
                             ncp = lambda
                           )),
                       alpha = .5,
                       slab_color = "black",
                       slab_size = .5)
    } else {
      t_plot = t_plot +
        stat_dist_slab(data = points_l,
                       aes(fill = after_stat(x < crit_l_left),
                           dist = dist_student_t(
                             mu = mu,
                             df = param,
                             sigma = sigma,
                             ncp = lambda
                           )),
                       alpha = .5,
                       slab_color = "black",
                       slab_size = .5) +
        stat_dist_slab(data = points_u,
                       aes(fill = after_stat(x > crit_u_right),
                           dist = dist_student_t(
                             mu = mu,
                             df = param,
                             sigma = sigma,
                             ncp = lambda
                           )),
                       alpha = .5,
                       slab_color = "black",
                       slab_size = .5)
    }

    t_plot = t_plot +
      geom_point(data = data.frame(y = -.1,
                                   x = points$point),
                 aes(x = x, y = y),
                 size = 3) +
      annotate("segment",
               x = points$ci_low,
               xend = points$ci_high,
               y = -.1, yend = -.1,
               size = 1.5,
               colour = "black") +
      scale_fill_manual(values = c("gray85", "green")) +
      geom_vline(aes(xintercept = low_eqt),
                 linetype = "dashed") +
      geom_vline(aes(xintercept = high_eqt),
                 linetype = "dashed") +
      facet_wrap(~ x_label) +
      labs(caption = paste0("Note: green indicates one-sided rejection region (",
                            ifelse(is_equivalence, "equivalence", "minimal effect"),
                            " test)")) +
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
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        legend.background = element_rect(fill = "transparent", colour = NA)
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
