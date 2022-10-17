plot_tost_jam <- function(x,
                          type = "cd",
                          estimates = c("raw", "SMD"),
                          ci_lines,
                          ci_shades,
                          ggtheme){

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


  if(type == "c"){

    if("boot" %in% names(x)){
      warning("Consonance plots from bootstrapped result based on estimates not bootstrap samples")
    }

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
      ggtheme
    #theme_tidybayes() +
    #theme(strip.text = element_text(face = "bold",
    #                                size = 10),
    #      axis.title.x = element_blank())

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
      ggtheme
    #theme_tidybayes() +
    #theme(strip.text = element_text(face = "bold",
    #                                size = 10),
    #      axis.title.x = element_blank())

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

    if(!missing(ci_lines) && length(ci_lines)>1){
      warning("Multiple CI lines provided only first element will be used.")
    }

    if(!("boot" %in% names(x))){


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
        ggtheme +
        theme(
          legend.position = "top",
          #strip.text = element_text(face = "bold", size = 11),
          #legend.text = element_text(face = "bold", size = 11),
          #legend.title = element_text(face = "bold", size = 11),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          #axis.text.x = element_text(face = "bold", size = 11),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          #panel.background = element_rect(fill = "transparent",colour = NA),
          #plot.background = element_rect(fill = "transparent",colour = NA),
          #legend.background = element_rect(fill = "transparent",colour = NA)
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
        na.rm = TRUE,
        slab_color = "black",
        slab_size = .5) +
        #scale_fill_brewer(direction = -1,
        #                  na.translate = FALSE) +
        #scale_fill_viridis_d(option = "D",
        #                     direction = -1,
        #                     na.translate = FALSE) +
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
        ggtheme +
        scale_x_continuous(sec.axis = dup_axis(breaks = c(
          round(low_eqt, round_t),
          round(high_eqt, round_t)
        ))) +
        theme(
          legend.position = "top",
          #strip.text = element_text(face = "bold", size = 11),
          #legend.text = element_text(face = "bold", size = 11),
          #legend.title = element_text(face = "bold", size = 11),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          #axis.text.x = element_text(face = "bold", size = 11),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          #panel.background = element_rect(fill = "transparent",colour = NA),
          #plot.background = element_rect(fill = "transparent",colour = NA),
          #legend.background = element_rect(fill = "transparent",colour = NA)
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
        ggtheme +
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
        ggtheme +
        scale_x_continuous(sec.axis = dup_axis(breaks = c(
          round(low_eqd, 2),
          round(high_eqd, 2)
        )))
    }

    # extract the legend from one of the plots
    legend <- get_legend(d_plot)

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
      warning("Multiple CI lines provided only first element will be used.")
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
      ggtheme +
      scale_x_continuous(sec.axis = dup_axis(breaks = c(
        round(low_eqt, round_t),
        round(high_eqt, round_t)
      )))
    return(t_plot)
  }

}
