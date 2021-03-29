
tostplot <- function(image, ggtheme = NULL, theme = NULL) {
  ciw <- 90

  # something wrong with what is being passed in the points dataframe
  data <- image$state


  data2 <- data.frame(
    m=0,
    degree_f = 17,
    SE=.5,
    cil=-5,
    ciu=5,
    low=-.5,
    high=.5,
    stringsAsFactors=FALSE)

  plot <- ggplot(data=data2, aes_string(x=0, y='m')) +
    geom_hline(yintercept=data$low,  colour=theme$color[1]) +
    geom_hline(yintercept=data$high, colour=theme$color[1]) +
    geom_text(aes(5, data$low,  vjust=-.9, hjust=1), label='Lower bound', colour=theme$color[1]) +
    geom_text(aes(5, data$high, vjust=1.9, hjust=1), label='Upper bound', colour=theme$color[1]) +
    geom_errorbar(aes_string(x=0, ymin='cil', ymax='ciu', width=.4), size=.8, colour=theme$color[1]) +
    geom_point(aes_string(x=0, y='m'), shape=21, fill=theme$fill[1], size=3, colour=theme$color[1]) +
    labs(x='', y='') +
    expand_limits(x=c(-2, 5), y=0) #+
  #  ggtheme +
  #  theme(
  #    axis.text.x=element_blank(),
  #    axis.title.x=element_blank(),
  #    axis.ticks.x=element_blank(),
  #    axis.title.y=element_blank(),)

  plot2 <- ggplot(data=data2,
                 aes_string(y=0)) +
  #geom_hline(yintercept=data$low,  colour=theme$color[1]) +
  #geom_hline(yintercept=data$high, colour=theme$color[1]) +
  #geom_text(aes(5, data$low,  vjust=-.9, hjust=1),
  #          label='Lower bound', colour=theme$color[1]) +
  #geom_text(aes(5, data$high, vjust=1.9, hjust=1),
  #          label='Upper bound', colour=theme$color[1]) +
  #geom_errorbar(aes_string(x=0, ymin='cil', ymax='ciu', width=.4),
  #              size=.8, colour=theme$color[1]) +
  #geom_point(aes_string(x=0, y='m'), shape=21, fill=theme$fill[1],
  #           size=3, colour=theme$color[1]) +
    ggdist::stat_dist_halfeye(aes(
      dist = distributional::dist_student_t(
        df = data$degree_f,
        mu = data$m,
        sigma = data$SE
      )
    ))

  suppressWarnings(print(plot))
}
