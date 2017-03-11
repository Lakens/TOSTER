
tostplot <- function(image) {
  ciw <- 90

  data <- image$state

  plot <- ggplot(data=data, aes_string(x=0, y='m')) +
    geom_hline(yintercept=data$low,  colour='#333333') +
    geom_hline(yintercept=data$high, colour='#333333') +
    geom_text(aes(5, data$low,  vjust=-.9, hjust=1), label='Lower bound') +
    geom_text(aes(5, data$high, vjust=1.9, hjust=1), label='Upper bound') +
    geom_errorbar(aes_string(x=0, ymin='cil', ymax='ciu', width=.4), size=.8, colour='#333333') +
    geom_point(aes_string(x=0, y='m'), shape=21, fill='white', size=3, colour='#333333') +
    labs(x='', y='') +
    expand_limits(x=c(-2, 5), y=0) +
    theme(
      text=element_text(size=16, colour='#333333'),
      plot.background=element_rect(fill='transparent', color=NA),
      panel.background=element_rect(fill='#E8E8E8'),
      panel.grid.minor.x = element_blank(),
      axis.text.x=element_blank(),
      axis.title.x=element_blank(),
      axis.ticks.x=element_blank())

  suppressWarnings(print(plot))
}
