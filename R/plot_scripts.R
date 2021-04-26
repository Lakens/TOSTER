# Plot scripts

t_curv = function (TOST_res,
                   steps = 10000) {
  intrvls <- (0:steps)/steps

  # Confidence interval of the SMD from Goulet-Pelletier & Cousineau
  results <- suppressWarnings({lapply(intrvls, FUN = function(i) t_CI(TOST_res, 1-i))})

  df <- data.frame(do.call(rbind, results))
  intrvl.limit <- c("lower.limit", "upper.limit")
  colnames(df) <- intrvl.limit
  df$intrvl.width <- (abs((df$upper.limit) - (df$lower.limit)))
  df$intrvl.level <- intrvls
  df$cdf <- (abs(df$intrvl.level/2)) + 0.5
  df$pvalue <- 1 - intrvls
  df$svalue <- -log2(df$pvalue)
  df <- head(df, -1)
  class(df) <- c("data.frame", "concurve")
  densdf <- data.frame(c(df$lower.limit, df$upper.limit))
  colnames(densdf) <- "x"
  densdf <- head(densdf, -1)
  class(densdf) <- c("data.frame", "concurve")

  return(list(df, densdf))

}

d_curv = function (TOST_res,
                   steps = 10000) {
  intrvls <- (0:steps)/steps

  # Confidence interval of the SMD from Goulet-Pelletier & Cousineau
  results <- suppressWarnings({lapply(intrvls, FUN = function(i) d_CI(d = TOST_res$smd$d,
                                                                      df = TOST_res$smd$d_df,
                                                                      lambda = TOST_res$smd$d_lambda,
                                                                      1-i))})

  df <- data.frame(do.call(rbind, results))
  intrvl.limit <- c("lower.limit", "upper.limit")
  colnames(df) <- intrvl.limit
  df$intrvl.width <- (abs((df$upper.limit) - (df$lower.limit)))
  df$intrvl.level <- intrvls
  df$cdf <- (abs(df$intrvl.level/2)) + 0.5
  df$pvalue <- 1 - intrvls
  df$svalue <- -log2(df$pvalue)
  df <- head(df, -1)
  class(df) <- c("data.frame", "concurve")
  densdf <- data.frame(c(df$lower.limit, df$upper.limit))
  colnames(densdf) <- "x"
  densdf <- head(densdf, -1)
  class(densdf) <- c("data.frame", "concurve")

  return(list(df, densdf))

}

plot_smd_curv = function(d,
                         df,
                         lambda,
                         ntilde,
                         r12,
                         smd_label = "Cohen's d",
                         type = "two",
                         ci_shades = c(.5, .90, .95, .99),
                         ci_line = .90){

  # lambda = delta*sqrt(ntilde/2)
  if (type == "two") {
    mult = sqrt(2 / ntilde)
  }
  if (type == "paired") {
    mult = sqrt((2 * (1 - r12)) / ntilde)
  }
  if (type == "one") {
    mult = sqrt(1 / ntilde)
  }

  ci_shade1 = sort(ci_shades, decreasing = TRUE)
  x_dt <- suppressWarnings({seq(suppressWarnings({qt(.0001, df, lambda)}),
                                suppressWarnings({qt(.9999, df, lambda)}),
                               by = 0.001)})
  y_dt <- suppressWarnings({dt(x_dt, df = df, ncp = lambda)})
  x <- mult * x_dt # mult = sqrt(1/n1+1/n2) for Cohen's d
  y <- y_dt / mult # mult = sqrt(1/n1+1/n2) Cohen's d

  x.dens  <- density(x)
  df.dens <- data.frame(x=x.dens$x, y=x.dens$y)

  ci_linerange = d_CI(d = d,
                      df = df,
                      lambda = lambda,
                      alpha = 1-ci_line)
  ci_shaderange1 = d_CI(d = d,
                        df = df,
                        lambda = lambda,
                        alpha = 1-ci_shade1[1])

  df_xy = data.frame(x = x,
                     y = y,
                     smd_label = smd_label)
  p1 = ggplot(df_xy,
              aes(x = x, y = y)) +
    #geom_line(color="white") +
    geom_area(subset(df.dens, x >= ci_shaderange1[1] & x <= ci_shaderange1[2]),
              aes(x = x, y = y, fill = as.character(ci_shade1[1]))) +
    labs(x = '', y = '',
         fill = "Confidence Interval")

  if(length(ci_shade1 > 1)){
    ci_shaderange2 = d_CI(d = d,
                          df = df,
                          lambda = lambda,
                          alpha = 1-ci_shade1[2])
    p2 = p1 +
      geom_area(data = subset(df.dens, x >= ci_shaderange2[1] & x <= ci_shaderange2[2]),
                aes(x = x, y = y, fill = as.character(ci_shade1[2])))
  } else {
    p2 = p1
  }

  if(length(ci_shade1 > 2)){
    ci_shaderange3 = d_CI(d = d,
                          df = df,
                          lambda = lambda,
                          alpha = 1-ci_shade1[3])
    p2 = p2 +
      geom_area(data = subset(df.dens, x >= ci_shaderange3[1] & x <= ci_shaderange3[2]),
                aes(x = x, y = y, fill = as.character(ci_shade1[3])))
  }

  if(length(ci_shade1 > 3)){
    ci_shaderange4 = d_CI(d = d,
                          df = df,
                          lambda = lambda,
                          alpha = 1-ci_shade1[4])
    p2 = p2 +
      geom_area(data = subset(df.dens, x >= ci_shaderange4[1] & x <= ci_shaderange4[2]),
                aes(x = x, y = y, fill = as.character(ci_shade1[4])))
  }

  p2 = p2 +
    annotate("segment",
             x = ci_linerange[1],
             xend = ci_linerange[2],
             y = 0, yend = 0,
             size = 1.5,
             colour = "black") +
    geom_point(data = data.frame(y = 0,
                                 x = d),
               aes(x = x, y = y),
               size = 3) +
    scale_fill_brewer(direction = -1,
                      na.translate = FALSE) +
    theme_tidybayes() +
    facet_grid(~smd_label) +
    labs(x = "") +
    theme(legend.position="top",
          strip.text = element_text(face="bold", size=10),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  return(p2)
}

plot_smd_cdf = function(cdf_dat,
                        d,
                        df,
                        lambda ,
                        ci_shades = c(.5, .90, .95, .99),
                        ci_line = .90){
  ci_shade1 = sort(ci_shades, decreasing = TRUE)

  ci_linerange = d_CI(d = d,
                      df = df,
                      lambda = lambda,
                      alpha = 1-ci_line)
  ci_shaderange1 = d_CI(d = d,
                        df = df,
                        lambda = lambda,
                        alpha = 1-ci_shade1[1])

  cdf_dat2 = cdf_dat$x

  x.dens  <- density(cdf_dat2)
  df.dens <- data.frame(x=x.dens$x, y=x.dens$y)

  p1 = ggplot(data = cdf_dat) +
    geom_density(aes(x = x, y = ..density..),
                 color = "white") +
    geom_area(data = subset(df.dens, x >= ci_shaderange1[1] & x <= ci_shaderange1[2]),
              aes(x = x, y = y, fill = as.character(ci_shade1[1]))) +
    scale_fill_brewer(direction = -1,
                      na.translate = FALSE) +
    labs(x = '', y = '',
         fill = "Confidence Interval")

  if(length(ci_shade1 > 1)){
    ci_shaderange2 = d_CI(d = d,
                          df = df,
                          lambda = lambda,
                          alpha = 1-ci_shade1[2])
    p2 = p1 +
      geom_area(data = subset(df.dens, x >= ci_shaderange2[1] & x <= ci_shaderange2[2]),
                aes(x = x, y = y, fill = as.character(ci_shade1[2])))
  } else {
    p2 = p1
  }

  if(length(ci_shade1 > 2)){
    ci_shaderange3 = d_CI(d = d,
                          df = df,
                          lambda = lambda,
                          alpha = 1-ci_shade1[3])
    p2 = p2 +
      geom_area(data = subset(df.dens, x >= ci_shaderange3[1] & x <= ci_shaderange3[2]),
                aes(x = x, y = y, fill = as.character(ci_shade1[3])))
  }

  if(length(ci_shade1 > 3)){
    ci_shaderange4 = d_CI(d = d,
                          df = df,
                          lambda = lambda,
                          alpha = 1-ci_shade1[4])
    p2 = p2 +
      geom_area(data = subset(df.dens, x >= ci_shaderange4[1] & x <= ci_shaderange4[2]),
                aes(x = x, y = y, fill = as.character(ci_shade1[4])))
  }

  p2 = p2 +
    geom_point(data = data.frame(y = 0,
                                 x = d),
               aes(x = x, y = y),
               size = 3) +
    annotate("segment",
             x = ci_linerange[1],
             xend = ci_linerange[2],
             y = 0, yend = 0,
             size = 1.5,
             colour = "black")


  return(p2)
}


# RMD Check

utils::globalVariables(c("y","subtitle","..density.."))
