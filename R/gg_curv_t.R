#' @importFrom tidyr pivot_longer

gg_curv_t <- function(data_list,
                      type = c("c","cd"),
                      levels = c(.5,.90,.95,.999),
                      position = "pyramid",
                      xaxis = expression(theta == ~ "Range of Values"),
                      yaxis1 = expression(paste("two-tailed ", italic(p),
                                                "-value")),
                      yaxis2 = "Confidence Interval (%)",
                      color = "black",
                      fill = "skyblue",
                      alpha_shade = .5
) {

  data = data_list[[1]]
  if (ncol(data) != 7) {
    stop("Error: 'data' or 'list' must be from 'concurve'.")
  }

  if (is.character(position) != TRUE) {
    stop("Error: 'position' must be a string such as 'pyramid' or 'inverted'.")
  }

  if (is.character(fill) != TRUE) {
    stop("Error: 'fill' must be a string for the color.")
  }

  ci_shade1 = sort(levels, decreasing = TRUE)
  interval <- lapply(
    ci_shade1,
    FUN = function(i)
      (c(i,
         data[which(abs(data$intrvl.level -
                          i) == min(abs(data$intrvl.level - i))), ][, 1],
         data[which(abs(data$intrvl.level -
                          i) == min(abs(data$intrvl.level - i))), ][, 2]))
  )

  interval <- data.frame(do.call(rbind, interval))
  interval <- pivot_longer(interval, X2:X3, names_to = "levels", values_to = "limits")
  interval <- interval[, -2]
  colum_names <- c("levels", "limits")
  colnames(interval) <- colum_names

  # Consonance Curve -----------------------------------------------------

  if ("c" %in% type) {
    # Plotting Intervals ------------------------------------------------------
    p_c = ggplot(data = data) +
      geom_line(aes(x = lower.limit, y = pvalue),
                color = color
      ) +
      geom_line(aes(x = upper.limit, y = pvalue),
                color = color
      ) +
      geom_ribbon(aes(x = lower.limit, ymin = min(pvalue), ymax = pvalue),
                  fill = fill, alpha = alpha_shade) +
      geom_ribbon(aes(x = upper.limit, ymin = min(pvalue), ymax = pvalue),
                  fill = fill, alpha = alpha_shade) +
      geom_point(data = interval,
                 mapping = aes(x = limits, y = 1 - levels),
                 size = 1.75, shape = 18) +
      geom_line(data = interval,
                mapping = aes(x = limits, y = 1 - levels, group = levels),
                size = .30) +
      labs(
        x = xaxis,
        y = yaxis1
      ) +
      theme_bw() +
      theme(
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        text = element_text(size = 11)
      ) +
      {
        if (position == "inverted") {
          scale_y_reverse(
            expand = expansion(mult = c(0.01, 0.025)),
            breaks = seq(0, 1, 0.10),
            sec.axis = sec_axis(~ (1 - .) * 100, name = yaxis2, breaks = seq(0, 100, 10))
          )
        }
      } +
      {
        if (position == "pyramid") {
          scale_y_continuous(
            expand = expansion(mult = c(0.01, 0.025)),
            breaks = seq(0, 1, 0.10),
            sec.axis = sec_axis(~ (1 - .) * 100, name = yaxis2, breaks = seq(0, 100, 10))
          )
        }
      }

    # Surprisal Curve ------------------------------------------------------
  }

  if ("cd" %in% type) {
    cdf_dat = data_list[[2]]
    cdf_dat2 = cdf_dat$x

    x.dens  <- density(cdf_dat2)
    df.dens <- data.frame(x=x.dens$x, y=x.dens$y)
    ci_shade1 = sort(levels, decreasing = TRUE)
    interval2 <- lapply(
      ci_shade1,
      FUN = function(i)
        (c(i,
           data[which(abs(data$intrvl.level -
                            i) == min(abs(data$intrvl.level - i))), ][, 1],
           data[which(abs(data$intrvl.level -
                            i) == min(abs(data$intrvl.level - i))), ][, 2]))
    )

    interval2 <- data.frame(do.call(rbind, interval2))
    colnames(interval2) = c("lvl","li","ui")
    #interval <- pivot_longer(interval, X2:X3, names_to = "levels", values_to = "limits")

    p_cd1 = ggplot(data = cdf_dat, mapping = aes(x = x)) +
      geom_density(color = color,
                   alpha = 0.20) +
      geom_area(data = subset(df.dens, x >= interval2$li[1]  & x <= interval2$ui[1]),
                aes(x = x, y = y, fill = as.character(ci_shade1[1]))) +
      scale_fill_brewer(direction = -1,
                        na.translate = FALSE) +
      labs(x = '', y = '',
           fill = "Confidence Interval")

    if(length(ci_shade1 > 1)){

      p_cd2 = p_cd1 +
        geom_area(data = subset(df.dens, x >= interval2$li[2]  & x <= interval2$ui[2]),
                  aes(x = x, y = y, fill = as.character(ci_shade1[2])))
    } else {
      p_cd2 = p_cd1
    }

    if(length(ci_shade1 > 2)){

      p_cd2 = p_cd2 +
        geom_area(data = subset(df.dens, x >= interval2$li[3]  & x <= interval2$ui[3]),
                  aes(x = x, y = y, fill = as.character(ci_shade1[3])))
    }

    if(length(ci_shade1 > 3)){

      p_cd2 = p_cd2 +
        geom_area(data = subset(df.dens, x >= interval2$li[4]  & x <= interval2$ui[4]),
                  aes(x = x, y = y, fill = as.character(ci_shade1[4])))
    }
    p_cd2 = p_cd2 +
      labs(
        x = xaxis,
        y = "Density"
      ) +
      theme_bw() +
      theme(
        legend.position="top",
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        text = element_text(size = 11)
      ) +
      scale_y_continuous(expand = expansion(mult = c(0.01, 0.05)))


    # Relative Likelihood Function -----------------------------------------------------
  }
  if("cd" %in% type && "c" %in% type){
    p1 = plot_grid(p_cd2,
                   p_c,
                   ncol = 1)
  } else
    if ("cd" %in% type){
      p1 = p_cd2
    } else
      if("c" %in% type){
        p1 = p_c
      }

  return(p1)
}

# RMD Check
utils::globalVariables(c("df", "lower.limit", "upper.limit", "intrvl.width", "intrvl.level", "cdf", "pvalue", "svalue"))
utils::globalVariables(c("X2", "X3", "limits", "x"))
