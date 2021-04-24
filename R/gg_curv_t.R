gg_curv_t <- function(data,
                      type = "c",
                      measure = "default",
                      levels = 0.95,
                      nullvalue = NULL,
                      position = "pyramid",
                      xaxis = expression(theta == ~ "Range of Values"),
                      yaxis1 = expression(paste("two-tailed ", italic(p),
                                                "-value")),
                      yaxis2 = "Confidence Interval (%)",
                      color = "black",
                      fill = "skyblue",
                      alpha = .5
) {


  # Consonance Curve -----------------------------------------------------

  if (type == "c") {
    if (ncol(data) != 7) {
      stop("Error: 'data' or 'list' must be from 'concurve'.")
    }
    if (is.character(measure) != TRUE) {
      stop("Error: 'measure' must be a string such as 'default' or 'ratio'.")
    }
    #  if (is.logical(nullvalue) != TRUE | is.numeric(nullvalue) != TRUE) {
    #   stop("Error: 'nullvalue' must be a logical statement such as 'TRUE' or 'FALSE' or a numeric vector.")
    #  }
    if (is.character(position) != TRUE) {
      stop("Error: 'position' must be a string such as 'pyramid' or 'inverted'.")
    }

    if (is.character(fill) != TRUE) {
      stop("Error: 'fill' must be a string for the color.")
    }


    # Plotting Intervals ------------------------------------------------------

    interval <- lapply(levels, FUN = function(i) (c(i,
                                                        subset(data,
                                                               intrvl.level == i)[, 1],
                                                        subset(data,
                                                               intrvl.level == i)[, 2])))
    interval <- data.frame(do.call(rbind, interval))
    interval <- pivot_longer(interval, X2:X3, names_to = "levels", values_to = "limits")
    interval <- interval[, -2]
    colum_names <- c("levels", "limits")
    colnames(interval) <- colum_names


    p1 = ggplot(data = data) +
      geom_line(aes(x = lower.limit, y = pvalue),
                color = color
      ) +
      geom_line(aes(x = upper.limit, y = pvalue),
                color = color
      ) +
      geom_ribbon(aes(x = lower.limit, ymin = min(pvalue), ymax = pvalue),
                  fill = fill, alpha = alpha) +
      geom_ribbon(aes(x = upper.limit, ymin = min(pvalue), ymax = pvalue),
                  fill = fill, alpha = alpha) +
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
      theme_minimal() +
      theme(
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        text = element_text(size = 11)
      ) +
      {
        if (measure == "ratio") scale_x_log10(breaks = scales::pretty_breaks(n = 5))
      } +
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

  if (type == "cd") {
    if (ncol(data) != 1) {
      stop("Error: 'data' must be a data frame from the curve_boot function in 'concurve'.")
    }
    if (is.character(title) != TRUE) {
      stop("Error: 'title' must be a string.")
    }
    if (is.character(subtitle) != TRUE) {
      stop("Error: 'subtitle' must be a string.")
    }
    if (is.character(fill) != TRUE) {
      stop("Error: 'fill' must be a string for the color.")
    }

    ggplot(data = data, mapping = aes(x = x)) +
      geom_density(color = color,
                   alpha = 0.20) +
      labs(
        x = xaxis,
        y = "Density"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14),
        plot.subtitle = element_text(size = 12),
        plot.caption = element_text(size = 8),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        text = element_text(size = 15)
      ) +
      scale_y_continuous(expand = expansion(mult = c(0.01, 0.05)), breaks = scales::pretty_breaks(n = 10))


    # Relative Likelihood Function -----------------------------------------------------
  }
return(p1)
}

# RMD Check
utils::globalVariables(c("df", "lower.limit", "upper.limit", "intrvl.width", "intrvl.level", "cdf", "pvalue", "svalue"))
utils::globalVariables(c("X2", "X3", "limits", "x"))
