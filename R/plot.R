
cut_cdf_qi = function(p, .width = c(.66, .95, 1), labels = NULL) {
  .width = sort(.width)

  if (is.function(labels)) {
    labels = labels(.width)
  } else if (is.null(labels)) {
    labels = .width
  }

  cut(abs(1 - p*2), labels = labels,
      breaks = c(0, .width), include.lowest = TRUE, ordered_result = TRUE)
}

