#' Methods for TOSTt objects
#'
#' Methods defined for objects returned from the agree functions.
#'
#' @param x object of class \code{TOSTt} as returned from the reli_stats function
#' @param digits Number of digits to print for p-values
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
#' @export

plot.TOSTt <- function(x,  ...){

  points = data.frame(
    type = c("Mean Difference", x$smd$smd_label),
    mu = c(x$effsize$estimate[1], 0),
    param = c(round(unname(x$TOST$df[1]), 0), round(unname(x$smd$d_df), 0)),
    sigma = c(x$TOST$SE[1], x$smd$d_sigma),
    lambda = c(0, x$smd$d_lambda),
    low = c(x$eqb$low_eq[1], x$eqb$low_eq[2]),
    high = c(x$eqb$high_eq[1], x$eqb$high_eq[2]),
    alpha = c(x$alpha, x$alpha),
    stringsAsFactors = FALSE
  )

  c1 = 1-x$alpha
  c2 = 1-x$alpha*2
  if(c1 < .999 && c2 > .5){
    sets = c(.5,c2,c1,.999)
  } else if(c2 <=.5 && c1 < .999) {
    sets = c(c2,c1,.999)
  } else {
    sets = c(.5,c2,c1)
  }
  #if
  #sets = ci.cuts
  p1 = ggplot(data = points,
                aes_string(y = 0)) +
    stat_dist_halfeye(aes(
      dist = dist_student_t(
        mu = mu,
        df = param,
        sigma = sigma,
        ncp = lambda
      ),
      fill = stat(cut_cdf_qi(p=cdf,
                             .width = sets))
    ),
    .width = c(c2, c1)) +
    scale_fill_brewer(direction = -1,
                      na.translate = FALSE) +
    labs(x = '', y = '',
         fill = "Confidence Interval") +
    geom_vline(aes(xintercept = low),
               linetype="dashed") +
    geom_vline(aes(xintercept = high),
               linetype="dashed") +
    geom_text(aes(y=1.5, x=low,
                  vjust=-.9, hjust=1),
              angle = 90,
              label='Lower Bound') +
    geom_text(aes(y=1.5, x=high, vjust=1.5, hjust=1),
              angle = 90,
              label='Upper Bound') +
    theme_tidybayes() +
    theme(legend.position="top",
          strip.text = element_text(face="bold", size=10),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    facet_wrap(~type,
               ncol = 1,
               scales = "free")

  return(p1)

}
