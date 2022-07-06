

model.frame(extra ~ 1, data = data1)
compare_t = function(x1,
                     y1 = NULL,
                     x2,
                     y2 = NULL,
                     paired = FALSE,
                     alternative = "two.sided",
                     R = 1999,
                     alpha = 0.05){

  if(!missing(null) && (length(null) != 1 || is.na(null))) {
    stop("'null' must be a single number")
  }


  if(!missing(alpha) && (length(alpha) != 1 || !is.finite(alpha) ||
                         alpha < 0 || alpha > 1)) {
    stop("'alpha' must be a single number between 0 and 1")
  }

  if(TOST){
    conf_level = 1-alpha*2
  } else {
    conf_level = 1-alpha
  }

  if(is.null(y1) && is.null(y2)){
    paired = TRUE
  }

  if(paired){
    if(is.null(y1) && is.null(y2)){
      df1 = na.omit(data.frame(z = x1))
      df2 = na.omit(data.frame(z = x2))
    } else {
      df1 = na.omit(data.frame(z = x1 - y1))
      df2 = na.omit(data.frame(z = x2 - y2))
    }
  } else {
    x1 = na.omit(x1)
    y1 = na.omit(y1)
    x2 = na.omit(x2)
    y2 = na.omit(y2)
    df1 = data.frame(y = c(x1,y1),
                     group = c(rep("x",length(x1)),
                               rep("y",length(y1))))
    df2 = data.frame(y = c(x2,y2),
                     group = c(rep("x",length(x2)),
                               rep("y",length(y2))))
  }

  d_vec <- rep(NA, times=length(R)) # smd vector
  m_vec <- rep(NA, times=length(R)) # mean difference vector
}



