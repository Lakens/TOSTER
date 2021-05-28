# Plot scripts

t_curv = function (TOST_res,
                   steps = 5000) {
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
                   steps = 5000) {
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

d_curv_raw = function (d,
                       df,
                       lambda,
                       steps = 5000) {
  intrvls <- (0:steps)/steps

  # Confidence interval of the SMD from Goulet-Pelletier & Cousineau
  results <- suppressWarnings({lapply(intrvls, FUN = function(i) d_CI(d = d,
                                                                      df = df,
                                                                      lambda = lambda,
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

