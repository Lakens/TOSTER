CIr <- function (r, n, conf.level = 0.95){
  z <- r2z(r)
  uciz <- CIz(z, n, conf.level)[2]
  lciz <- CIz(z, n, conf.level)[1]
  ur <- z2r(uciz)
  lr <- z2r(lciz)
  mat <- list(lr, ur)
  return(as.numeric(mat))
}

CIz = function (z, n, conf.level = 0.95){
  noma <- 1 - conf.level
  sez <- SEz(n)
  zs <- -qnorm(noma/2)
  mez <- zs * sez
  lcl <- z - mez
  ucl <- z + mez
  mat <- list(lcl, ucl)
  return(as.numeric(mat))
}

z2r = function (x) {
  (exp(2 * x) - 1)/(exp(2 * x) + 1)
}

corr_curv = function (r,
                     n,
                     steps = 5000) {
  intrvls <- (0:steps)/steps
  intrvls = subset(intrvls,intrvls>0 & intrvls<1)

  # Confidence interval of the SMD from Goulet-Pelletier & Cousineau
  results <-
    suppressWarnings({
      lapply(
        intrvls,
        FUN = function(i)
          CIr(
            r = r,
            n = n,
            conf.level = i
          )
      )
    })

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
