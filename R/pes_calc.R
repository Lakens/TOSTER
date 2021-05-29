# pes calculations

pes_ci <- function(Fstat,
                   df1,
                   df2,
                   conf.level = .95){

  pes = Fstat * df1 / (Fstat*df1+df2)

  F_limits <- conf.limits.ncf(F.value = Fstat,
                              df.1 = df1,
                              df.2 = df2,
                              conf.level = conf.level)
  if(all(is.na(F_limits))){
    F_limits$Lower.Limit = 0
    F_limits$Upper.Limit = 0
  }

  LL_lambda <- F_limits$Lower.Limit
  UL_lambda <- F_limits$Upper.Limit

  LL_partial_eta2 <- LL_lambda / (LL_lambda + df1 + df2 + 1)
  UL_partial_eta2 <- UL_lambda / (UL_lambda + df1 + df2 + 1)


  if (is.na(LL_partial_eta2)) {
    LL_partial_eta2 <- 0
  }

  if (is.na(UL_partial_eta2)) {
    UL_partial_eta2 <- 1
  }

  cint = c(LL_partial_eta2 ,UL_partial_eta2)
  return(cint)
}

conf.limits.ncf = function (F.value = NULL,
                            conf.level = 0.95,
                            df.1 = NULL, df.2 = NULL,
                            alpha.lower = NULL,
                            alpha.upper = NULL,
                            tol = 1e-09, Jumping.Prop = 0.1)
{
  if (Jumping.Prop <= 0 | Jumping.Prop >= 1)
    stop("The Jumping Proportion ('Jumping.Prop') must be between zero and one.")
  if (is.null(F.value))
    stop("Your 'F.value' is not correctly specified.")
  if (F.value < 0)
    stop("Your 'F.value' is not correctly specified.")
  if (is.null(df.1) | is.null(df.2))
    stop("You must specify the degrees of freedom ('df.1' and 'df.2').")
  if (is.null(alpha.lower) & is.null(alpha.upper) & is.null(conf.level))
    stop("You need to specify the confidence interval parameters.")
  if ((!is.null(alpha.lower) | !is.null(alpha.upper)) & !is.null(conf.level))
    stop("You must specify only one method of defining the confidence limits.")
  if (!is.null(conf.level)) {
    if (conf.level >= 1 | conf.level <= 0)
      stop("Your confidence level ('conf.level') must be between 0 and 1.")
    alpha.lower <- alpha.upper <- (1 - conf.level)/2
  }
  if (alpha.lower == 0)
    alpha.lower <- NULL
  if (alpha.upper == 0)
    alpha.upper <- NULL
  FAILED <- NULL
  if (!is.null(alpha.lower)) {
    LL.0 <- qf(p = alpha.lower * 5e-04, df1 = df.1, df2 = df.2)
    Diff <- pf(q = F.value, df1 = df.1, df2 = df.2, ncp = LL.0) -
      (1 - alpha.lower)
    if (pf(q = F.value, df1 = df.1, df2 = df.2, ncp = LL.0) <
        (1 - alpha.lower)) {
      FAILED <- if (pf(q = F.value, df1 = df.1, df2 = df.2,
                       ncp = 0) < 1 - alpha.lower)
        LL.0 <- 1e-08
      if (pf(q = F.value, df1 = df.1, df2 = df.2, ncp = LL.0) <
          1 - alpha.lower)
        FAILED <- TRUE
    }
    if (is.null(FAILED)) {
      LL.1 <- LL.2 <- LL.0
      while (Diff > tol) {
        LL.2 <- LL.1 * (1 + Jumping.Prop)
        Diff <- pf(q = F.value, df1 = df.1, df2 = df.2,
                   ncp = LL.2) - (1 - alpha.lower)
        LL.1 <- LL.2
      }
      LL.1 <- LL.2/(1 + Jumping.Prop)
      LL.Bounds <- c(LL.1, (LL.1 + LL.2)/2, LL.2)
      Diff <- pf(q = F.value, df1 = df.1, df2 = df.2, ncp = LL.Bounds[2]) -
        (1 - alpha.lower)
      while (abs(Diff) > tol) {
        Diff.1 <- pf(q = F.value, df1 = df.1, df2 = df.2,
                     ncp = LL.Bounds[1]) - (1 - alpha.lower) > tol
        Diff.2 <- pf(q = F.value, df1 = df.1, df2 = df.2,
                     ncp = LL.Bounds[2]) - (1 - alpha.lower) > tol
        Diff.3 <- pf(q = F.value, df1 = df.1, df2 = df.2,
                     ncp = LL.Bounds[3]) - (1 - alpha.lower) > tol
        if (Diff.1 == TRUE & Diff.2 == TRUE & Diff.3 ==
            FALSE) {
          LL.Bounds <- c(LL.Bounds[2], (LL.Bounds[2] +
                                          LL.Bounds[3])/2, LL.Bounds[3])
        }
        if (Diff.1 == TRUE & Diff.2 == FALSE & Diff.3 ==
            FALSE) {
          LL.Bounds <- c(LL.Bounds[1], (LL.Bounds[1] +
                                          LL.Bounds[2])/2, LL.Bounds[2])
        }
        Diff <- pf(q = F.value, df1 = df.1, df2 = df.2,
                   ncp = LL.Bounds[2]) - (1 - alpha.lower)
      }
      LL <- LL.Bounds[2]
    }
  }
  if (!is.null(FAILED))
    LL <- NA
  if (!is.null(alpha.upper)) {
    FAILED.Up <- NULL
    UL.0 <- qf(p = 1 - alpha.upper * 5e-04, df1 = df.1, df2 = df.2)
    Diff <- pf(q = F.value, df1 = df.1, df2 = df.2, ncp = UL.0) -
      alpha.upper
    if (Diff < 0)
      UL.0 <- 1e-08
    Diff <- pf(q = F.value, df1 = df.1, df2 = df.2, ncp = UL.0) -
      alpha.upper
    if (Diff < 0) {
      FAILED.Up <- TRUE
    }
    if (is.null(FAILED.Up)) {
      UL.1 <- UL.2 <- UL.0
      while (Diff > tol) {
        UL.2 <- UL.1 * (1 + Jumping.Prop)
        Diff <- pf(q = F.value, df1 = df.1, df2 = df.2,
                   ncp = UL.2) - alpha.upper
        UL.1 <- UL.2
      }
      UL.1 <- UL.2/(1 + Jumping.Prop)
      UL.Bounds <- c(UL.1, (UL.1 + UL.2)/2, UL.2)
      Diff <- pf(q = F.value, df1 = df.1, df2 = df.2, ncp = UL.Bounds[2]) -
        alpha.upper
      while (abs(Diff) > tol) {
        Diff.1 <- pf(q = F.value, df1 = df.1, df2 = df.2,
                     ncp = UL.Bounds[1]) - alpha.upper > tol
        Diff.2 <- pf(q = F.value, df1 = df.1, df2 = df.2,
                     ncp = UL.Bounds[2]) - alpha.upper > tol
        Diff.3 <- pf(q = F.value, df1 = df.1, df2 = df.2,
                     ncp = UL.Bounds[3]) - alpha.upper > tol
        if (Diff.1 == TRUE & Diff.2 == TRUE & Diff.3 ==
            FALSE) {
          UL.Bounds <- c(UL.Bounds[2], (UL.Bounds[2] +
                                          UL.Bounds[3])/2, UL.Bounds[3])
        }
        if (Diff.1 == TRUE & Diff.2 == FALSE & Diff.3 ==
            FALSE) {
          UL.Bounds <- c(UL.Bounds[1], (UL.Bounds[1] +
                                          UL.Bounds[2])/2, UL.Bounds[2])
        }
        Diff <- pf(q = F.value, df1 = df.1, df2 = df.2,
                   ncp = UL.Bounds[2]) - alpha.upper
      }
      UL <- UL.Bounds[2]
    }
    if (!is.null(FAILED.Up))
      UL <- NA
  }
  if (!is.null(alpha.lower) & !is.null(alpha.upper))
    return(list(Lower.Limit = LL, Prob.Less.Lower = 1 - pf(q = F.value,
                                                           df1 = df.1, df2 = df.2, ncp = LL), Upper.Limit = UL,
                Prob.Greater.Upper = pf(q = F.value, df1 = df.1,
                                        df2 = df.2, ncp = UL)))
  if (is.null(alpha.lower) & !is.null(alpha.upper))
    return(list(Upper.Limit = UL, Prob.Greater.Upper = pf(q = F.value,
                                                          df1 = df.1, df2 = df.2, ncp = UL)))
  if (!is.null(alpha.lower) & is.null(alpha.upper))
    return(list(Lower.Limit = LL, Prob.Less.Lower = 1 - pf(q = F.value,
                                                           df1 = df.1, df2 = df.2, ncp = LL)))
}

pes_curv = function (Fstat,
                     df1,
                     df2,
                     steps = 5000) {
  intrvls <- (0:steps)/steps
  intrvls = subset(intrvls,intrvls>0 & intrvls<1)

  # Confidence interval of the SMD from Goulet-Pelletier & Cousineau
  results <-
    suppressWarnings({
      lapply(
        intrvls,
        FUN = function(i)
          pes_ci(
            Fstat = Fstat,
            df1 =
              df1,
            df2 =
              df2,
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
