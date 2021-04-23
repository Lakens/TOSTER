
harm_mean = function(n1,n2){
ntilde = (2*n1*n2)/(n1+n2)
return(ntilde)
}

d_est_pair <- function(n,
                       m1,
                       m2,
                       sd1,
                       sd2,
                       r12,
                       type = "g",
                       denom = "z",
                       alpha = .05){

  sdif <- sqrt(sd1 ^ 2 + sd2 ^ 2 - 2 * r12 * sd1 * sd2)
  if(denom == "z"){
    d_denom = sdif
  } else if(denom == "rm"){
    d_denom = sdif * sqrt(2*(1-r12))
  }

  df <- n-1

  cohend = abs(m2-m1) / d_denom

  d_df = 2*(n)-2
  #J <- gamma(df / 2) / (sqrt(df / 2) * gamma((df - 1) / 2))
  J = gamma(d_df / 2) / (sqrt(d_df / 2) * gamma((d_df - 1) / 2))

  if(denom == "z"){
    if (type == 'g') {
      cohend <-  cohend * J
      smd_label = "Hedges' g(z)"
    } else {
      smd_label = "Cohen's d(z)"
    }
  } else if(denom == "rm"){
    if (type == 'g') {
      cohend <-  cohend * J
      smd_label = "Hedges' g(rm)"
    } else {
      smd_label = "Cohen's d(rm)"
    }
  }

  d_lambda <- cohend * sqrt(n / (2*(1 - r12)))

  # Equation 4b Goulet-Pelletier and Cousineau, 2018
  # ((2*(1-r12))/n)
  d_sigma = sqrt((d_df/(d_df-2)) * ((2*(1-r12))/n)*(1+cohend^2*(n/(2*(1-r12)))) - cohend^2/J^2)
  if(denom == "z"){
    d_sigma = d_sigma * sqrt(2*(1-r12))
  }
  # Get cohend confidence interval
  tlow <- suppressWarnings(qt(1 / 2 - (1-alpha*2) / 2,
                              df = d_df,
                              ncp = d_lambda))
  thigh <- suppressWarnings(qt(1 / 2 + (1-alpha*2) / 2,
                               df = d_df,
                               ncp = d_lambda))
  if (m1 == m2) {
    dlow <- (tlow * (2 * n)) / (sqrt(d_df) * sqrt(2 * n))
    dhigh <- (thigh * (2 * n)) / (sqrt(d_df) * sqrt(2 * n))
  } else {
    dlow <- tlow / d_lambda * cohend
    dhigh <- thigh / d_lambda * cohend
  }

  if (m2-m1 < 0) {
    cohend <- cohend * -1
    tdlow <- dlow
    dlow <- dhigh * -1
    dhigh <- tdlow * -1
    d_lambda <- cohend * sqrt(n / (2*(1 - r12)))
  }

  return(list(
    d = cohend,
    d_df = d_df,
    dlow = dlow,
    dhigh = dhigh,
    d_sigma = d_sigma,
    d_lambda = d_lambda,
    smd_label = smd_label,
    J = J,
    d_denom = d_denom,
    ntilde = n,
    r12 = r12
  ))
}


d_est_ind <- function(n1,
                      n2,
                      m1,
                      m2,
                      sd1,
                      sd2,
                      type = "g",
                      var.equal = TRUE,
                      alpha = .05){


  if (var.equal) {
    denomSD <- sqrt((((n1 - 1)*(sd1^2)) + (n2 - 1)*(sd2^2))/((n1+n2)-2)) #calculate sd pooled
    d_df = n1 + n2 - 2
  } else {
    denomSD <- sqrt((sd1^2 + sd2^2)/2) #calculate sd root mean squared for Welch's t-test
    d_df1 = (n1 - 1)*(n2 - 1)*(sd1^2+sd2^2)^2
    d_df2 = (n2-1)*sd1^4+(n1-1)*sd2^4
    d_df = d_df1/d_df2
  }

  denomSD[is.na(denomSD)] <- NaN

  #denomSD <- jmvcore::tryNaN(sqrt(((n1-1)*v[1]+(n2-1)*v[2])/(n1+n2-2)))
  d <- abs(m1-m2)/denomSD # Cohen's d

  d[is.na(d)] <- NaN

  cohend = d
  ntilde <- harm_mean(n1,n2)

  # Compute unbiased Hedges' g
  # Use the lgamma function, and update to what Goulet-Pelletier & Cousineau used; works with larger inputs
  J <- gamma(d_df/2)/(sqrt(d_df/2)*gamma((d_df-1)/2))

  if(var.equal == TRUE){
    if(type == 'g') {
      cohend <-  cohend * J
      smd_label = "Hedges' g"
    } else {
      smd_label = "Cohen's d"
    }
  } else {
    if(type == 'g') {
      cohend <-  cohend * J
      smd_label = "Hedges' g(av)"
    } else {
      smd_label = "Cohen's d(av)"
    }
  }


  # add options for cohend here
  d_lambda <- cohend * sqrt(ntilde/2)
  #d_sigma = sqrt((n1+n2)/(n1*n2)+(cohend^2/(2*(n1+n2))))
  d_sigma = sqrt((d_df/(d_df-2)) * (2/ntilde) *(1+cohend^2*(ntilde/2)) - cohend^2/J^2)

  # Confidence interval of the SMD from Goulet-Pelletier & Cousineau
  tlow <- qt(1 / 2 - (1-alpha*2) / 2, df = d_df, ncp = d_lambda)
  thigh <- qt(1 / 2 + (1-alpha*2) / 2, df = d_df, ncp = d_lambda)
  if(m1 == m2) {
    dlow <- (tlow*(n1+n2)) / (sqrt(d_df) * sqrt(n1*n2))
    dhigh <- (thigh*(n1+n2)) / (sqrt(d_df) * sqrt(n1*n2))
  } else {
    dlow <- tlow / d_lambda * cohend
    dhigh <- thigh / d_lambda * cohend
  }

  # The function provided by Goulet-Pelletier & Cousineau works with +mdiff.
  # Now we fix the signs and directions for negative differences
  if ((m1-m2) < 0) {
    cohend <- cohend * -1
    tdlow <- dlow
    dlow <- dhigh * -1
    dhigh <- tdlow * -1
    d_lambda <- cohend * sqrt(ntilde/2)
  }


  return(list(
    d = cohend,
    d_df = d_df,
    dlow = dlow,
    dhigh = dhigh,
    d_sigma = d_sigma,
    d_lambda = d_lambda,
    smd_label = smd_label,
    J = J,
    d_denom = denomSD,
    ntilde = ntilde
  ))
}

d_CI = function(d,
                df,
                lambda,
                alpha){
  d_lambda = lambda
  d_df = df
  cohend = d
  tlow <- suppressWarnings({
    qt(1 / 2 - (1 - alpha) / 2,
       df = d_df,
       ncp = d_lambda)
  })
  thigh <- suppressWarnings({
    qt(1 / 2 + (1 - alpha) / 2,
       df = d_df,
       ncp = d_lambda)
  })

  dlow <- tlow / d_lambda * cohend
  dhigh <- thigh / d_lambda * cohend

  return(c(dlow,dhigh))
}

d_est_one <- function(n,
                      mu,
                      sd,
                      testValue,
                      type = "g",
                      alpha = .05){

  cohend <- abs(mu-testValue)/sd # Cohen's d
  df <- n-1
  d_df = df
  # Compute unbiased Hedges' g
  # Use the lgamma function, and update to what Goulet-Pelletier & cousineau used; works with larger inputs
  J <- gamma(df/2)/(sqrt(df/2)*gamma((df-1)/2))

  if(type == 'g') {
    cohend <-  cohend * J
    smd_label = "Hedges' g"
  } else {
    smd_label = "Cohen's d"
  }

  d_lambda <- cohend * sqrt(n)
  #d_sigma = sqrt((df + 1)/(df - 1)*(2/n)*(1 + cohend^2/8))
  d_sigma = sqrt((df/(df-2)) * (1/n) *(1+cohend^2*(n/1)) - cohend^2/J^2)

  # Confidence interval of the SMD from Goulet-Pelletier & Cousineau
  tlow <- suppressWarnings(qt(1 / 2 - (1-alpha*2) / 2,
                              df =   d_df,
                              ncp = d_lambda))
  thigh <- suppressWarnings(qt(1 / 2 + (1-alpha*2) / 2,
                               df =   d_df,
                               ncp = d_lambda))
  if((mu-testValue) == 0) {
    dlow <- (tlow*(2*n)) / (sqrt(d_df) * sqrt(2*n))
    dhigh <- (thigh*(2*n)) / (sqrt(d_df) * sqrt(2*n))
  } else {
    dlow <- tlow / d_lambda * cohend
    dhigh <- thigh / d_lambda * cohend
  }

  if ((mu-testValue) < 0) {
    cohend <- cohend * -1
    tdlow <- dlow
    dlow <- dhigh * -1
    dhigh <- tdlow * -1
  }



  return(list(
    d = cohend,
    d_df = d_df,
    dlow = dlow,
    dhigh = dhigh,
    d_sigma = d_sigma,
    d_lambda = d_lambda,
    smd_label = smd_label,
    J = J,
    d_denom = sd,
    ntilde = n
  ))
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

t_CI = function(TOST_res,
                alpha){
  est = TOST_res$effsize$estimate[1]
  SE = TOST_res$effsize$SE[1]
  df = TOST_res$TOST$df[1]
  tlow = est - stats::qt(1-alpha/2,df)*SE
  thigh = est + stats::qt(1-alpha/2,df)*SE
  return(c(tlow,thigh))
}


