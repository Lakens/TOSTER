
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
                       alpha = .05,
                       smd_ci = "goulet"){

  sdif <- sqrt(sd1 ^ 2 + sd2 ^ 2 - 2 * r12 * sd1 * sd2)
  if(denom == "z"){
    d_denom = sdif
  } else if(denom == "rm"){
    d_denom =  (sqrt(2*(1-r12)) / sdif)^-1
  } else if (denom == "glass1"){
    d_denom = sd1
  } else if (denom == "glass2"){
    d_denom = sd2
  }

  df <- n-1
  hn <- 1 / n
  cohend = abs(m1-m2) / d_denom
  if(smd_ci == "goulet"){
    d_df = 2*(n)-2
  } else {
    d_df = n-1
  }

  #J <- gamma(df / 2) / (sqrt(df / 2) * gamma((df - 1) / 2))
  if(type == "g"){
    J = hedge_J(d_df)
  } else {
    J = 1
  }
  cohend <-  cohend * J

  if(denom == "z"){
    if (type == 'g') {
      #cohend <-  cohend * J
      smd_label = "Hedges's g(z)"
    } else {
      smd_label = "Cohen's d(z)"
    }
  } else if(denom == "rm"){
    if (type == 'g') {
      #cohend <-  cohend * J
      smd_label = "Hedges's g(rm)"
    } else {
      smd_label = "Cohen's d(rm)"
    }
  } else if(denom %in% c("glass1","glass2")){
    if (type == 'g') {
      #cohend <-  cohend * J
      smd_label = "Glass's delta(g)"
    } else {
      smd_label = "Glass's delta(d)"
    }
  }

  d_lambda <- cohend * sqrt(n / (2*(1 - r12)))

  # Equation 4b Goulet-Pelletier and Cousineau, 2018
  # ((2*(1-r12))/n)
  d_sigma = sqrt((d_df/(d_df-2)) * ((2*(1-r12))/n)*(1+cohend^2*(n/(2*(1-r12)))) - cohend^2/J^2)
  if(denom == "z"){
    d_sigma = d_sigma * sqrt(2*(1-r12))
  }
  if(smd_ci == "nct" && denom != "rm"){
    #d_sigma = d_denom / sqrt(n)
    d_sigma = sqrt(1/n + (cohend^2/(2*n)))
  }
  if(denom %in% c("glass1","glass2")){
    #sep1 = (n-1)/(n*(n-3))
    #sep2 = (2*(1-r12)+cohend^2*n)
    #sep3 = cohend^2/(J)^2
    # Borenstein 2009 --- adopted from metafor
    d_s1 = J^2*(2*(((1-r12)/n)+((cohend^2*J^(-1))/(2*n))))
      #sqrt(sep1*sep2-sep3)
    d_sigma = sqrt(d_s1)
  }

  if(smd_ci == "goulet"){

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
  }

  if(smd_ci == "nct"){
    if(denom == "rm"){
      t_stat <- (abs(m1 - m2) / (sqrt((sd1^2 + sd2^2)-(2*r12*sd1*sd2))/sqrt(n))) * sqrt(2*(1-r12))
    }else{
      SE1 <- d_denom / sqrt(n)
      t_stat = abs(m1 - m2) / SE1
    }

    ts <- get_ncp_t2(t_stat, d_df, conf.level = 1-alpha*2)
    dlow <- ts[1] * sqrt(hn) * J
    dhigh <- ts[2] * sqrt(hn) * J
  } else{
    t_stat = NULL
    hn = NULL
  }

  if(smd_ci == "t"){
    dlow <- cohend - qt(1-alpha,d_df)*d_sigma
    dhigh <- cohend + qt(1-alpha,d_df)*d_sigma
  }

  if(smd_ci == "z"){
    dlow <- cohend - qnorm(1-alpha)*d_sigma
    dhigh <- cohend + qnorm(1-alpha)*d_sigma
  }

  if (m1-m2 < 0) {
    cohend <- cohend * -1
    tdlow <- dlow
    dlow <- dhigh * -1
    dhigh <- tdlow * -1
    d_lambda <- cohend * sqrt(n / (2*(1 - r12)))
    t_stat = -1*t_stat
  }
  if(smd_ci != "goulet"){
    d_lambda = hn
  }
  return(list(
    d = cohend,
    d_df = d_df,
    dlow = dlow,
    dhigh = dhigh,
    d_sigma = d_sigma,
    d_lambda = d_lambda,
    #hn = hn,
    smd_label = smd_label,
    J = J,
    d_denom = d_denom,
    ntilde = n,
    r12 = r12,
    t_stat = t_stat,
    smd_ci = smd_ci
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
                      alpha = .05,
                      denom = "d",
                      smd_ci = "goulet"){


  if (var.equal) {
    denomSD <- sqrt((((n1 - 1)*(sd1^2)) + (n2 - 1)*(sd2^2))/((n1+n2)-2)) #calculate sd pooled
    d_df = n1 + n2 - 2
    hn <- (1 / n1 + 1 / n2)
  } else {
    denomSD <- sqrt((sd1^2 + sd2^2)/2) #calculate sd root mean squared for Welch's t-test
    d_df1 = (n1 - 1)*(n2 - 1)*(sd1^2+sd2^2)^2
    d_df2 = (n2-1)*sd1^4+(n1-1)*sd2^4
    d_df = d_df1/d_df2
    hn <- (2 * (n2 * sd1^2 + n1 * sd2^2)) / (n1 * n2 * (sd1^2 + sd2^2))
  }

  if (denom == "glass1"){
    denomSD <- sd1
    d_df = n1 -1
    hn <- 1 / n2 + denomSD^2 / (n1 * denomSD^2)
    n_glass = n1
    nn_glass = n2
    sdn_glass = sd2
  } else if (denom == "glass2"){
    denomSD <- sd2
    hn <- 1 / n2 + denomSD^2 / (n1 * denomSD^2)
    d_df = n2 - 1
    n_glass = n2
    nn_glass = n1
    sdn_glass = sd1
  }


  denomSD[is.na(denomSD)] <- NaN

  #denomSD <- jmvcore::tryNaN(sqrt(((n1-1)*v[1]+(n2-1)*v[2])/(n1+n2-2)))
  d <- abs(m1-m2)/denomSD # Cohen's d

  #d[is.na(d)] <- NaN

  cohend = d
  ntilde <- harm_mean(n1,n2)

  # Compute unbiased Hedges's g
  # Use the lgamma function, and update to what Goulet-Pelletier & Cousineau used; works with larger inputs
  if(type == "g"){
    J = hedge_J(d_df)
  } else {
    J = 1
  }


  if(denom %in% c("glass1","glass2")){
    if(type == 'g') {
      cohend <-  cohend * J
      smd_label = "Glass's delta(g)"
    } else {
      smd_label = "Glass's delta(d)"
    }
  } else if(var.equal == TRUE){
    if(type == 'g') {
      cohend <-  cohend * J
      smd_label = "Hedges's g"
    } else {
      smd_label = "Cohen's d"
    }
  } else {
    if(type == 'g') {
      cohend <-  cohend * J
      smd_label = "Hedges's g(av)"
    } else {
      smd_label = "Cohen's d(av)"
    }
  }



  if(var.equal == TRUE && !(denom %in% c("glass1","glass2"))){
    mult_lamb = sqrt((n1*n2*(sd1^2 + sd2^2))/(2*(n2*sd1^2 + n1*sd2^2)))
    d_lambda = cohend * mult_lamb
  } else if(denom %in% c("glass1","glass2")){
    d_lambda <- cohend * sqrt(ntilde/2)
  } else {
    d_lambda <- cohend * sqrt(ntilde/2)
  }

  # add options for cohend here
  if(smd_ci == "goulet"){
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
  } else{
    if (denom %in% c("glass1", "glass2")) {
      N = n1 + n2
      # morris and deshon 2002
      # d_sigma = sqrt((1 / ntilde) * ((N - 2) / (N - 4)) * (1 + ntilde *
      #                                                        cohend ^ 2) - cohend ^ 2 / J ^ 2)
      # Algina, Keselman, and Penfield (2006) from Delacre et al 2021
      d_sigma2 = d_df / (d_df -2) * (1/n_glass + sdn_glass^2/(nn_glass*denomSD^2))+ cohend^2 * (d_df/(d_df-2)-J^2)
      d_sigma = sqrt(d_sigma2)
    } else {
      if (var.equal) {
        d_sigma = sqrt(((n1 + n2) / (n1 * n2) + d ^ 2 / (2 * (n1 + n2))) * J ^ 2)
      } else{
        par1 = 2*(sd1^2/n1+sd2^2/n2)/(sd1^2+sd2^2)
        par2 = d_df/(d_df-2)-J^2
        d_sigma = sqrt(d_df/(d_df-2)*par1+cohend^2*par2)
      }
    }

  }

  if(smd_ci == "nct"){
    if( !(denom %in% c("glass1","glass2"))){
      #d_sigma = denomSD * sqrt(1 / n1 + 1 / n2)
      if(var.equal){
        SE1 = denomSD * sqrt(1 / n1 + 1 / n2)
      } else {
        se1 <- sqrt(sd1^2 / n1)
        se2 <- sqrt(sd2^2 / n2)
        SE1 <- sqrt(se1^2 + se2^2)
      }

      t_stat = abs(m1-m2)/SE1
      ts <- get_ncp_t2(t_stat, d_df, conf.level = 1-alpha*2)
    } else {
      SE1 = (denomSD * sqrt(1 / n_glass + sdn_glass^2 / (nn_glass * denomSD^2)))
      d_df <- n1+n2 - 2
      t_stat = abs(m1-m2)/SE1
      ts <- get_ncp_t2(t_stat, d_df, conf.level = 1-alpha*2)
    }
    dlow <- ts[1] * sqrt(hn) * J
    dhigh <- ts[2] * sqrt(hn) * J
  } else {
    t_stat = NULL
    hn = NULL
  }

  if(smd_ci == "t"){
    dlow <- cohend - qt(1-alpha,d_df)*d_sigma
    dhigh <- cohend + qt(1-alpha,d_df)*d_sigma
  }

  if(smd_ci == "z"){
    dlow <- cohend - qnorm(1-alpha)*d_sigma
    dhigh <- cohend + qnorm(1-alpha)*d_sigma
  }

  if ((m1-m2) < 0) {
    cohend <- cohend * -1
    tdlow <- dlow
    dlow <- dhigh * -1
    dhigh <- tdlow * -1
    t_stat = -1*t_stat
    if(var.equal == TRUE && !(denom %in% c("glass1","glass2"))){
      mult_lamb = sqrt((n1*n2*(sd1^2 + sd2^2))/(2*(n2*sd1^2 + n1*sd2^2)))
      d_lambda = cohend * mult_lamb
    } else if(denom %in% c("glass1","glass2")){
      d_lambda <- cohend * sqrt(ntilde/2)
    } else {
      d_lambda <- cohend * sqrt(ntilde/2)
    }
  }

  if(smd_ci != "goulet"){
    d_lambda = hn
  }
  return(list(
    d = cohend,
    d_df = d_df,
    dlow = dlow,
    dhigh = dhigh,
    d_sigma = d_sigma,
    d_lambda = d_lambda,
    #hn = hn,
    smd_label = smd_label,
    J = J,
    d_denom = denomSD,
    ntilde = ntilde,
    t_stat = t_stat,
    smd_ci = smd_ci
  ))
}


d_CI = function(d,
                df,
                lambda,
                sigma,
                t_stat,
                #hn,
                alpha,
                smd_ci = "goulet"){
  d_lambda = lambda
  d_df = df
  cohend = d
  d_sigma = sigma
  t_stat = t_stat
  #hn = hn

  if(smd_ci == "goulet"){
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
  } else{
    hn = d_lambda
  }

  if(smd_ci == "nct"){
    ts <- get_ncp_t2(t_stat, d_df, conf.level = 1 / 2 + (1-alpha*2) / 2)
    #MBESS::conf.limits.nct(t_stat, d_df, conf.level = (1 - alpha),
    #                       sup.int.warns = TRUE)
    dlow <- ts[1] * sqrt(hn)
    dhigh <- ts[2] * sqrt(hn)
  }

  if(smd_ci == "t"){
    dlow <- cohend - qt(1-alpha/2,d_df)*d_sigma
    dhigh <- cohend + qt(1-alpha/2,d_df)*d_sigma
  }

  if(smd_ci == "z"){
    dlow <- cohend - qnorm(1-alpha/2)*d_sigma
    dhigh <- cohend + qnorm(1-alpha/2)*d_sigma
  }

  return(c(dlow,dhigh))
}

d_est_one <- function(n,
                      mu,
                      sd,
                      testValue,
                      type = "g",
                      alpha = .05,
                      smd_ci = "goulet"){

  cohend <- abs(mu-testValue)/sd # Cohen's d
  df <- n-1
  d_df = df
  hn <- 1 / n
  # Compute unbiased Hedges' g
  # Use the lgamma function, and update to what Goulet-Pelletier & cousineau used; works with larger inputs
  if(type == "g"){
    J = hedge_J(d_df)
  } else {
    J = 1
  }

  if(type == 'g') {
    cohend <-  cohend * J
    smd_label = "Hedges's g"
  } else {
    smd_label = "Cohen's d"
  }

  d_lambda <- cohend * sqrt(n)
  if(smd_ci == "goulet"){
    d_sigma = sqrt((df/(df-2)) * (1/n) *(1+cohend^2*(n/1)) - cohend^2/J^2)
  } else {
    d_sigma = sqrt(1/n + (cohend^2/(2*n)))
  }
  if(smd_ci == "goulet"){
  #d_sigma = sqrt((df + 1)/(df - 1)*(2/n)*(1 + cohend^2/8))

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
  }

  if(smd_ci == "nct"){
    SE1 <- sd / sqrt(n)
    t_stat = abs(mu-testValue) / SE1
    ts <- get_ncp_t2(t_stat, d_df, conf.level = 1-alpha*2)
    dlow <- ts[1] * sqrt(hn)*J
    dhigh <- ts[2] * sqrt(hn)*J

  } else {
    t_stat = NULL
    hn = NULL
  }

  if(smd_ci == "t"){
    dlow <- cohend - qt(1-alpha/2,d_df)*d_sigma
    dhigh <- cohend + qt(1-alpha/2,d_df)*d_sigma
  }

  if(smd_ci == "z"){
    dlow <- cohend - qnorm(1-alpha/2)*d_sigma
    dhigh <- cohend + qnorm(1-alpha/2)*d_sigma
  }

  if ((mu-testValue) < 0) {
    cohend <- cohend * -1
    tdlow <- dlow
    dlow <- dhigh * -1
    dhigh <- tdlow * -1
    t_stat = -1*t_stat
  }

  if(smd_ci != "goulet"){
    d_lambda = hn
  }

  return(list(
    d = cohend,
    d_df = d_df,
    dlow = dlow,
    dhigh = dhigh,
    d_sigma = d_sigma,
    d_lambda = d_lambda,
    #hn = hn,
    smd_label = smd_label,
    J = J,
    d_denom = sd,
    ntilde = n,
    t_stat = t_stat,
    smd_ci = smd_ci
  ))
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

hedge_J <- function(df) {
  # compute unbiasing factor; works for small or large df;
  # thanks to Robert Calin-Jageman
  # Source: https://doi.org/10.31234/osf.io/s2597
  exp ( lgamma(df/2) - log(sqrt(df/2)) - lgamma((df-1)/2) )
  }

poolSD = function(x,y){
  n1 = length(x)
  n2 = length(y)
  sd1 = sd(x)
  sd2 = sd(y)
  sqrt((((n1 - 1)*(sd1^2)) + (n2 - 1)*(sd2^2))/((n1+n2)-2)) #calculate sd pooled
}

get_ncp_t <- function(t, df_error, alpha = 0.05) {
  # # Note: these aren't actually needed - all t related functions would fail earlier
  # if (!is.finite(t) || !is.finite(df_error)) {
  #   return(c(NA, NA))
  # }

  probs <- c(alpha / 2, 1 - alpha / 2)

  ncp <- suppressWarnings(optim(
    par = 1.1 * rep(t, 2),
    fn = function(x) {
      p <- pt(q = t, df = df_error, ncp = x)

      abs(max(p) - probs[2]) +
        abs(min(p) - probs[1])
    },
    control = list(abstol = 1e-09)
  ))
  t_ncp <- unname(sort(ncp$par))

  return(t_ncp)
}

get_ncp_t2 = function (ncp, df, conf.level = 0.95,
                       t.value, tol = 1e-09)
{
  sup.int.warns = TRUE
  alpha.lower = NULL
  alpha.upper = NULL
  if (missing(ncp)) {
    if (missing(t.value))
      stop("You need to specify either 'ncp' or its alias, 't.value,' you have not specified either")
    ncp <- t.value
  }
  if (df <= 0)
    stop("The degrees of freedom must be some positive value.",
         call. = FALSE)
  #if (abs(ncp) > 37.62)
  #  print("The observed noncentrality parameter of the noncentral t-distribution has exceeded 37.62 in magnitude (R's limitation for accurate probabilities from the noncentral t-distribution) in the function's iterative search for the appropriate value(s). The results may be fine, but they might be inaccurate; use caution.")
  if (sup.int.warns == TRUE)
    Orig.warn <- options()$warn
  options(warn = -1)
  if (!is.null(conf.level) & is.null(alpha.lower) & !is.null(alpha.upper))
    stop("You must choose either to use 'conf.level' or define the 'lower.alpha' and 'upper.alpha' values; here, 'upper.alpha' is specified but 'lower.alpha' is not",
         call. = FALSE)
  if (!is.null(conf.level) & !is.null(alpha.lower) & is.null(alpha.upper))
    stop("You must choose either to use 'conf.level' or define the 'lower.alpha' and 'upper.alpha' values; here, 'lower.alpha' is specified but 'upper.alpha' is not",
         call. = FALSE)
  if (!is.null(conf.level) & is.null(alpha.lower) & is.null(alpha.upper)) {
    alpha.lower <- (1 - conf.level)/2
    alpha.upper <- (1 - conf.level)/2
  }
  .conf.limits.nct.M1 <- function(ncp, df, conf.level = NULL,
                                  alpha.lower, alpha.upper, tol = 1e-09, sup.int.warns = TRUE,
                                  ...) {
    if (sup.int.warns == TRUE)
      Orig.warn <- options()$warn
    options(warn = -1)
    min.ncp = min(-150, -5 * ncp)
    max.ncp = max(150, 5 * ncp)
    .ci.nct.lower <- function(val.of.interest, ...) {
      (qt(p = alpha.lower, df = df, ncp = val.of.interest,
          lower.tail = FALSE, log.p = FALSE) - ncp)^2
    }
    .ci.nct.upper <- function(val.of.interest, ...) {
      (qt(p = alpha.upper, df = df, ncp = val.of.interest,
          lower.tail = TRUE, log.p = FALSE) - ncp)^2
    }
    if (alpha.lower != 0) {
      if (sup.int.warns == TRUE)
        Low.Lim <- suppressWarnings(optimize(f = .ci.nct.lower,
                                             interval = c(min.ncp, max.ncp), alpha.lower = alpha.lower,
                                             df = df, ncp = ncp, maximize = FALSE, tol = tol))
      if (sup.int.warns == FALSE)
        Low.Lim <- optimize(f = .ci.nct.lower, interval = c(min.ncp,
                                                            max.ncp), alpha.lower = alpha.lower, df = df,
                            ncp = ncp, maximize = FALSE, tol = tol)
    }
    if (alpha.upper != 0) {
      if (sup.int.warns == TRUE)
        Up.Lim <- suppressWarnings(optimize(f = .ci.nct.upper,
                                            interval = c(min.ncp, max.ncp), alpha.upper = alpha.upper,
                                            df = df, ncp = ncp, maximize = FALSE, tol = tol))
      if (sup.int.warns == FALSE)
        Up.Lim <- optimize(f = .ci.nct.upper, interval = c(min.ncp,
                                                           max.ncp), alpha.upper = alpha.upper, df = df,
                           ncp = ncp, maximize = FALSE, tol = tol)
    }
    if (alpha.lower == 0)
      Result <- list(Lower.Limit = -Inf, Prob.Less.Lower = 0,
                     Upper.Limit = Up.Lim$minimum, Prob.Greater.Upper = pt(q = ncp,
                                                                           ncp = Up.Lim$minimum, df = df))
    if (alpha.upper == 0)
      Result <- list(Lower.Limit = Low.Lim$minimum, Prob.Less.Lower = pt(q = ncp,
                                                                         ncp = Low.Lim$minimum, df = df, lower.tail = FALSE),
                     Upper.Limit = Inf, Prob.Greater.Upper = 0)
    if (alpha.lower != 0 & alpha.upper != 0)
      Result <- list(Lower.Limit = Low.Lim$minimum, Prob.Less.Lower = pt(q = ncp,
                                                                         ncp = Low.Lim$minimum, df = df, lower.tail = FALSE),
                     Upper.Limit = Up.Lim$minimum, Prob.Greater.Upper = pt(q = ncp,
                                                                           ncp = Up.Lim$minimum, df = df))
    if (sup.int.warns == TRUE)
      options(warn = Orig.warn)
    return(Result)
  }
  .conf.limits.nct.M2 <- function(ncp, df, conf.level = NULL,
                                  alpha.lower, alpha.upper, tol = 1e-09, sup.int.warns = TRUE,
                                  ...) {
    .ci.nct.lower <- function(val.of.interest, ...) {
      (qt(p = alpha.lower, df = df, ncp = val.of.interest,
          lower.tail = FALSE, log.p = FALSE) - ncp)^2
    }
    .ci.nct.upper <- function(val.of.interest, ...) {
      (qt(p = alpha.upper, df = df, ncp = val.of.interest,
          lower.tail = TRUE, log.p = FALSE) - ncp)^2
    }
    if (sup.int.warns == TRUE) {
      Low.Lim <- suppressWarnings(nlm(f = .ci.nct.lower,
                                      p = ncp, ...))
      Up.Lim <- suppressWarnings(nlm(f = .ci.nct.upper,
                                     p = ncp, ...))
    }
    if (sup.int.warns == FALSE) {
      Low.Lim <- nlm(f = .ci.nct.lower, p = ncp, ...)
      Up.Lim <- nlm(f = .ci.nct.upper, p = ncp, ...)
    }
    if (alpha.lower == 0)
      Result <- list(Lower.Limit = -Inf, Prob.Less.Lower = 0,
                     Upper.Limit = Up.Lim$estimate, Prob.Greater.Upper = pt(q = ncp,
                                                                            ncp = Up.Lim$estimate, df = df))
    if (alpha.upper == 0)
      Result <- list(Lower.Limit = Low.Lim$estimate, Prob.Less.Lower = pt(q = ncp,
                                                                          ncp = Low.Lim$estimate, df = df, lower.tail = FALSE),
                     Upper.Limit = Inf, Prob.Greater.Upper = 0)
    if (alpha.lower != 0 & alpha.upper != 0)
      Result <- list(Lower.Limit = Low.Lim$estimate, Prob.Less.Lower = pt(q = ncp,
                                                                          ncp = Low.Lim$estimate, df = df, lower.tail = FALSE),
                     Upper.Limit = Up.Lim$estimate, Prob.Greater.Upper = pt(q = ncp,
                                                                            ncp = Up.Lim$estimate, df = df))
    return(Result)
  }
  Res.M1 <- Res.M2 <- NULL
  try(Res.M1 <- .conf.limits.nct.M1(ncp = ncp, df = df, conf.level = NULL,
                                    alpha.lower = alpha.lower, alpha.upper = alpha.upper,
                                    tol = tol, sup.int.warns = sup.int.warns), silent = TRUE)
  if (length(Res.M1) != 4)
    Res.M1 <- NULL
  try(Res.M2 <- .conf.limits.nct.M2(ncp = ncp, df = df, conf.level = NULL,
                                    alpha.lower = alpha.lower, alpha.upper = alpha.upper,
                                    tol = tol, sup.int.warns = sup.int.warns), silent = TRUE)
  if (length(Res.M2) != 4)
    Res.M2 <- NULL
  Low.M1 <- Res.M1$Lower.Limit
  Prob.Low.M1 <- Res.M1$Prob.Less.Lower
  Upper.M1 <- Res.M1$Upper.Limit
  Prob.Upper.M1 <- Res.M1$Prob.Greater.Upper
  Low.M2 <- Res.M2$Lower.Limit
  Prob.Low.M2 <- Res.M2$Prob.Less.Lower
  Upper.M2 <- Res.M2$Upper.Limit
  Prob.Upper.M2 <- Res.M2$Prob.Greater.Upper
  Min.for.Best.Low <- min((c(Prob.Low.M1, Prob.Low.M2) - alpha.lower)^2)
  if (!is.null(Res.M1)) {
    if (Min.for.Best.Low == (Prob.Low.M1 - alpha.lower)^2)
      Best.Low <- 1
  }
  if (!is.null(Res.M2)) {
    if (Min.for.Best.Low == (Prob.Low.M2 - alpha.lower)^2)
      Best.Low <- 2
  }
  Min.for.Best.Up <- min((c(Prob.Upper.M1, Prob.Upper.M2) -
                            alpha.upper)^2)
  if (!is.null(Res.M1)) {
    if (Min.for.Best.Up == (Prob.Upper.M1 - alpha.upper)^2)
      Best.Up <- 1
  }
  if (!is.null(Res.M2)) {
    if (Min.for.Best.Up == (Prob.Upper.M2 - alpha.upper)^2)
      Best.Up <- 2
  }
  if (is.null(Res.M1)) {
    Low.M1 <- NA
    Prob.Low.M1 <- NA
    Upper.M1 <- NA
    Prob.Upper.M1 <- NA
  }
  if (is.null(Res.M2)) {
    Low.M2 <- NA
    Prob.Low.M2 <- NA
    Upper.M2 <- NA
    Prob.Upper.M2 <- NA
  }
  Result <- c(c(Low.M1, Low.M2)[Best.Low],
              c(Upper.M1, Upper.M2)[Best.Up])
  return(Result)
}

utils::globalVariables(c("sd1"))
