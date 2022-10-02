
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

  cohend = abs(m1-m2) / d_denom

  d_df = 2*(n)-2
  #J <- gamma(df / 2) / (sqrt(df / 2) * gamma((df - 1) / 2))
  if(type = "g"){
    J = hedge_J(d_df)
  } else {
    J = 1
  }


  if(denom == "z"){
    if (type == 'g') {
      cohend <-  cohend * J
      smd_label = "Hedges's g(z)"
    } else {
      smd_label = "Cohen's d(z)"
    }
  } else if(denom == "rm"){
    if (type == 'g') {
      cohend <-  cohend * J
      smd_label = "Hedges's g(rm)"
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

  if (m1-m2 < 0) {
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
  } else if (denom == "glass2"){
    denomSD <- sd2
    hn <- 1 / n2 + denomSD^2 / (n1 * denomSD^2)
    d_df = n2 - 1
  }


  denomSD[is.na(denomSD)] <- NaN

  #denomSD <- jmvcore::tryNaN(sqrt(((n1-1)*v[1]+(n2-1)*v[2])/(n1+n2-2)))
  d <- abs(m1-m2)/denomSD # Cohen's d

  #d[is.na(d)] <- NaN

  cohend = d
  ntilde <- harm_mean(n1,n2)

  # Compute unbiased Hedges's g
  # Use the lgamma function, and update to what Goulet-Pelletier & Cousineau used; works with larger inputs
  if(type = "g"){
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

  # The function provided by Goulet-Pelletier & Cousineau works with +mdiff.
  # Now we fix the signs and directions for negative differences

  }

  if(smd_ci == "nct"){
    if(var.equal == TRUE && !(denom %in% c("glass1","glass2"))){
      d_sigma = denomSD * sqrt(1 / n1 + 1 / n2)
      t_stat = cohend/d_sigma
      ts <- get_ncp_t(t_stat, d_df, alpha)
    } else if(denom %in% c("glass1","glass2")){
      d_sigma = (sd2 * sqrt(1 / n2 + sd1^2 / (n1 * sd2^2)))
      t_stat = cohend/d_sigma
      ts <- get_ncp_t(t_stat, d_df, alpha)
    } else {
      se1 <- sqrt(sd1^2 / n1)
      se2 <- sqrt(sd2^2 / n2)
      d_sigma <- sqrt(se1^2 + se2^2)
      t_stat = cohend/d_sigma
      ts <- get_ncp_t(t_stat, d_df, alpha)
    }

    d_low <- ts[1] * sqrt(hn)
    d_high <- ts[2] * sqrt(hn)
  }

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
    hn = hn,
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
                      alpha = .05,
                      smd_ci = "goulet"){

  cohend <- abs(mu-testValue)/sd # Cohen's d
  df <- n-1
  d_df = df
  hn <- 1 / n
  # Compute unbiased Hedges' g
  # Use the lgamma function, and update to what Goulet-Pelletier & cousineau used; works with larger inputs
  if(type = "g"){
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
