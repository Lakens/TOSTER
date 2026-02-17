
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
                       smd_ci = "goulet",
                       tr = 0){

  # Trimming adjustments
  cg <- trim_rescale(tr)
  h <- trim_h(n, tr)

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

  if (tr > 0) {
    df <- h - 1
  } else {
    df <- n - 1
  }
  hn <- 1 / n
  cohend = cg * abs(m1-m2) / d_denom
  if(smd_ci == "goulet"){
    d_df = 2*(n)-2
  } else {
    if (tr > 0) {
      d_df = h - 1
    } else {
      d_df = n - 1
    }
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
      smd_label = "Hedges's g(z)"
    } else {
      smd_label = "Cohen's d(z)"
    }
  } else if(denom == "rm"){
    if (type == 'g') {
      smd_label = "Hedges's g(rm)"
    } else {
      smd_label = "Cohen's d(rm)"
    }
  } else if(denom %in% c("glass1","glass2")){
    if (type == 'g') {
      smd_label = "Glass's delta(g)"
    } else {
      smd_label = "Glass's delta(d)"
    }
  }

  if (tr > 0) {
    # For trimmed paired, use h-based SE and noncentrality
    if (denom == "z") {
      # Yuen-adjusted SE for differences
      SE_yuen <- d_denom / sqrt(h) * sqrt((n - 1) / (h - 1))
      d_lambda <- abs(m1 - m2) / SE_yuen
      d_unscaled <- cohend / (cg * J)
      d_sigma <- cg * sqrt(1 / h * ((n - 1) / (h - 1))) *
        sqrt(d_unscaled^2 / (2 * (h - 1)) + 1)
      d_sigma <- d_sigma * J
    } else if (denom %in% c("glass1", "glass2")) {
      SE_yuen <- d_denom / sqrt(h) * sqrt((n - 1) / (h - 1))
      d_lambda <- abs(m1 - m2) / SE_yuen
      d_unscaled <- cohend / (cg * J)
      d_s1 <- sdif^2 / (d_denom^2 * (h - 1)) + d_unscaled^2 / (2 * (h - 1))
      d_sigma <- cg * J * sqrt(d_s1)
    }
  } else {
    d_lambda <- cohend * sqrt(n / (2*(1 - r12)))

    # Equation 4b Goulet-Pelletier and Cousineau, 2018
    d_sigma = sqrt((d_df/(d_df-2)) * ((2*(1-r12))/n)*(1+cohend^2*(n/(2*(1-r12)))) - cohend^2/J^2)
    if(denom == "z"){
      d_sigma = d_sigma * sqrt(2*(1-r12))
    }
    if(smd_ci == "nct" && denom != "rm"){
      d_sigma2 =  1/n + (1 - (d_df-2)/(d_df*J^2)) * cohend^2
      d_sigma = sqrt(d_sigma2)
    }
    if(denom %in% c("glass1","glass2")){
      d_s1 = sdif^2/(d_denom^2*(df)) + cohend^2 / (2*(df))
      d_sigma = sqrt(d_s1)
    }
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
    if (tr > 0) {
      if (denom == "z") {
        SE1 <- d_denom / sqrt(h) * sqrt((n - 1) / (h - 1))
      } else {
        SE1 <- d_denom / sqrt(h) * sqrt((n - 1) / (h - 1))
      }
      t_stat <- abs(m1 - m2) / SE1
      ts <- get_ncp_t2(t_stat, d_df, conf.level = 1-alpha*2)
      conv_factor <- cg * (1 / sqrt(h)) * sqrt((n - 1) / (h - 1))
      dlow <- ts[1] * conv_factor * J
      dhigh <- ts[2] * conv_factor * J
    } else {
      if(denom == "rm"){
        t_stat <- (abs(m1 - m2) / (sqrt((sd1^2 + sd2^2)-(2*r12*sd1*sd2))/sqrt(n))) * sqrt(2*(1-r12))
      }else{
        SE1 <- d_denom / sqrt(n)
        t_stat = abs(m1 - m2) / SE1
      }

      ts <- get_ncp_t2(t_stat, d_df, conf.level = 1-alpha*2)
      dlow <- ts[1] * sqrt(hn) * J
      dhigh <- ts[2] * sqrt(hn) * J
    }
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
    if (tr > 0) {
      d_lambda <- -d_lambda
    } else {
      d_lambda <- cohend * sqrt(n / (2*(1 - r12)))
    }
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
                      smd_ci = "goulet",
                      tr = 0){

  # Trimming adjustments
  cg <- trim_rescale(tr)
  h1 <- trim_h(n1, tr)
  h2 <- trim_h(n2, tr)
  N <- n1 + n2

  if (var.equal) {
    denomSD <- sqrt((((n1 - 1)*(sd1^2)) + (n2 - 1)*(sd2^2))/((n1+n2)-2)) #calculate sd pooled
    if (tr > 0) {
      d_df <- h1 + h2 - 2
      hn <- (1 / h1 + 1 / h2) * ((N - 2) / (h1 + h2 - 2))
    } else {
      d_df = n1 + n2 - 2
      hn <- (1 / n1 + 1 / n2)
    }
  } else {
    denomSD <- sqrt((sd1^2 + sd2^2)/2) #calculate sd root mean squared for Welch's t-test
    if (tr > 0) {
      # Welch-Satterthwaite with Winsorized variances and h
      d_df1 = (sd1^2/h1 + sd2^2/h2)^2
      d_df2 = (sd1^2/h1)^2/(h1-1) + (sd2^2/h2)^2/(h2-1)
      d_df = d_df1/d_df2
      hn <- (2 * (h2 * sd1^2 + h1 * sd2^2)) / (h1 * h2 * (sd1^2 + sd2^2))
    } else {
      d_df1 = (n1 - 1)*(n2 - 1)*(sd1^2+sd2^2)^2
      d_df2 = (n2-1)*sd1^4+(n1-1)*sd2^4
      d_df = d_df1/d_df2
      hn <- (2 * (n2 * sd1^2 + n1 * sd2^2)) / (n1 * n2 * (sd1^2 + sd2^2))
    }
  }

  if (denom == "glass1"){
    denomSD <- sd1
    if (tr > 0) {
      d_df = h1 - 1
    } else {
      d_df = n1 - 1
    }
    hn <- 1 / n2 + denomSD^2 / (n1 * denomSD^2)
    n_glass = n1
    nn_glass = n2
    sdn_glass = sd2
  } else if (denom == "glass2"){
    denomSD <- sd2
    hn <- 1 / n2 + denomSD^2 / (n1 * denomSD^2)
    if (tr > 0) {
      d_df = h2 - 1
    } else {
      d_df = n2 - 1
    }
    n_glass = n2
    nn_glass = n1
    sdn_glass = sd1
  }


  denomSD[is.na(denomSD)] <- NaN

  d <- cg * abs(m1-m2)/denomSD # Cohen's d (or robust version)

  cohend = d
  ntilde <- harm_mean(n1,n2)

  # Compute unbiased Hedges's g
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


  if (tr > 0) {
    # Trimmed SE and noncentrality
    d_unscaled <- cohend / (cg * J)
    if (var.equal && !(denom %in% c("glass1", "glass2"))) {
      # Pooled Winsorized: Yuen-Dixon adjusted variance
      S_tilde2 <- (N - 2) * denomSD^2 / (h1 + h2 - 2)
      S_tilde <- sqrt(S_tilde2)
      SE1 <- S_tilde * sqrt(1/h1 + 1/h2)
      d_lambda <- abs(m1 - m2) / SE1
      conv <- cg * sqrt((h1 + h2) * (N - 2) / (h1 * h2 * (h1 + h2 - 2)))
      d_sigma <- conv * sqrt(d_unscaled^2 / (2 * (h1 + h2 - 2)) + (h1 + h2) / (h1 * h2))
      d_sigma <- d_sigma * J
    } else if (denom %in% c("glass1", "glass2")) {
      h_glass <- if (denom == "glass1") h1 else h2
      SE1 <- denomSD / sqrt(h_glass) * sqrt((n_glass - 1) / (h_glass - 1))
      d_lambda <- abs(m1 - m2) / SE1
      d_sigma <- cg * J * sqrt((sdn_glass^2/denomSD^2)/(nn_glass-1) + 1/(h_glass-1) + d_unscaled^2/(2*(h_glass-1)))
    } else {
      # avg (Welch) with trimming
      se1 <- sqrt(sd1^2 / h1)
      se2 <- sqrt(sd2^2 / h2)
      SE1 <- sqrt(se1^2 + se2^2)
      d_lambda <- abs(m1 - m2) / SE1
      d_sigma2 <- cg^2 * (d_unscaled^2 * (sd1^4 / (h1-1) + sd2^4 / (h2-1)) / (8*denomSD^4) +
        (sd1^2 / (h1-1) + sd2^2 / (h2-1)) / denomSD^2)
      d_sigma <- sqrt(d_sigma2) * J
    }
  } else {
    if(var.equal == TRUE && !(denom %in% c("glass1","glass2"))){
      mult_lamb = sqrt((n1*n2*(sd1^2 + sd2^2))/(2*(n2*sd1^2 + n1*sd2^2)))
      d_lambda = cohend * mult_lamb
    } else if(denom %in% c("glass1","glass2")){
      d_lambda <- cohend * sqrt(ntilde/2)
    } else {
      d_lambda <- cohend * sqrt(ntilde/2)
    }
  }

  # add options for cohend here
  if(smd_ci == "goulet"){
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
  } else if (tr == 0) {
    if (denom %in% c("glass1", "glass2")) {
      N = n1 + n2
      d_sigma2 = (sdn_glass^2/denomSD^2)/(nn_glass-1) + 1/(n_glass-1) + cohend^2/(2*(n_glass-1))
      d_sigma = sqrt(d_sigma2)
    } else {
      if (var.equal) {
        d_sigma2 =  1/n1 + 1/n2 + (1 - (d_df-2)/(d_df*J^2)) * cohend^2
        d_sigma = sqrt(d_sigma2)
      } else{
        d_sigma2 = cohend^2 * (sd1^4 / (n1-1) + sd2^4 / (n2-1)) / (8*denomSD^4) +
          (sd1^2 / (n1-1) + sd2^2 / (n2-1)) / denomSD^2
        d_sigma = sqrt(d_sigma2)
      }
    }
  }

  if(smd_ci == "nct"){
    if (tr > 0) {
      # SE1 and d_lambda already computed in trimming block above
      t_stat <- abs(m1 - m2) / SE1
      ts <- get_ncp_t2(t_stat, d_df, conf.level = 1-alpha*2)
      if (var.equal && !(denom %in% c("glass1", "glass2"))) {
        conv <- cg * sqrt((h1 + h2) * (N - 2) / (h1 * h2 * (h1 + h2 - 2)))
        dlow <- ts[1] * conv * J
        dhigh <- ts[2] * conv * J
      } else if (denom %in% c("glass1", "glass2")) {
        h_glass <- if (denom == "glass1") h1 else h2
        conv <- cg * (1 / sqrt(h_glass)) * sqrt((n_glass - 1) / (h_glass - 1))
        dlow <- ts[1] * conv * J
        dhigh <- ts[2] * conv * J
      } else {
        # Welch avg
        conv <- cg * SE1 / denomSD
        dlow <- ts[1] * conv * J
        dhigh <- ts[2] * conv * J
      }
    } else {
      if( !(denom %in% c("glass1","glass2"))){
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
    }
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
    if (tr == 0) {
      if(var.equal == TRUE && !(denom %in% c("glass1","glass2"))){
        mult_lamb = sqrt((n1*n2*(sd1^2 + sd2^2))/(2*(n2*sd1^2 + n1*sd2^2)))
        d_lambda = cohend * mult_lamb
      } else if(denom %in% c("glass1","glass2")){
        d_lambda <- cohend * sqrt(ntilde/2)
      } else {
        d_lambda <- cohend * sqrt(ntilde/2)
      }
    } else {
      d_lambda <- -d_lambda
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
                      smd_ci = "goulet",
                      tr = 0){

  # Trimming adjustments
  cg <- trim_rescale(tr)
  h <- trim_h(n, tr)

  cohend <- cg * abs(mu-testValue)/sd # Cohen's d (or robust version)
  if (tr > 0) {
    df <- h - 1
  } else {
    df <- n - 1
  }
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

  if (tr > 0) {
    # Noncentrality: lambda = (trimmed mean diff / Winsorized SE)
    # SE_trimmed = sd / sqrt(h) * sqrt((n-1)/(h-1)) for Yuen adjustment
    SE_yuen <- sd / sqrt(h) * sqrt((n - 1) / (h - 1))
    d_lambda <- abs(mu - testValue) / SE_yuen
    # SE of d_R: uses the noncentral-t approximation with trimming
    d_unscaled <- cohend / (cg * J)  # remove cg and J to get raw ratio
    d_sigma <- cg * sqrt(1 / h * ((n - 1) / (h - 1))) *
      sqrt(d_unscaled^2 / (2 * (h - 1)) + 1)
    d_sigma <- d_sigma * J
  } else {
    d_lambda <- cohend * sqrt(n)
    if(smd_ci == "goulet"){
      d_sigma = sqrt((df/(df-2)) * (1/n) *(1+cohend^2*(n/1)) - cohend^2/J^2)
    } else {
      d_sigma = sqrt(1/n + (cohend^2/(2*n)))
    }
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
    if (tr > 0) {
      SE1 <- sd / sqrt(h) * sqrt((n - 1) / (h - 1))
      t_stat <- abs(mu - testValue) / SE1
      ts <- get_ncp_t2(t_stat, d_df, conf.level = 1-alpha*2)
      # Convert noncentrality to delta_R
      conv_factor <- cg * (1 / sqrt(h)) * sqrt((n - 1) / (h - 1))
      dlow <- ts[1] * conv_factor * J
      dhigh <- ts[2] * conv_factor * J
    } else {
      SE1 <- sd / sqrt(n)
      t_stat = abs(mu-testValue) / SE1
      ts <- get_ncp_t2(t_stat, d_df, conf.level = 1-alpha*2)
      dlow <- ts[1] * sqrt(hn)*J
      dhigh <- ts[2] * sqrt(hn)*J
    }
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

# Internal helper: Winsorized variance
# Replaces the g smallest and g largest observations with the nearest remaining
# values, then computes the variance of the Winsorized sample.
winvar <- function(x, tr = 0.2) {
  n <- length(x)
  g <- floor(tr * n)
  if (g == 0) return(var(x))
  y <- sort(x)
  y[seq_len(g)] <- y[g + 1L]
  y[seq.int(n - g + 1L, n)] <- y[n - g]
  var(y)
}

# Internal helper: Rescaling constant for trimmed SMD
# Returns c(gamma) such that c(gamma) * (trimmed mean diff / Winsorized SD)
# equals Cohen's delta under normality.
trim_rescale <- function(tr) {
  if (tr == 0) return(1)
  a <- qnorm(1 - tr)
  sqrt(1 - 2 * tr + 2 * tr * a^2 - 2 * a * dnorm(a))
}

# Internal helper: Effective sample size after trimming
# h = n - 2 * floor(tr * n)
trim_h <- function(n, tr) {
  n - 2L * floor(tr * n)
}

# Internal helper to resolve denom into concrete arguments
# Returns a list with resolved values of: glass, rm_correction, var.equal
# and a character vector of messages to emit (if any)
resolve_denom <- function(denom,
                          sample_type,
                          var.equal,
                          rm_correction,
                          glass,
                          var.equal_explicit = FALSE,
                          rm_correction_explicit = FALSE,
                          glass_explicit = FALSE) {

  msgs <- character(0)

  if (denom == "auto") {
    # Return current values unchanged, no messages
    return(list(
      var.equal = var.equal,
      rm_correction = rm_correction,
      glass = glass,
      messages = msgs
    ))
  }

  # --- Design-validity checks ---
  paired_only <- c("rm")
  paired_or_one <- c("z")
  ind_only <- c("pooled", "avg")
  needs_two_samples <- c("glass1", "glass2", "rm")

  if (denom %in% ind_only && sample_type != "Two Sample") {
    stop(paste0("denom = '", denom, "' is only valid for independent samples designs."))
  }
  if (denom %in% paired_only && sample_type != "Paired Sample") {
    stop(paste0("denom = '", denom, "' is only valid for paired samples designs."))
  }
  if (denom %in% paired_or_one && sample_type == "Two Sample") {
    stop(paste0("denom = '", denom, "' is not valid for independent samples designs."))
  }
  if (denom %in% needs_two_samples && sample_type == "One Sample") {
    stop(paste0("denom = '", denom, "' is not valid for one-sample designs."))
  }

  # --- Remapping with conflict detection ---
  new_glass <- glass
  new_rm <- rm_correction
  new_var.equal <- var.equal

  if (denom == "z") {
    if (rm_correction_explicit && isTRUE(rm_correction)) {
      msgs <- c(msgs, "denom = 'z' overrides rm_correction to FALSE.")
    }
    if (glass_explicit && !is.null(glass)) {
      msgs <- c(msgs, "denom = 'z' overrides glass argument.")
    }
    new_rm <- FALSE
    new_glass <- NULL

  } else if (denom == "rm") {
    if (rm_correction_explicit && !isTRUE(rm_correction)) {
      msgs <- c(msgs, "denom = 'rm' overrides rm_correction to TRUE.")
    }
    if (glass_explicit && !is.null(glass)) {
      msgs <- c(msgs, "denom = 'rm' overrides glass argument.")
    }
    new_rm <- TRUE
    new_glass <- NULL

  } else if (denom == "pooled") {
    if (var.equal_explicit && !isTRUE(var.equal)) {
      msgs <- c(msgs, "denom = 'pooled' overrides var.equal to TRUE.")
    }
    if (glass_explicit && !is.null(glass)) {
      msgs <- c(msgs, "denom = 'pooled' overrides glass argument.")
    }
    new_var.equal <- TRUE
    new_glass <- NULL
    new_rm <- FALSE

  } else if (denom == "avg") {
    if (var.equal_explicit && isTRUE(var.equal)) {
      msgs <- c(msgs, "denom = 'avg' overrides var.equal to FALSE.")
    }
    if (glass_explicit && !is.null(glass)) {
      msgs <- c(msgs, "denom = 'avg' overrides glass argument.")
    }
    new_var.equal <- FALSE
    new_glass <- NULL
    new_rm <- FALSE

  } else if (denom %in% c("glass1", "glass2")) {
    if (glass_explicit && !is.null(glass) && glass != denom) {
      msgs <- c(msgs, paste0("denom = '", denom, "' overrides glass argument."))
    }
    if (rm_correction_explicit && isTRUE(rm_correction)) {
      msgs <- c(msgs, paste0("denom = '", denom, "' overrides rm_correction to FALSE."))
    }
    new_glass <- denom
    new_rm <- FALSE
  }

  return(list(
    var.equal = new_var.equal,
    rm_correction = new_rm,
    glass = new_glass,
    messages = msgs
  ))
}

utils::globalVariables(c("sd1"))
