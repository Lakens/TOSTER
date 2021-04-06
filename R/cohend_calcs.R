

d_est_pair <- function(n,
                       m1,
                       m2,
                       sd1,
                       sd2,
                       r12,
                       type = "g",
                       denom = "z",
                       alpha = .05){

  #n <- nrow(data)
  #i1 <- data$i1
  #i2 <- data$i2
  #m1 <- base::mean(i1)
  #m2 <- base::mean(i2)
  #sd1  <- stats::sd(i1)
  #sd2  <- stats::sd(i2)
  #r12 <- stats::cor(i1, i2)

  sdif <- sqrt(sd1 ^ 2 + sd2 ^ 2 - 2 * r12 * sd1 * sd2)
  if(denom == "z"){
    d_denom = sdif
  } else if(denom == "rm"){
    d_denom = sdif * sqrt(2*(1-r12))
  }


  #res <- t.test(
  #  y = data$i1,
  #  x = data$i2,
  #  paired = TRUE,
  #  mu = 0,
  #  alternative = "two.sided"
  #)

  #t <- unname(res$statistic)

  df <- n-1

  cohend = abs(m2-m1) / d_denom
  #if(denom == "rm"){
  #  cohend = cohend/sqrt(2*(1-r12)) # drm
  #}
  cohen_df = 2*(n)-2
  #J <- gamma(df / 2) / (sqrt(df / 2) * gamma((df - 1) / 2))
  J = gamma(cohen_df / 2) / (sqrt(cohen_df / 2) * gamma((cohen_df - 1) / 2))

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
  #d_sigma = sqrt((df + 1) / (df - 1) * ((2*(1-r12))/n) * (1 + cohend ^ 2 /8))
  d_sigma = sqrt((cohen_df/(cohen_df-2)) * ((2*(1-r12))/n)*(1+cohend^2*(n/(2*(1-r12)))) - cohend^2/J^2)
  if(denom == "z"){
    d_sigma = d_sigma * sqrt(2*(1-r12))
  }
  # Get cohend confidence interval
  tlow <- suppressWarnings(qt(alpha,
                              df = cohen_df,
                              ncp = d_lambda))
  thigh <- suppressWarnings(qt(alpha,
                               df = cohen_df,
                               ncp = d_lambda))
  if (m1 == m2) {
    dlow <- (tlow * (2 * n)) / (sqrt(cohen_df) * sqrt(2 * n))
    dhigh <- (thigh * (2 * n)) / (sqrt(cohen_df) * sqrt(2 * n))
  } else {
    dlow <- tlow / d_lambda * cohend
    dhigh <- thigh / d_lambda * cohend
  }

  if (unname(res$estimate) < 0) {
    cohend <- cohend * -1
    tdlow <- dlow
    dlow <- dhigh * -1
    dhigh <- tdlow * -1
    d_lambda <- cohend * sqrt(n / (2*(1 - r12)))
  }

  return(list(
    cohend = cohend,
    cohen_df = cohen_df,
    dlow = dlow,
    dhigh = dhigh,
    d_sigma = d_sigma,
    d_lambda = d_lambda,
    smd_label = smd_label,
    J = J
  ))
}


d_est_ind <- function(data,
                      type = "g",
                      var.equal = TRUE,
                      alpha = .05){

  dataTTest = data

  v <- tapply(dataTTest$dep, dataTTest$group, function(x) jmvcore::tryNaN(var(x)))
  n <- tapply(dataTTest$dep, dataTTest$group, length)
  m <- tapply(dataTTest$dep, dataTTest$group, function(x) jmvcore::tryNaN(mean(x)))
  med <- tapply(dataTTest$dep, dataTTest$group, function(x) jmvcore::tryNaN(median(x)))
  se <- sqrt(v/n)
  sd <- sqrt(v)

  #tresult <- t.test(dep ~ group,
  #                  dataTTest,
  #                  paired = FALSE,
  #                  alternative = "two.sided",
  #                  mu = 0,
  #                  var.equal=var.equal)

  #t <- unname(tresult$statistic)
  #p <- unname(tresult$p.value)
  #df <- unname(tresult$parameter)

  if (var.equal) {
    denomSD <- sqrt((((n[1] - 1)*(sd[1]^2)) + (n[2] - 1)*(sd[2]^2))/((n[1]+n[2])-2)) #calculate sd pooled
    cohen_df = n[1] + n[1] - 2
  } else {
    denomSD <- sqrt((sd[1]^2 + sd[2]^2)/2) #calculate sd root mean squared for Welch's t-test
    cohen_df1 = (n[1] - 1)*(n[2] - 1)*(sd[1]^2+sd[2]^2)^2
    cohen_df2 = (n[2]-1)*sd[1]^4+(n[1]-1)*sd[2]^4
    cohen_df = cohen_df1/cohen_df2
  }

  denomSD[is.na(denomSD)] <- NaN

  #denomSD <- jmvcore::tryNaN(sqrt(((n[1]-1)*v[1]+(n[2]-1)*v[2])/(n[1]+n[2]-2)))
  d <- abs(m[1]-m[2])/denomSD # Cohen's d

  d[is.na(d)] <- NaN

  cohend = d
  ntilde <- harmonic.mean(c(n[1],n[2]))

  # Compute unbiased Hedges' g
  # Use the lgamma function, and update to what Goulet-Pelletier & Cousineau used; works with larger inputs
  J <- gamma(cohen_df/2)/(sqrt(cohen_df/2)*gamma((cohen_df-1)/2))

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
  d_sigma = sqrt((n[1]+n[2])/(n[1]*n[2])+(cohend^2/(2*(n[1]+n[2]))))

  # Confidence interval of the SMD from Goulet-Pelletier & Cousineau
  tlow <- qt(1 / 2 - (1-alpha*2) / 2, df = cohen_df, ncp = d_lambda)
  thigh <- qt(1 / 2 + (1-alpha*2) / 2, df = cohen_df, ncp = d_lambda)
  if(m[1] == m[2]) {
    dlow <- (tlow*(n[1]+n[2])) / (sqrt(cohen_df) * sqrt(n[1]*n[2]))
    dhigh <- (thigh*(n[1]+n[2])) / (sqrt(cohen_df) * sqrt(n[1]*n[2]))
  } else {
    dlow <- tlow / d_lambda * cohend
    dhigh <- thigh / d_lambda * cohend
  }

  # The function provided by Goulet-Pelletier & Cousineau works with +mdiff.
  # Now we fix the signs and directions for negative differences
  if ((m[1]-m[2]) < 0) {
    cohend <- cohend * -1
    tdlow <- dlow
    dlow <- dhigh * -1
    dhigh <- tdlow * -1
    d_lambda <- cohend * sqrt(ntilde/2)
  }


  return(list(
    cohend = cohend,
    cohen_df = cohen_df,
    dlow = dlow,
    dhigh = dhigh,
    d_sigma = d_sigma,
    d_lambda = d_lambda,
    smd_label = smd_label,
    J = J
  ))
}


d_est_one <- function(t_result,
                      sd,
                      type = "g",
                      alpha = .05){

  cohend <- abs(t_result$estimate)/sd # Cohen's d
  df <- t_result$parameter
  n = df+1
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
  tlow <- suppressWarnings(qt(alpha,
                              df = df,
                              ncp = d_lambda))
  thigh <- suppressWarnings(qt(alpha,
                               df = df,
                               ncp = d_lambda))
  if(dif == 0) {
    dlow <- (tlow*(2*n)) / (sqrt(df) * sqrt(2*n))
    dhigh <- (thigh*(2*n)) / (sqrt(df) * sqrt(2*n))
  } else {
    dlow <- tlow / d_lambda * cohend
    dhigh <- thigh / d_lambda * cohend
  }

  if ((dif) < 0) {
    cohend <- cohend * -1
    tdlow <- dlow
    dlow <- dhigh * -1
    dhigh <- tdlow * -1
  }



  return(list(
    cohend = cohend,
    cohen_df = cohen_df,
    dlow = dlow,
    dhigh = dhigh,
    d_sigma = d_sigma,
    d_lambda = d_lambda,
    smd_label = smd_label,
    J = J
  ))
}
