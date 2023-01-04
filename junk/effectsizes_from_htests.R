#' @rdname htest-helpers
#' @export

smd_from_htest = function(htest,
                          alpha = .05,
                          bias_correction = FALSE){
  if(!("htest" %in% class(htest)) || !grepl("t-test",htest$method,
                                            ignore.case=TRUE)){
    stop("htest must be an htest from a t-test.")

  }
  smd_type = ifelse(bias_correction, "g", "d")

  if(grepl("two",htest$method,ignore.case=TRUE)){
    mu = unname(htest$estimate[1])- unname(htest$estimate[2])
  } else {
    mu = unname(htest$estimate)
  }

  if(grepl("paired",htest$method,ignore.case=TRUE) || grepl("one",htest$method,ignore.case=TRUE)){
    mult = 1
    sample_size = unname(htest$parameter) + 1
    denom_sd = sqrt(sample_size) * unname(htest$stderr)

    smd_vals = d_est_one(n = sample_size,
                         mu = mu,
                         sd = denom_sd,
                         testValue = 0,
                         type = smd_type,
                         alpha = alpha,
                         smd_ci = "nct")
    res = data.frame(row.names = smd_vals$smd_label,
                     estimate = smd_vals$d,
                     lower.ci = smd_vals$dlow,
                     upper.ci = smd_vals$dhigh,
                     conf.level = 1-alpha)

  }

  if(grepl("two",htest$method,ignore.case=TRUE)){
    mult = 2
    if(is.null(sample_size)){
      sample_size = unname(htest$parameter) + 2
      if(grepl("welch",htest$method,ignore.case=TRUE)){
        stop("sample_size argument must be provided if Welch's t-test provided.")
      }
    }
    if(bias_correction){
      smd_label = "Hedges's g"
      J = TOSTER:::hedge_J(unname(htest$parameter))
    } else {
      smd_label = "Cohen's d"
      J = 1
    }
    tstat = mu / unname(htest$stderr)
    ts = TOSTER:::get_ncp_t2(tstat,
                        unname(htest$parameter),
                        1-alpha)
    ci <- 2 * ts / sqrt(unname(htest$parameter))

    smd <- unname(mult * htest$statistic / sqrt(htest$parameter))

    res = data.frame(row.names = smd_label,
                     estimate = smd*J,
                     lower.ci = min(ci)*J,
                     upper.ci = max(ci)*J,
                     conf.level = 1-alpha)

  }

  return(res)
}

#' @rdname htest-helpers
#' @export

ses_from_htest = function(htest){
  if(!("htest" %in% class(htest)) || !grepl("wilcoxon",htest$method,
                                            ignore.case=TRUE)){
    stop("htest must be an htest from a Wilcoxon-Mann-Whitney test.")
  }
}


d_t_ind <- function(mu,
                    se,
                    sd,
                    df,
                    type = "g",
                    alpha = .05){
  #message("Estimations from htest may differ slightly from smd_calc function")
  if(type == "g"){
    J = hedge_J(d_df)
    smd_label = "Hedges's g"
  } else {
    J = 1
    smd_label = "Cohen's d"
  }
  d <- mu/sd # Cohen's d

  cohend = J*d


  d_df = df
  ntot = df+2
  n1 = ntot/2
  n2 = ntot/2
  sd1 = sd
  sd2 = sd1
  hn <- (1 / n1 + 1 / n2)
  t_stat = mu/se

  d_sigma = sqrt(((n1 + n2) / (n1 * n2) + d ^ 2 / (2 * (n1 + n2))) * J ^ 2)

  ts <- TOSTER:::get_ncp_t2(t_stat, d_df, conf.level = 1-alpha*2)

  dlow <- ts[1] * sqrt(hn) * J
  dhigh <- ts[2] * sqrt(hn) * J

  return(list(
    d = cohend,
    d_df = d_df,
    dlow = dlow,
    dhigh = dhigh,
    d_sigma = d_sigma,
    d_lambda = NA,
    #hn = hn,
    smd_label = smd_label,
    J = J,
    d_denom = sd,
    ntilde = n1,
    t_stat = t_stat,
    smd_ci = "nct"
  ))
}
