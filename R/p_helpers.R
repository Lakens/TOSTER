#' Alternative values to complement p-values
#'
#' Functions to help convert p-values that can be reported alongside test results.
#'
#' @param p.value The p-value. Must be numeric and between 0 and 1.
#' @param calc The calculative approach.
#' @param log_base The base for the logarithm. Default is log base 2.
#' @param bf Bayes factor.
#' @param df degrees of freedom.
#' @param num_df numerator degrees of freedom.
#' @param den_df denominator degrees of freedom.
#' @examples
#' # simple example with t-test
#' tres = t.test(extra ~ group, data = sleep)
#'
#' # As a data frame
#' s_value(tres$p.value)
#'
#' @name p_calibrate
#'


#' @family htest
#' @name p_calibrate
#' @export
#'

s_value = function(p.value,
                   log_base = c("2","10","exp(1)")){
  base = match.arg(log_base)
  if(min(p.value) <=0 || max(p.value) >1){
    stop("All elements of p.value must lie in (0,1]!")
  }

    vec = -1*log(x=p.value, base=s_base(base))
    if(base == "exp(1)"){
      base = "e"
    }
    names(vec) = paste0("S","_",base)
    return(vec)
  }


s_base = function(base = c("2","10","exp(1)","none")){
  eval(str2lang(match.arg(base)))
}


#' @family htest
#' @name p_calibrate
#' @export

minbf_pvalue = function(p.value,
                        calc = c("plog","qlog")){
  calc = match.arg(calc)
  if(min(p.value) <=0 || max(p.value) >1){
    stop("All elements of p.value must lie in (0,1]!")
  }

  if(calc == "plog"){
    if(p.value < 1/exp(1)){
      val = (-exp(1) * p.value * log(p.value, base = exp(1)))
    } else {
      val = 1
    }
  }

  if(calc = "qlog"){
    if(p.value < 1-1/exp(1)){
      val = (-exp(1) * (1-p.value) * log(1-p.value, base = exp(1)))
    } else {
      val = 1
    }
  }

  names(val) = "minBF[null]"
  return(val)
}

#' @family htest
#' @name p_calibrate
#' @export

minbf_t = function(p.value, df, alternative = "two.sided"){
  if(min(p.value) <=0 || max(p.value) >1)
    stop("All elements of p.value must lie in (0,1]!")

  if(alternative != "two.sided"){
    k <- length(p)
    m <- length(n)
    log_minBF <- matrix(NA, ncol=k, nrow=m)
    for(j in 1:m){
      t <- qt(p.value, df=df[j], lower.tail = FALSE)
      log_minBF[j, ] <- ifelse(t <= 0, 0, -(df[j]+1)/2*log(1+1/(df[j])*t^2))
    }
    if(k==1 | m==1)
      log_minBF <- as.vector(log_minBF)
  }else {
    k <- length(p)
    m <- length(n)
    log_minBF <- matrix(NA, ncol=k, nrow=m)
    # df <- n-2
    for(j in 1:m){
      for(i in 1:k){
        t.star <- qt(p.value[i]/2, df=df[j], lower.tail = FALSE)
        res <- optimize(bf.ft, p=p[i], df=df[j], log=TRUE, lower=0, upper=t.star+1, maximum=TRUE)
        log_minBF[j,i] <- - res$objective
      }
    }
    if(k==1 | m==1)
      log_minBF <- as.vector(log_minBF)
  }

  result <- exp(log_minBF)
  names(result) = "minBF[null]"
  return(result)
}


#' @family htest
#' @name p_calibrate
#' @export

# approx. Bayes factor to approx. posterior probability (of the hypothesis in the *numerator* of the BF)
pp_from_bf <- function(bf, prior_pr = 0.5){
  prior_odds <- prior_pr/(1-prior_pr)
  log_post_odds <- log(prior_odds) + log(bf)
  post_prob <- exp(log_post_odds)/(1+ exp(log_post_odds))
  names(post_prob) = "approximate posterior probability"
  return (post_prob)
}

