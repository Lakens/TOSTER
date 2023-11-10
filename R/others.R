

tost_decision = function(hypothesis = "EQU",
                         mu_text,
                         pvalue,
                         pTOST,
                         alpha){
  if (hypothesis == "EQU"){
    if(pvalue <= alpha && pTOST <= alpha){
      combined_outcome <- paste0("NHST: reject null significance hypothesis that the effect is equal to ",mu_text," \n",
                                 "TOST: reject null equivalence hypothesis")
    }
    if(pvalue < alpha && pTOST > alpha){
      combined_outcome <- paste0("NHST: reject null significance hypothesis that the effect is equal to ",mu_text," \n",
                                 "TOST: don't reject null equivalence hypothesis")
      # paste0("statistically different from ",mu_text," and not statistically equivalent")
    }
    if(pvalue > alpha && pTOST <= alpha){
      combined_outcome <- paste0("NHST: don't reject null significance hypothesis that the effect is equal to ",mu_text," \n",
                                 "TOST: reject null equivalence hypothesis")
      #paste0("statistically not different from ",mu_text," and statistically equivalent")
    }
    if(pvalue > alpha && pTOST > alpha){
      combined_outcome <- paste0("NHST: don't reject null significance hypothesis that the effect is equal to ",mu_text," \n",
                                 "TOST: don't reject null equivalence hypothesis")
      #paste0("statistically not different from ",mu_text," and not statistically equivalent")
    }
  } else {
    if(pvalue <= alpha && pTOST <= alpha){
      combined_outcome <- paste0("NHST: reject null significance hypothesis that the effect is equal to ",mu_text," \n",
                                 "TOST: reject null MET hypothesis")
      #paste0("statistically different from ",mu_text," and statistically greater than the minimal effect threshold")
    }
    if(pvalue < alpha && pTOST > alpha){
      combined_outcome <- paste0("NHST: reject null significance hypothesis that the effect is equal to ",mu_text," \n",
                                 "TOST: don't reject null MET hypothesis")
      #paste0("statistically different from ",mu_text," but not statistically greater than the minimal effect threshold")
    }
    if(pvalue > alpha && pTOST <= alpha){
      combined_outcome <- paste0("NHST: don't reject null significance hypothesis that the effect is equal to ",mu_text," \n",
                                 "TOST: reject null MET hypothesis")
      #paste0("statistically not different from ",mu_text," and statistically greater than the minimal effect threshold")
    }
    if(pvalue > alpha && pTOST > alpha){
      combined_outcome <- paste0("NHST: don't reject null significance hypothesis that the effect is equal to ",mu_text," \n",
                                 "TOST: don't reject null MET hypothesis")
      #paste0("statistically not different from ",mu_text," and not statistically greater than the minimal effect threshold")
    }
  }
  return(combined_outcome)
}

# Bootstrap CI functions ------

bca <- function(vector, alpha = .05){
  conf.level = 1-alpha
  if(var(vector)==0){
    lower <- mean(vector)
    upper <- mean(vector)
    return(c(lower, upper))
  }

  if(max(vector)==Inf | min(vector)==-Inf){
    stop("bca bootstrap CIs do not work when some values are infinite")
  }

  low <- (1 - conf.level)/2
  high <- 1 - low
  sims <- length(vector)
  z.inv <- length(vector[vector < mean(vector)])/sims
  z <- qnorm(z.inv)
  U <- (sims - 1) * (mean(vector, na.rm=TRUE) - vector)
  top <- sum(U^3)
  under <- 6 * (sum(U^2))^{3/2}
  a <- top / under
  lower.inv <-  pnorm(z + (z + qnorm(low))/(1 - a * (z + qnorm(low))))
  lower <- quantile(vector, lower.inv, names=FALSE)
  upper.inv <-  pnorm(z + (z + qnorm(high))/(1 - a * (z + qnorm(high))))
  upper <- quantile(vector, upper.inv, names=FALSE)
  return(c(lower, upper))
}

perc <- function(vector, alpha = 0.05){
  conf.level = 1-alpha
  if(var(vector)==0){
    lower <- mean(vector)
    upper <- mean(vector)
    return(c(lower, upper))
  }

  low <- (1 - conf.level)/2
  high <- 1 - low

  lower <- quantile(vector, low, names=FALSE)
  upper <- quantile(vector, high, names=FALSE)
  return(c(lower, upper))
}

