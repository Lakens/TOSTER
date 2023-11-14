

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

## only an approximation... rather useless
bca <- function(vector, alpha = 0.05){
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


basic <- function(vector, t0, alpha){
  conf = 1-alpha
  qq <- norm.inter(vector, (1 + c(conf, -conf))/2)
  c((2 *  t0 - qq[, 2L]))
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

stud <- function(vector, se, se0, t0, alpha){
  conf = 1-alpha
  z <- (vector - t0)/(se)
  qq <- norm.inter(z, (1 + c(conf, -conf))/2)
  c( ((t0 - (se0) * qq[, 2L])))
}

norm.inter = function (t, alpha)
{
  t <- t[is.finite(t)]
  R <- length(t)
  rk <- (R + 1) * alpha
  if (!all(rk > 1 & rk < R))
    warning("extreme order statistics used as endpoints")
  k <- trunc(rk)
  inds <- seq_along(k)
  out <- inds
  kvs <- k[k > 0 & k < R]
  tstar <- sort(t, partial = sort(union(c(1, R), c(kvs, kvs +
                                                     1))))
  ints <- (k == rk)
  if (any(ints))
    out[inds[ints]] <- tstar[k[inds[ints]]]
  out[k == 0] <- tstar[1L]
  out[k == R] <- tstar[R]
  not <- function(v) xor(rep(TRUE, length(v)), v)
  temp <- inds[not(ints) & k != 0 & k != R]
  temp1 <- qnorm(alpha[temp])
  temp2 <- qnorm(k[temp]/(R + 1))
  temp3 <- qnorm((k[temp] + 1)/(R + 1))
  tk <- tstar[k[temp]]
  tk1 <- tstar[k[temp] + 1L]
  out[temp] <- tk + (temp1 - temp2)/(temp3 - temp2) * (tk1 -
                                                         tk)
  cbind(round(rk, 2), out)
}

