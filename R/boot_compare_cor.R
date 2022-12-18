boot_compare_cor <- function(x1, y1, x2, y2,
                     alternative = c("two.sided", "less", "greater"),
                     method = c("pearson", "kendall", "spearman","winsorized","bendpercent"),
                     alpha = 0.05,
                     null = 0,
                     TOST = FALSE,
                     R = 1999){
  DNAME <- paste(deparse(substitute(x1)), "and", deparse(substitute(y1)),
                 "vs.",
                 deparse(substitute(x2)), "and", deparse(substitute(y2)))
  nboot = R
  null.value = null
  if(!is.vector(x1) || !is.vector(x2) || !is.vector(y1) || !is.vector(y2)){
    stop("x1, x2, y1, y2  must be vectors.")
  }
  if(length(x1)!=length(y1)){
    stop("x1 and y1 do not have equal lengths.")
  }
  if(length(x2)!=length(y2)){
    stop("x2 and y2 do not have equal lengths.")
  }

  if(TOST && null <=0){
    stop("positive value for null must be supplied if using TOST.")
  }
  if(TOST){
    alternative = "less"
  }

  if(alternative != "two.sided"){
    ci = 1 - alpha*2
    intmult = c(1,1)
  } else {
    ci = 1 - alpha
    if(TOST){
      intmult = c(1,1)
    } else if(alternative == "less"){
      intmult = c(1,NA)
    } else {
      intmult = c(NA,1)
    }
  }

  alternative = match.arg(alternative)
  method = match.arg(method)

  df <- cbind(x1,y1)
  df <- df[complete.cases(df), ]
  n1 <- nrow(df)
  x1 <- df[,1]
  y1 <- df[,2]
  df <- cbind(x2,y2)
  df <- df[complete.cases(df), ]
  n2 <- nrow(df)
  x2 <- df[,1]
  y2 <- df[,2]

  if(method %in% c("bendpercent","winsorized")){
    if(method == "bendpercent"){
      r1 <- pbcor(x1, y1, ...)
      r2 <- pbcor(x2, y2, ...)
      data <- matrix(sample(n, size=n*nboot, replace=TRUE), nrow=nboot)
      bvec <- apply(data, 1, .corboot_pbcor, x, y, ...) # get bootstrap results corr
      # bootstrap
      data1 <- matrix(sample(n1, size=n1*nboot, replace=TRUE), nrow=nboot)
      bvec1 <- apply(data1, 1, .corboot_pbcor, x1, y1,   ...) # A 1 by nboot matrix.
      data2 <- matrix(sample(n2, size=n2*nboot, replace=TRUE), nrow=nboot)
      bvec2 <- apply(data2, 1, .corboot_pbcor, x2, y2,  ...) # A 1 by nboot matrix.
    }

    if(method == "winsorized"){
      r1 <- wincor(x1, y1, ...)
      r2 <- wincor(x2, y2, ...)
      data <- matrix(sample(n, size=n*nboot, replace=TRUE), nrow=nboot)
      bvec <- apply(data, 1, .corboot_wincor, x, y, ...) # get bootstrap results corr
      # bootstrap
      data1 <- matrix(sample(n1, size=n1*nboot, replace=TRUE), nrow=nboot)
      bvec1 <- apply(data1, 1, .corboot_wincor, x1, y1,   ...) # A 1 by nboot matrix.
      data2 <- matrix(sample(n2, size=n2*nboot, replace=TRUE), nrow=nboot)
      bvec2 <- apply(data2, 1, .corboot_wincor, x2, y2,  ...) # A 1 by nboot matrix
    }


  } else {
    # correlations
    r1 <- cor(x1,y1,method = method)
    r2 <- cor(x2,y2,method = method)
    # bootstrap
    data1 <- matrix(sample(n1, size=n1*nboot, replace=TRUE), nrow=nboot)
    bvec1 <- apply(data1, 1, .corboot, x1, y1, cor, method = method, ...) # A 1 by nboot matrix.
    data2 <- matrix(sample(n2, size=n2*nboot, replace=TRUE), nrow=nboot)
    bvec2 <- apply(data2, 1, .corboot, x2, y2, cor, method = method, ...) # A 1 by nboot matrix.
  }


  bvec <- bvec1 - bvec2
  bsort <- sort(bvec)
  boot.cint = quantile(bvec, c((1 - ci) / 2, 1 - (1 - ci) / 2))
  attr(boot.cint, "conf.level") <- ci
  # p value
  if(alternative == "two.sided"){
    phat <- (sum(bvec < null.value)+.5*sum(bvec==0))/nboot
    sig <- 2 * min(phat, 1 - phat)
  }
  if(alternative == "greater"){
    sig <- 1 - sum(bvec >= null.value)/nboot
  }
  if(alternative == "less"){
    sig <- 1 - sum(bvec <= null.value)/nboot
  }
  if(TOST){
    sig2 <- 1 - sum(bvec >= -1*null.value)/nboot
    sig = max(sig,sig2)
  }


  if (method == "pearson") {
    # Pearson # Fisher
    method2 <- "Bootstrapped difference in Pearson's correlation"
    names(null.value) = "difference in correlation"
    rfinal = c(cor = r1-r2)
  }
  if (method == "spearman") {
    method2 <- "Bootstrapped difference in Spearman's rho"
    #  # Fieller adjusted
    rfinal = c(rho = r1-r2)
    names(null.value) = "difference in rho"
  }
  if (method == "kendall") {
    method2 <- "Bootstrapped difference in Kendall's tau"
    # # Fieller adjusted
    rfinal = c(tau = r1-r2)
    names(null.value) = "difference in tau"
  }
  if (method == "bendpercent") {
    method2 <- "Bootstrapped difference in percentage bend correlation pb"
    # # Fieller adjusted
    rfinal = c(pb = est)
    names(null.value) = "difference in pb"
  }
  if (method == "winsorized") {
    method2 <- "Bootstrapped difference in Winsorized correlation wincor"
    # # Fieller adjusted
    rfinal = c(win = est)
    names(null.value) = "differnce in win"
  }
  # Store as htest
  rval <- list(p.value = sig,
               conf.int = boot.cint,
               estimate = rfinal,
               stderr = sd(bvec,na.rm=TRUE),
               null.value = null.value,
               alternative = alternative,
               method = method2,
               data.name = DNAME,
               boot = list(diff = bvec,
                           r1 = bvec1,
                           r2 = bvec2),
               call = match.call())
  class(rval) <- "htest"
  return(rval)
}
