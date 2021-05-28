#' @title The Cohen's d independent distribution
#' @description Density, distribution function, quantile function and random generation


#' @describeIn dist_smd_ind Density function
#' @export

ddist_smd_ind <- function(x, at, ...){
  ncp <- x[[1]]$ncp
  sigma <- x[[1]]$sigma
  df <- x[[1]]$df

  stats::dt((at - df)/sigma, df, ncp) / sigma
}

#' @describeIn dist_smd_ind Quantile function
#' @export
qdist_smd_ind <- function(x, p, ...){
  ncp <- x[[1]]$ncp
  df = x[[1]]$df
  n1 = x[[1]]$n1
  n2 = x[[1]]$n2

  stats::qt(p, df, ncp) * sqrt(1/n1+1/n2)
}

#' @describeIn dist_smd_ind Probability function
#' @export
pdist_smd_ind <- function(x, q, ...){
  ncp <- x[[1]]$ncp
  sigma <- x[[1]]$sigma
  df <- x[[1]]$df

  # sqrt((n1+n2)/(n1*n2))*rt(1000000,df=df,ncp=ncp)
  # sqrt((x[["n1"]]+x[["n2"]])/(x[["n1"]]*x[["n2"]]))

  stats::pt((q)/sigma, df, ncp)
}

#' @describeIn dist_smd_ind Random generator
#' @export
rdist_smd_ind <- function(x, times, ...){
  ncp <- x[[1]]$ncp
  sigma <- x[[1]]$sigma
  df <- x[[1]]$df
  n1 = x[[1]]$n1
  n2 = x[[1]]$n2
  # sqrt((n1+n2)/(n1*n2))*rt(1000000,df=df,ncp=ncp)
  # sqrt((x[["n1"]]+x[["n2"]])/(x[["n1"]]*x[["n2"]]))

  stats::rt(n = times, df = df, ncp = ncp) * sqrt((n1+n2)/(n1*n2))
}

