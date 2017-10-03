#' Power analysis for TOST for one-sample t-test (raw scores).
#' @param alpha alpha used for the test (e.g., 0.05)
#' @param sd population standard deviation
#' @param statistical_power desired power (e.g., 0.8)
#' @param N sample size (e.g., 108)
#' @param low_eqbound lower equivalence bounds (e.g., -0.5) expressed in raw scores
#' @param high_eqbound upper equivalence bounds (e.g., 0.5) expressed in raw scores
#' @return Calculate either achieved power, equivalence bounds, or required N.
#' Returns a string summarizing the power analysis, and a numeric variable for number of observations, equivalence bounds, or power.
#' @examples
#' ## Sample size for alpha = 0.05, 90% power, equivalence bounds of -0.3 and 0.3 in
#' ## raw units, assuming pooled standard deviation of 1, and assuming true effect = 0
#' powerTOSTone.raw(alpha=0.05, statistical_power=0.9, sd = 1, low_eqbound=-0.3, high_eqbound=0.3)
#'
#' ## Power for sample size of 121, alpha = 0.05, equivalence bounds of
#' ## -0.3 and 0.3 in raw units, assuming pooled standard deviation of 1, and assuming true effect = 0
#'
#' powerTOSTone.raw(alpha=0.05, N=121, sd = 1, low_eqbound=-0.3, high_eqbound=0.3)
#'
#' ## Power for sample size of 121, alpha = 0.05, statistical power of
#' ## 0.9, and assuming true effect = 0
#'
#' powerTOSTone.raw(alpha=0.05, N=121, statistical_power=.9)
#' @section References:
#' Chow, S.-C., Wang, H., & Shao, J. (2007). Sample Size Calculations in Clinical Research, Second Edition - CRC Press Book. Formula 3.1.9
#' @importFrom stats pnorm pt qnorm qt
#' @importFrom graphics abline plot points segments title
#' @export


powerTOSTone.raw<-function(alpha, statistical_power, sd, low_eqbound, high_eqbound){
  NT1<-(qnorm(1-alpha)+qnorm(1-((1-statistical_power)/2)))^2/(low_eqbound/sd)^2
  NT2<-(qnorm(1-alpha)+qnorm(1-((1-statistical_power)/2)))^2/(high_eqbound/sd)^2
  N<-ceiling(max(NT1,NT2))
  message(cat("The required sample size to achieve",100*statistical_power,"% power with equivalence bounds of",low_eqbound,"and",high_eqbound,"is",N))
  return(N)
}

powerTOSTone.raw<-function(alpha, statistical_power, N, sd, low_eqbound, high_eqbound){
  if(missing(N)) {
    NT1<-(qnorm(1-alpha)+qnorm(1-((1-statistical_power)/2)))^2/(low_eqbound/sd)^2
    NT2<-(qnorm(1-alpha)+qnorm(1-((1-statistical_power)/2)))^2/(high_eqbound/sd)^2
    N<-ceiling(max(NT1,NT2))
    message(cat("The required sample size to achieve",100*statistical_power,"% power with equivalence bounds of",low_eqbound,"and",high_eqbound,"is",N))
    return(N)
  }
  if(missing(statistical_power)) {
    statistical_power1<-2*(pnorm((abs(low_eqbound)/sd*sqrt(N))-qnorm(1-alpha))+pnorm(-(abs(low_eqbound)/sd*sqrt(N))-qnorm(1-alpha)))-1
    statistical_power2<-2*(pnorm((abs(high_eqbound)/sd*sqrt(N))-qnorm(1-alpha))+pnorm(-(abs(high_eqbound)/sd*sqrt(N))-qnorm(1-alpha)))-1
    statistical_power<-min(statistical_power1,statistical_power2)
    if(statistical_power<0) {statistical_power<-0}
    message(cat("The statistical power is",round(100*statistical_power,2),"% for equivalence bounds of",low_eqbound,"and",high_eqbound,"."))
    return(statistical_power)
  }
  if(missing(low_eqbound) && missing(high_eqbound)) {
    low_eqbound<--sqrt((qnorm(1-alpha)+qnorm(1-((1-statistical_power)/2)))^2/N)
    high_eqbound<-sqrt((qnorm(1-alpha)+qnorm(1-((1-statistical_power)/2)))^2/N)
    message(cat("The equivalence bounds to achieve",100*statistical_power,"% power with N =",N,"are",round(low_eqbound,2),"and",round(high_eqbound,2),"."))
    bounds<-c(low_eqbound,high_eqbound)
    return(bounds)
  }
}
