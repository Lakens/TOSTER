#' Power analysis for TOST for dependent t-test (raw scores).
#' @param alpha alpha used for the test (e.g., 0.05)
#' @param statistical_power desired power (e.g., 0.8)
#' @param N number of pairs (e.g., 96)
#' @param low_eqbound lower equivalence bounds (e.g., -0.5) expressed in raw mean difference
#' @param high_eqbound upper equivalence bounds (e.g., 0.5) expressed in raw mean difference
#' @param sdif standard deviation of the difference scores
#' @return Calculate either achieved power, equivalence bounds, or required N, assuming a true effect size of 0.
#' Returns a string summarizing the power analysis, and a numeric variable for number of observations, equivalence bounds, or power.
#' @examples
#' ## Sample size for alpha = 0.05, 80% power, equivalence bounds of -3 and 3 in raw units
#' ## and assuming a standard deviation of the difference scores of 10, and assuming a true effect = 0
#' powerTOSTpaired.raw(alpha=0.05,statistical_power=0.8,low_eqbound=-3, high_eqbound=3, sdif=10)
#'
#' ## Sample size for alpha = 0.05, N = 96 pairs, equivalence bounds of -3 and 3 in raw units
#' ## and assuming a standard deviation of the difference scores of 10, and assuming a true effect = 0
#' powerTOSTpaired.raw(alpha=0.05,N=96,low_eqbound=-3, high_eqbound=3, sdif=10)
#'
#' ## Equivalence bounds for alpha = 0.05, N = 96 pairs, statistical power of 0.8
#' ## and assuming a standard deviation of the difference scores of 10, and assuming a true effect = 0
#' powerTOSTpaired.raw(alpha=0.05,N=96, statistical_power=0.8, sdif=10)
#' @section References:
#' Chow, S.-C., Wang, H., & Shao, J. (2007). Sample Size Calculations in Clinical Research, Second Edition - CRC Press Book. Formula 3.1.9
#' @importFrom stats pnorm pt qnorm qt
#' @importFrom graphics abline plot points segments title
#' @export

powerTOSTpaired.raw<-function(alpha, statistical_power, low_eqbound, high_eqbound, sdif){
  message("Note: this function is defunct. Please use power_t_TOST instead")
  NT1<-(qnorm(1-alpha)+qnorm(1-((1-statistical_power)/2)))^2/(low_eqbound/sdif)^2
  NT2<-(qnorm(1-alpha)+qnorm(1-((1-statistical_power)/2)))^2/(high_eqbound/sdif)^2
  N<-ceiling(max(NT1,NT2))
  message(cat("The required sample size to achieve",100*statistical_power,"% power with equivalence bounds of",low_eqbound,"and",high_eqbound,"is",N,"pairs"))
  return(N)
}

powerTOSTpaired.raw<-function(alpha, statistical_power, N, sdif, low_eqbound, high_eqbound){
  if(missing(N)) {
    NT1<-(qnorm(1-alpha)+qnorm(1-((1-statistical_power)/2)))^2/(low_eqbound/sdif)^2
    NT2<-(qnorm(1-alpha)+qnorm(1-((1-statistical_power)/2)))^2/(high_eqbound/sdif)^2
    N<-max(NT1,NT2)
    message(cat("The required sample size to achieve",100*statistical_power,"% power with equivalence bounds of",low_eqbound,"and",high_eqbound,"is",ceiling(N),"pairs"))
    return(N)
  }
  if(missing(statistical_power)) {
    statistical_power1<-2*(pnorm((abs(low_eqbound)/sdif*sqrt(N))-qnorm(1-alpha))+pnorm(-(abs(low_eqbound)/sdif*sqrt(N))-qnorm(1-alpha)))-1
    statistical_power2<-2*(pnorm((abs(high_eqbound)/sdif*sqrt(N))-qnorm(1-alpha))+pnorm(-(abs(high_eqbound)/sdif*sqrt(N))-qnorm(1-alpha)))-1
    statistical_power<-min(statistical_power1,statistical_power2)
    if(statistical_power<0) {statistical_power<-0}
    message(cat("The statistical power is",round(100*statistical_power,2),"% for equivalence bounds of",low_eqbound,"and",high_eqbound,"."))
    return(statistical_power)
  }
  if(missing(low_eqbound) && missing(high_eqbound)) {
    low_eqbound<--sqrt((qnorm(1-alpha)+qnorm(1-((1-statistical_power)/2)))^2/N)*sdif
    high_eqbound<-sqrt((qnorm(1-alpha)+qnorm(1-((1-statistical_power)/2)))^2/N)*sdif
    message(cat("The equivalence bounds to achieve",100*statistical_power,"% power with N =",N,"are",round(low_eqbound,2),"and",round(high_eqbound,2),"."))
    bounds<-c(low_eqbound,high_eqbound)
    return(bounds)
  }
}
