#' Power analysis for TOST for correlations.
#' @param alpha alpha used for the test (e.g., 0.05)
#' @param statistical_power desired power (e.g., 0.8)
#' @param N number of pairs (e.g., 96)
#' @param low_eqbound_r lower equivalence bounds (e.g., -0.3) expressed in a correlation effect size
#' @param high_eqbound_r upper equivalence bounds (e.g., 0.3) expressed in a correlation effect size
#' @return Calculate either achieved power, equivalence bounds, or required N.
#' Returns a string summarizing the power analysis, and a numeric variable for number of observations, equivalence bounds, or power.
#' @examples
#' ## Sample size for alpha = 0.05, 90% power, equivalence bounds of
#' ## r = -0.1 and r = 0.1, assuming true effect = 0
#' powerTOSTr(alpha=0.05, statistical_power=0.9, low_eqbound_r=-0.1, high_eqbound_r=0.1)
#'
#' ## Sample size for alpha = 0.05, N=536, equivalence bounds of
#' ## r = -0.1 and r = 0.1, assuming true effect = 0
#' powerTOSTr(alpha=0.05, N=536, low_eqbound_r=-0.1, high_eqbound_r=0.1)
#'
#' ## Equivalence bounds for alpha = 0.05, N=536, statistical power of
#' ## 0.9, assuming true effect = 0
#' powerTOSTr(alpha=0.05, N=536, statistical_power=0.9)
#' @importFrom stats pnorm pt qnorm qt
#' @importFrom graphics abline plot points segments title
#' @export

powerTOSTr<-function(alpha, statistical_power, N, low_eqbound_r, high_eqbound_r){
  if(missing(N)) {
    NT1<-2*(qnorm(1-alpha)+qnorm(1-((1-statistical_power)/2)))^2/(((2*low_eqbound_r)/sqrt(1-low_eqbound_r^2)))^2
    NT2<-2*(qnorm(1-alpha)+qnorm(1-((1-statistical_power)/2)))^2/(((2*high_eqbound_r)/sqrt(1-high_eqbound_r^2)))^2
    N<-max(NT1,NT2)
    message(cat("The required sample size to achieve",100*statistical_power,"% power with equivalence bounds of",low_eqbound_r,"and",high_eqbound_r,"is",ceiling(N),"pairs"))
    return(N)
  }
  if(missing(statistical_power)) {
    statistical_power1<-2*(pnorm((abs((((2*low_eqbound_r)/sqrt(1-low_eqbound_r^2))))*sqrt(N/2))-qnorm(1-alpha))+pnorm(-(abs((((2*low_eqbound_r)/sqrt(1-low_eqbound_r^2))))*sqrt(N/2))-qnorm(1-alpha)))-1
    statistical_power2<-2*(pnorm((abs(((2*high_eqbound_r)/sqrt(1-high_eqbound_r^2)))*sqrt(N/2))-qnorm(1-alpha))+pnorm(-(abs(((2*high_eqbound_r)/sqrt(1-high_eqbound_r^2)))*sqrt(N/2))-qnorm(1-alpha)))-1
    statistical_power<-min(statistical_power1,statistical_power2)
    if(statistical_power<0) {statistical_power<-0}
    message(cat("The statistical power is",round(100*statistical_power,2),"% for equivalence bounds of",low_eqbound_r,"and",high_eqbound_r,"."))
    return(statistical_power)
  }
  if(missing(low_eqbound_r) && missing(high_eqbound_r)) {
    low_eqbound_r<-(-sqrt(2*(qnorm(1-alpha)+qnorm(1-((1-statistical_power)/2)))^2/N))/2
    high_eqbound_r<-(sqrt(2*(qnorm(1-alpha)+qnorm(1-((1-statistical_power)/2)))^2/N))/2
    message(cat("The equivalence bounds to achieve",100*statistical_power,"% power with N =",N,"are",round(low_eqbound_r,2),"and",round(high_eqbound_r,2),"."))
    bounds<-c(low_eqbound_r,high_eqbound_r)
    return(bounds)
  }
}
