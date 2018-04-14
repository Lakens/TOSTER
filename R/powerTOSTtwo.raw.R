#' Power analysis for TOST for independent t-test (raw scores).
#' @param alpha alpha used for the test (e.g., 0.05)
#' @param statistical_power desired power (e.g., 0.8)
#' @param N sample size per group (e.g., 108)
#' @param low_eqbound lower equivalence bounds (e.g., -0.5) expressed in raw scale units (e.g., scalepoints)
#' @param high_eqbound upper equivalence bounds (e.g., 0.5) expressed in raw scale units (e.g., scalepoints)
#' @param sdpooled specify the pooled standard deviation
#' @return Calculate either achieved power, equivalence bounds, or required N.
#' Returns a string summarizing the power analysis, and a numeric variable for number of observations, equivalence bounds, or power.
#' @examples
#' ## Sample size for alpha = 0.05, 80% power, equivalence bounds of -200 and 200 in raw
#' ## units, assuming pooled standard deviation of 350, and assuming true effect = 0
#' powerTOSTtwo.raw(alpha=0.05,statistical_power=0.8,low_eqbound=-200,high_eqbound=200,sdpooled=350)
#'
#' ## Power for alpha = 0.05, N = 53 per group, equivalence bounds of
#' ## -200 and 200 in raw units, assuming sdpooled = 350 and true effect = 0
#' powerTOSTtwo.raw(alpha=0.05, N=53, low_eqbound=-200, high_eqbound=200, sdpooled=350)
#'
#' ## Equivalence bounds for alpha = 0.05, N = 108 per group, statistical power of
#' ## 0.8, assuming true effect = 0
#' powerTOSTtwo.raw(alpha=0.05, N=53, statistical_power=0.8, sdpooled=350)

#' @section References:
#' Chow, S.-C., Wang, H., & Shao, J. (2007). Sample Size Calculations in Clinical Research, Second Edition - CRC Press Book. Formula 3.2.4 with k = 1
#' @importFrom stats pnorm pt qnorm qt
#' @importFrom graphics abline plot points segments title
#' @export

powerTOSTtwo.raw<-function(alpha, statistical_power, N, sdpooled, low_eqbound, high_eqbound){
  if(missing(N)) {
    NT1<-2*sdpooled^2*(qnorm(1-alpha)+qnorm(1-((1-statistical_power)/2)))^2/(low_eqbound)^2
    NT2<-2*sdpooled^2*(qnorm(1-alpha)+qnorm(1-((1-statistical_power)/2)))^2/(high_eqbound)^2
    N<-max(NT1,NT2)
    message(cat("The required sample size to achieve",100*statistical_power,"% power with equivalence bounds of",low_eqbound,"and",high_eqbound,"is",N,"per group, or", 2*ceiling(N),"in total."))
    return(N)
  }
  if(missing(statistical_power)) {
    statistical_power1<-2*(pnorm(abs(low_eqbound)/sdpooled*sqrt(N/2)-qnorm(1-alpha))+pnorm(-abs(low_eqbound)/sdpooled*sqrt(N/2)-qnorm(1-alpha)))-1
    statistical_power2<-2*(pnorm(abs(high_eqbound)/sdpooled*sqrt(N/2)-qnorm(1-alpha))+pnorm(-abs(high_eqbound)/sdpooled*sqrt(N/2)-qnorm(1-alpha)))-1
    statistical_power<-min(statistical_power1,statistical_power2)
    if(statistical_power<0) {statistical_power<-0}
    message(cat("The statistical power is",round(100*statistical_power,2),"% for equivalence bounds of",low_eqbound,"and",high_eqbound,"."))
    return(statistical_power)
  }
  if(missing(low_eqbound) && missing(high_eqbound)) {
    low_eqbound<--sqrt(2*(qnorm(1-alpha)+qnorm(1-((1-statistical_power)/2)))^2/N)*sdpooled
    high_eqbound<-sqrt(2*(qnorm(1-alpha)+qnorm(1-((1-statistical_power)/2)))^2/N)*sdpooled
    message(cat("The equivalence bounds to achieve",100*statistical_power,"% power with N =",N,"are",round(low_eqbound,2),"and",round(high_eqbound,2),"."))
    bounds<-c(low_eqbound,high_eqbound)
    return(bounds)
  }
}
