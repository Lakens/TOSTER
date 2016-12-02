#' Power analysis for TOST for one-sample t-test (raw scores)
#' @param alpha alpha used for the test (e.g., 0.05)
#' @param statistical_power desired power (e.g., 0.8)
#' @param low_eqbound lower equivalence bounds (e.g., -0.5) expressed in raw scores
#' @param high_eqbound upper equivalence bounds (e.g., 0.5) expressed in raw scores
#' @return Returns a string summarizing the power analysis, and a numeric variable for the number of observations needed.
#' @examples 
#' powerTOSTone.raw(alpha=0.05, statistical_power=0.9, sd = 1, low_eqbound=-0.3, high_eqbound=0.3)
#' @section References:
#' Chow, S.-C., Wang, H., & Shao, J. (2007). Sample Size Calculations in Clinical Research, Second Edition - CRC Press Book. Formula 3.1.9
#' @export


powerTOSTone.raw<-function(alpha, statistical_power, sd, low_eqbound, high_eqbound){
  NT1<-(qnorm(1-alpha)+qnorm(1-((1-statistical_power)/2)))^2/(low_eqbound/sd)^2
  NT2<-(qnorm(1-alpha)+qnorm(1-((1-statistical_power)/2)))^2/(high_eqbound/sd)^2
  N<-ceiling(max(NT1,NT2))
  message(cat("The required sample size to achieve",100*statistical_power,"% power with equivalence bounds of",low_eqbound,"and",high_eqbound,"is",N))
  return(N)
}