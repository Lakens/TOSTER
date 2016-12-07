#' Power analysis for TOST for dependent t-test (raw scores)
#' @param alpha alpha used for the test (e.g., 0.05)
#' @param statistical_power desired power (e.g., 0.8)
#' @param low_eqbound lower equivalence bounds (e.g., -0.5) expressed in raw mean difference
#' @param high_eqbound upper equivalence bounds (e.g., 0.5) expressed in raw mean difference
#' @param sdif standard deviation of the difference scores
#' @return Returns a string summarizing the power analysis, and a numeric variable for the number of pairs needed
#' @examples 
#' ## Sample size for alpha = 0.05, 80% power, equivalence bounds of -3 and 3 in raw units, and assuming a standard deviation of the difference scores of 10, and assuming a true effect = 0
#' powerTOSTpaired.raw(alpha=0.05,statistical_power=0.8,low_eqbound=-3, high_eqbound=3, sdif=10)
#' @section References:
#' Chow, S.-C., Wang, H., & Shao, J. (2007). Sample Size Calculations in Clinical Research, Second Edition - CRC Press Book. Formula 3.1.9
#' @importFrom stats pnorm pt qnorm qt
#' @importFrom graphics abline plot points segments title
#' @export

powerTOSTpaired.raw<-function(alpha, statistical_power, low_eqbound, high_eqbound, sdif){
  NT1<-(qnorm(1-alpha)+qnorm(1-((1-statistical_power)/2)))^2/(low_eqbound/sdif)^2
  NT2<-(qnorm(1-alpha)+qnorm(1-((1-statistical_power)/2)))^2/(high_eqbound/sdif)^2
  N<-ceiling(max(NT1,NT2))
  message(cat("The required sample size to achieve",100*statistical_power,"% power with equivalence bounds of",low_eqbound,"and",high_eqbound,"is",N,"pairs"))
  return(N)
}