#' Power analysis for TOST for independent t-test (Cohen's d).
#' @param alpha alpha used for the test (e.g., 0.05)
#' @param statistical_power desired power (e.g., 0.8)
#' @param N sample size per group (e.g., 108)
#' @param low_eqbound_d lower equivalence bounds (e.g., -0.5) expressed in standardized mean difference (Cohen's d)
#' @param high_eqbound_d upper equivalence bounds (e.g., 0.5) expressed in standardized mean difference (Cohen's d)
#' @return Calculate either achieved power (given N, alpha, and equivalence bounds) or required N (given desired power, alpha, and equivalence bounds).
#' Returns a string summarizing the power analysis, and a numeric variable for the number of observations needed in each group or the power.
#' @examples
#' ## Sample size for alpha = 0.05, 80% power, equivalence bounds of
#' ## Cohen's d = -0.4 and Cohen's d = 0.4, assuming true effect = 0
#' powerTOSTtwo(alpha=0.05, statistical_power=0.8, low_eqbound_d=-0.4, high_eqbound_d=0.4)
#'
#' ## Sample size for alpha = 0.05, N = 108 per group, equivalence bounds of
#' ## Cohen's d = -0.4 and Cohen's d = 0.4, assuming true effect = 0
#' powerTOSTtwo(alpha=0.05, N=108, low_eqbound_d=-0.4, high_eqbound_d=0.4)
#' @section References:
#' Chow, S.-C., Wang, H., & Shao, J. (2007). Sample Size Calculations in Clinical Research, Second Edition - CRC Press Book. Formula 3.2.4 with k = 1
#' @importFrom stats pnorm pt qnorm qt
#' @importFrom graphics abline plot points segments title
#' @export

powerTOSTtwo<-function(alpha, statistical_power, N, low_eqbound_d, high_eqbound_d){
  if(missing(N)) {
    NT1<-2*(qnorm(1-alpha)+qnorm(1-((1-statistical_power)/2)))^2/(low_eqbound_d)^2
    NT2<-2*(qnorm(1-alpha)+qnorm(1-((1-statistical_power)/2)))^2/(high_eqbound_d)^2
    N<-ceiling(max(NT1,NT2))
    message(cat("The required sample size to achieve",100*statistical_power,"% power with equivalence bounds of",low_eqbound_d,"and",high_eqbound_d,"is",N,"per group, or", 2*N,"in total."))
    return(N)
  }
  if(missing(statistical_power)) {
    statistical_power1<-2*(pnorm(abs(low_eqbound_d)*sqrt(N/2)-qnorm(1-alpha))+pnorm(-abs(low_eqbound_d)*sqrt(N/2)-qnorm(1-alpha)))-1
    statistical_power2<-2*(pnorm(abs(high_eqbound_d)*sqrt(N/2)-qnorm(1-alpha))+pnorm(-abs(high_eqbound_d)*sqrt(N/2)-qnorm(1-alpha)))-1
    statistical_power<-min(statistical_power1,statistical_power2)
    if(statistical_power<0) {statistical_power<-0}
    message(cat("The statistical power is",round(100*statistical_power,2),"% for equivalence bounds of",low_eqbound_d,"and",high_eqbound_d,"."))
    return(statistical_power)
  }
}
