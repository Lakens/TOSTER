#' Power analysis for TOST for dependent t-test (Cohen's dz).
#' @param alpha alpha used for the test (e.g., 0.05)
#' @param statistical_power desired power (e.g., 0.8)
#' @param N number of pairs (e.g., 96)
#' @param low_eqbound_dz lower equivalence bounds (e.g., -0.5) expressed in standardized mean difference (Cohen's dz)
#' @param high_eqbound_dz upper equivalence bounds (e.g., 0.5) expressed in standardized mean difference (Cohen's dz)
#' @return Calculate either achieved power, equivalence bounds, or required N, assuming a true effect size of 0.
#' Returns a string summarizing the power analysis, and a numeric variable for number of observations, equivalence bounds, or power.
#' @examples
#' ## Sample size for alpha = 0.05, 80% power, equivalence bounds of
#' ## Cohen's dz = -0.3 and Cohen's d = 0.3, and assuming true effect = 0
#' powerTOSTpaired(alpha=0.05,statistical_power=0.8,low_eqbound_dz=-0.3,high_eqbound_dz=0.3)
#'
#' ## Sample size for alpha = 0.05, N = 96 pairs, equivalence bounds of
#' ## Cohen's dz = -0.3 and Cohen's d = 0.3, and assuming true effect = 0
#' powerTOSTpaired(alpha=0.05,N=96,low_eqbound_dz=-0.3,high_eqbound_dz=0.3)
#'
#' ## Equivalence bounds for alpha = 0.05, N = 96 pairs, statistical power of
#' ## 0.8, and assuming true effect = 0
#' powerTOSTpaired(alpha=0.05,N=96,statistical_power=0.8)
#' @section References:
#' Chow, S.-C., Wang, H., & Shao, J. (2007). Sample Size Calculations in Clinical Research, Second Edition - CRC Press Book. Formula 3.1.9
#' @importFrom stats pnorm pt qnorm qt
#' @importFrom graphics abline plot points segments title
#' @export

powerTOSTpaired<-function(alpha, statistical_power, N, low_eqbound_dz, high_eqbound_dz){
  message("Note: this function is defunct. Please use power_t_TOST instead")
  if(missing(N)) {
    NT1<-(qnorm(1-alpha)+qnorm(1-((1-statistical_power)/2)))^2/(low_eqbound_dz)^2
    NT2<-(qnorm(1-alpha)+qnorm(1-((1-statistical_power)/2)))^2/(high_eqbound_dz)^2
    N<-max(NT1,NT2)
    message(cat("The required sample size to achieve",100*statistical_power,"% power with equivalence bounds of",low_eqbound_dz,"and",high_eqbound_dz,"is",ceiling(N),"pairs"))
    return(N)
  }
  if(missing(statistical_power)) {
    statistical_power1<-2*(pnorm((abs(low_eqbound_dz)*sqrt(N))-qnorm(1-alpha))+pnorm(-(abs(low_eqbound_dz)*sqrt(N))-qnorm(1-alpha)))-1
    statistical_power2<-2*(pnorm((abs(high_eqbound_dz)*sqrt(N))-qnorm(1-alpha))+pnorm(-(abs(high_eqbound_dz)*sqrt(N))-qnorm(1-alpha)))-1
    statistical_power<-min(statistical_power1,statistical_power2)
    if(statistical_power<0) {statistical_power<-0}
    message(cat("The statistical power is",round(100*statistical_power,2),"% for equivalence bounds of",low_eqbound_dz,"and",high_eqbound_dz,"."))
    return(statistical_power)
  }
  if(missing(low_eqbound_dz) && missing(high_eqbound_dz)) {
    low_eqbound_dz<--sqrt((qnorm(1-alpha)+qnorm(1-((1-statistical_power)/2)))^2/N)
    high_eqbound_dz<-sqrt((qnorm(1-alpha)+qnorm(1-((1-statistical_power)/2)))^2/N)
    message(cat("The equivalence bounds to achieve",100*statistical_power,"% power with N =",N,"are",round(low_eqbound_dz,2),"and",round(high_eqbound_dz,2),"."))
    bounds<-c(low_eqbound_dz,high_eqbound_dz)
    return(bounds)
  }
}
