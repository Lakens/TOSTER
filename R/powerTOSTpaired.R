#' Power analysis for TOST for dependent t-test (Cohen's dz)
#' @param alpha alpha used for the test (e.g., 0.05)
#' @param statistical_power desired power (e.g., 0.8)
#' @param low_eqbound_dz lower equivalence bounds (e.g., -0.5) expressed in standardized mean difference (Cohen's dz)
#' @param high_eqbound_dz upper equivalence bounds (e.g., 0.5) expressed in standardized mean difference (Cohen's dz)
#' @return Returns a string summarizing the power analysis, and a numeric variable for the number of pairs needed
#' @examples 
#' ## Sample size for alpha = 0.05, 80% power, equivalence bounds of 
#' ## Cohen's dz = -0.3 and Cohen's d = 0.3, and assuming true effect = 0
#' powerTOSTpaired(alpha=0.05,statistical_power=0.8,low_eqbound_dz=-0.3,high_eqbound_dz=0.3)
#' @section References:
#' Chow, S.-C., Wang, H., & Shao, J. (2007). Sample Size Calculations in Clinical Research, Second Edition - CRC Press Book. Formula 3.1.9
#' @importFrom stats pnorm pt qnorm qt
#' @importFrom graphics abline plot points segments title
#' @export

powerTOSTpaired<-function(alpha, statistical_power, low_eqbound_dz, high_eqbound_dz){
  NT1<-(qnorm(1-alpha)+qnorm(1-((1-statistical_power)/2)))^2/(low_eqbound_dz)^2 
  NT2<-(qnorm(1-alpha)+qnorm(1-((1-statistical_power)/2)))^2/(high_eqbound_dz)^2
  N<-ceiling(max(NT1,NT2))
  message(cat("The required sample size to achieve",100*statistical_power,"% power with equivalence bounds of",low_eqbound_dz,"and",high_eqbound_dz,"is",N,"pairs"))
  return(N)
}