#' Power analysis for TOST for difference between two proportions using Z-test (pooled)
#' @param alpha alpha used for the test (e.g., 0.05)
#' @param statistical_power desired power (e.g., 0.8)
#' @param prop1 expected proportion in control condition
#' @param prop2 expected proportion in the experimental condition
#' @param N sample size (e.g., 108)
#' @param low_eqbound_prop lower equivalence bounds (e.g., -0.05) expressed in proportion
#' @param high_eqbound_prop upper equivalence bounds (e.g., 0.05) expressed in proportion
#' @return Calculate either achieved power, equivalence bounds, or required N.
#' Returns a string summarizing the power analysis, and a numeric variable for number of observations, equivalence bounds, or power.
#' @examples
#' ## Sample size for alpha = 0.05, 90% power, assuming true effect prop1 = prop 2 = 0.5,
#' equivalence bounds of 0.4 and 0.6 (so low_eqbound_prop = -0.1 and high_eqbound_prop = 0.1)
#'
#' powerTOSTproptwo(alpha=0.05, statistical_power=0.9, prop1 = 0.5, prop2 = 0.5, low_eqbound_prop=-0.1, high_eqbound_prop=0.1)
#'
#' ## Power for alpha = 0.05, N 542 , assuming true effect prop1 = prop 2 = 0.5,
#' equivalence bounds of 0.4 and 0.6 (so low_eqbound_prop = -0.1 and high_eqbound_prop = 0.1)
#'
#' powerTOSTproptwo(alpha=0.05, N=542, prop1 = 0.5, prop2 = 0.5, low_eqbound_prop=-0.1, high_eqbound_prop=0.1)
#'
#' ## Equivalence bounds for alpha = 0.05, N 542 , assuming true effect prop1 = prop 2 = 0.5,
#' and 90% power
#'
#' powerTOSTproptwo(alpha=0.05, statistical_power=0.9, N=542, prop1 = 0.5, prop2 = 0.5)
#'
#' #Example 4.2.4 from Chow, Wang, & Shao (2007, p. 93)
#' powerTOSTproptwo(alpha=0.05, statistical_power=0.8, prop1 = 0.75, prop2 = 0.8, low_eqbound_prop=-0.2, high_eqbound_prop=0.2)
#'
#' # Example 5 from Julious & Campbell (2012, p. 2932)
#' powerTOSTproptwo(alpha=0.025, statistical_power=0.9, prop1 = 0.8, prop2 = 0.8, low_eqbound_prop=-0.1, high_eqbound_prop=0.1)
#' # From Machin, D. (Ed.). (2008). Sample size tables for clinical studies (3rd ed). Chichester, West Sussex, UK ; Hoboken, NJ: Wiley-Blackwell.
#' # Example 9.4b equivalence of two proportions (p. 113) #
#' powerTOSTproptwo(alpha=0.010, statistical_power=0.8, prop1 = 0.5, prop2 = 0.5, low_eqbound_prop=-0.2, high_eqbound_prop=0.2)/2
#' @section References:
#' Silva, G. T. da, Logan, B. R., & Klein, J. P. (2008). Methods for Equivalence and Noninferiority Testing. Biology of Blood and Marrow Transplantation: Journal of the American Society for Blood and Marrow Transplantation, 15(1 Suppl), 120-127. https://doi.org/10.1016/j.bbmt.2008.10.004
#' Julious, S. A. & Camprop2ell, M. J. (2012). Tutorial in biostatistics: sample sizes for prop1rallel group clinical trials with binary data. Statistics in Medicine, 31:2904-2936.
#' Chow, S.-C., Wang, H., & Shao, J. (2007). Sample Size Calculations in Clinical Research, Second Edition (2 edition). Boca Raton: Chapman and Hall/CRC.
#' @importFrom stats pnorm pt qnorm qt
#' @export

powerTOSTproptwo<-function(alpha, statistical_power, prop1, prop2, N, low_eqbound_prop, high_eqbound_prop){
  if(missing(N)) {
    NT1<-(prop1*(1-prop1)+prop2*(1-prop2))*((qnorm(1-alpha)+qnorm(1-((1-statistical_power)/2)))/(abs(prop1-prop2)-abs(low_eqbound_prop)))^2
    NT2<-(prop1*(1-prop1)+prop2*(1-prop2))*((qnorm(1-alpha)+qnorm(1-((1-statistical_power)/2)))/(abs(prop1-prop2)-abs(high_eqbound_prop)))^2
    N<-max(NT1,NT2)
    message(cat("The required sample size to achieve",100*statistical_power,"% power with equivalence bounds of",low_eqbound_prop,"and",high_eqbound_prop,"is",ceiling(N)))
    return(N)
  }
  if(missing(statistical_power)) {
    statistical_power1<-2*(pnorm((abs(prop1-prop2)-low_eqbound_prop)/sqrt(prop1*(1-prop1)/N+prop2*(1-prop2)/N)-qnorm(1-alpha))+pnorm(-(abs(prop1-prop2)-low_eqbound_prop)/sqrt(prop1*(1-prop1)/N+prop2*(1-prop2)/N)-qnorm(1-alpha)))-1
    statistical_power2<-2*(pnorm((abs(prop1-prop2)-high_eqbound_prop)/sqrt(prop1*(1-prop1)/N+prop2*(1-prop2)/N)-qnorm(1-alpha))+pnorm(-(abs(prop1-prop2)-high_eqbound_prop)/sqrt(prop1*(1-prop1)/N+prop2*(1-prop2)/N)-qnorm(1-alpha)))-1
    statistical_power<-min(statistical_power1,statistical_power2)
    if(statistical_power<0) {statistical_power<-0}
    message(cat("The statistical power is",round(100*statistical_power,2),"% for equivalence bounds of",low_eqbound_prop,"and",high_eqbound_prop,"."))
    return(statistical_power)
  }
  if(missing(low_eqbound_prop) && missing(high_eqbound_prop)) {
    low_eqbound_prop<--sqrt((prop1*(1-prop1)+prop2*(1-prop2))*((qnorm(1-alpha)+qnorm(1-((1-statistical_power)/2))))^2/N)
    high_eqbound_prop<-sqrt((prop1*(1-prop1)+prop2*(1-prop2))*((qnorm(1-alpha)+qnorm(1-((1-statistical_power)/2))))^2/N)
    message(cat("The equivalence bounds to achieve",100*statistical_power,"% power with N =",N,"are",round(low_eqbound_prop,2),"and",round(high_eqbound_prop,2),"."))
    bounds<-c(low_eqbound_prop,high_eqbound_prop)
    return(bounds)
  }
}
