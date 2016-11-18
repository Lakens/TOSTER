#' Function to perform power analysis for TOST
#' @param alpha alpha used for the test (e.g., 0.05)
#' @param statistical_power desired power (e.g., 0.8)
#' @param low_eqbound_dz lower equivalence bounds (e.g., -0.5) expressed in standardized mean difference (Cohen's dz)
#' @param high_eqbound_dz upper equivalence bounds (e.g., 0.5) expressed in standardized mean difference (Cohen's dz)
#' @param epsilon By default, true effect is assumed to be 0. If you want to perform an equivalence test when expecting a non-zero effect, specify the expected standardized effect size in Cohen's dz as epsilon (e.g., 0.05)
#' @return Returns a string summarizing the power analysis, and a numeric variable for the number of pairs needed
#' @examples 
#' powerTOSTdz(alpha=0.05,statistical_power=0.8,low_eqbound_dz=-0.3, high_eqbound_dz=0.3, epsilon=0.1)
#' @export

powerTOSTdz<-function(alpha, statistical_power, low_eqbound_dz, high_eqbound_dz, epsilon){
  if(missing(epsilon)) {
    epsilon<-0
  }
  NT1<-ceiling((qnorm(1-alpha)+qnorm(1-((1-statistical_power)/2)))^2/(0-low_eqbound_dz-epsilon)^2)
  NT2<-ceiling((qnorm(1-alpha)+qnorm(1-((1-statistical_power)/2)))^2/(high_eqbound_dz-epsilon)^2)
  N<-max(NT1,NT2)
  message(cat("The required sample size to achieve",100*statistical_power,"% power with equivalence bounds of",low_eqbound_dz,"and",high_eqbound_dz,"is",N,"pairs"))
  return(N)
}