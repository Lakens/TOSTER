#' Power analysis for TOST for dependent t-test (raw scores)
#' @param alpha alpha used for the test (e.g., 0.05)
#' @param statistical_power desired power (e.g., 0.8)
#' @param low_eqbound lower equivalence bounds (e.g., -0.5) expressed in standardized mean difference (Cohen's dz)
#' @param high_eqbound upper equivalence bounds (e.g., 0.5) expressed in standardized mean difference (Cohen's dz)
#' @param sdif standard deviation of the difference scores
#' @param epsilon By default, true effect is assumed to be 0. If you want to perform an equivalence test when expecting a non-zero effect, specify the expected standardized effect size in Cohen's dz as epsilon (e.g., 0.05)
#' @return Returns a string summarizing the power analysis, and a numeric variable for the number of pairs needed
#' @examples 
#' powerTOSTpaired.raw(alpha=0.05,statistical_power=0.8,low_eqbound=-3, high_eqbound=3, sdif=10, epsilon=0)
#' @export

powerTOSTpaired.raw<-function(alpha, statistical_power, low_eqbound, high_eqbound, sdif, epsilon){
  if(missing(epsilon)) {
    epsilon<-0
  }
  NT1<-ceiling((qnorm(1-alpha)+qnorm(1-((1-statistical_power)/2)))^2/(low_eqbound/sdif-epsilon/sdif)^2)
  NT2<-ceiling((qnorm(1-alpha)+qnorm(1-((1-statistical_power)/2)))^2/(high_eqbound/sdif-epsilon/sdif)^2)
  N<-max(NT1,NT2)
  message(cat("The required sample size to achieve",100*statistical_power,"% power with equivalence bounds of",low_eqbound,"and",high_eqbound,"is",N,"pairs"))
  return(N)
}