#' Power analysis for TOST for independent t-test (raw scores)
#' @param alpha alpha used for the test (e.g., 0.05)
#' @param statistical_power desired power (e.g., 0.8)
#' @param low_eqbound lower equivalence bounds (e.g., -0.5) expressed in raw scale units (e.g., scalepoints)
#' @param high_eqbound upper equivalence bounds (e.g., 0.5) expressed in raw scale units (e.g., scalepoints)
#' @param sdpooled specify the pooled standard deviation
#' @return Returns a string summarizing the power analysis, and a numeric variable for the number of observations needed in each group
#' @examples 
#' powerTOSTtwo.raw(alpha=0.05, statistical_power=0.8, low_eqbound=-200, high_eqbound=200, sdpooled=350)
#' @section References:
#' Chow, S.-C., Wang, H., & Shao, J. (2007). Sample Size Calculations in Clinical Research, Second Edition - CRC Press Book. Formula 3.2.4 with k = 1
#' @export

powerTOSTtwo.raw<-function(alpha, statistical_power, low_eqbound, high_eqbound, sdpooled){
  NT1<-2*sdpooled^2*(qnorm(1-alpha)+qnorm(1-((1-statistical_power)/2)))^2/(low_eqbound)^2
  NT2<-2*sdpooled^2*(qnorm(1-alpha)+qnorm(1-((1-statistical_power)/2)))^2/(high_eqbound)^2
  N<-ceiling(max(NT1,NT2))
  message(cat("The required sample size to achieve",100*statistical_power,"% power with equivalence bounds of",low_eqbound,"and",high_eqbound,"is",N,"per group, or", 2*N,"in total."))
  return(N)
}