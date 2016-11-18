#' TOST raw function
#' @param m1 mean of group 1
#' @param m2 mean of group 2
#' @param sd1 standard deviation of group 1
#' @param sd2 standard deviation of group 2
#' @param n1 sample size in group 1
#' @param n2 sample size in group 2
#' @param low_eqbound lower equivalence bounds (e.g., -0.5) expressed in raw scale units (e.g., scalepoints) 
#' @param high_eqbound upper equivalence bounds (e.g., 0.5) expressed in raw scale units (e.g., scalepoints)
#' @param alpha alpha level (default = 0.05)
#' @param var.equal logical variable indicating whether equal variances assumption is assumed to be TRUE or FALSE. Defaults to FALSE.
#' @return Returns dataframe with t-value, degrees of freedom, p-value, and lower and upper limit of confidence interval
#' @examples
#' TOSTraw(m1=7.83,m2=7.98,sd1=1.21,sd2=1.29,n1=400,n2=400,low_eqbound=-0.25, high_eqbound=0.25,alpha=0.05)
#' @export


TOSTraw<-function(m1,m2,sd1,sd2,n1,n2,low_eqbound, high_eqbound, alpha, var.equal){
  if(missing(alpha)) {
    alpha<-0.05
  }
  if(missing(var.equal)) {
    var.equal<-FALSE
  }
  sdpooled<-sqrt((((n1 - 1)*(sd1^2)) + (n2 - 1)*(sd2^2))/((n1+n2)-2)) #calculate sd pooled
  if(var.equal==TRUE) {
    t1<-(abs(m1-m2)-high_eqbound)/(sdpooled*sqrt(1/n1 + 1/n2))
    degree_f<-n1+n2-2
    p1<-pt(t1, degree_f, lower=FALSE) 
    t2<-(abs(m1-m2)-low_eqbound)/(sdpooled*sqrt(1/n1 + 1/n2))
    p2<-pt(t2, degree_f, lower=TRUE) 
    LL<-(m1-m2)-qt(1-alpha, degree_f)*(sdpooled*sqrt(1/n1 + 1/n2))
    UL<-(m1-m2)+qt(1-alpha, degree_f)*(sdpooled*sqrt(1/n1 + 1/n2))
    t<-(m1-m2)/(sdpooled*sqrt(1/n1 + 1/n2))
    pttest<-2*pt(-abs(t), df=degree_f)
  } else {
    t1<-(abs(m1-m2)-low_eqbound)/sqrt(sd1^2/n1 + sd2^2/n2) #welch's t-test lower bound
    degree_f<-(sd1^2/n1+sd2^2/n2)^2/(((sd1^2/n1)^2/(n1-1))+((sd2^2/n2)^2/(n2-1))) #degrees of freedom for Welch's t-test
    p1<-pt(t1, degree_f, lower=FALSE) #p-value for Welch's TOST t-test
    t2<-(abs(m1-m2)-high_eqbound)/sqrt(sd1^2/n1 + sd2^2/n2) #welch's t-test upper bound
    p2<-pt(t2, degree_f, lower=TRUE) #p-value for Welch's TOST t-test
    t<-(m1-m2)/sqrt(sd1^2/n1 + sd2^2/n2) #welch's t-test NHST
    pttest<-2*pt(-abs(t), df=degree_f) #p-value for Welch's t-test
    LL<-(m1-m2)-qt(1-alpha, degree_f)*sqrt(sd1^2/n1 + sd2^2/n2) #Lower limit for CI Welch's t-test
    UL<-(m1-m2)+qt(1-alpha, degree_f)*sqrt(sd1^2/n1 + sd2^2/n2) #Upper limit for CI Welch's t-test
  }
  ptost<-max(p1,p2) #Get highest p-value for summary TOST result
  ttost<-ifelse(abs(t1) < abs(t2), t1, t2) #Get lowest t-value for summary TOST result
  results<-data.frame(ttost,degree_f,ptost,LL,UL)
  dif<-(m1-m2)
  df = data.frame(labels=c(paste(100*(1-2*alpha),"% CI mean", sep=""),"Equivalence Range"), mean=c(dif,0), lower=c(LL,low_eqbound), upper = c(UL,high_eqbound))
  plot(NA, xlim=c(.5,2.5), ylim=c(min(LL,low_eqbound)-0.5, max(UL,high_eqbound)+0.5), bty="l", xaxt="n", xlab="",ylab="Mean Difference")
  points(df$mean[1:2], pch=19)
  points(1,(m1-m2)-qnorm(0.975)*(sdpooled*sqrt(1/n1 + 1/n2)),pch=10)
  points(1,(m1-m2)+qnorm(0.975)*(sdpooled*sqrt(1/n1 + 1/n2)),pch=10)
  axis(1, 1:2, df$labels)
  segments(1:2,df$lower[1:2],1:2,df$upper[1:2])
  segments(1:1,df$upper[1:1],1:2,df$upper[1:1],lty=3)
  segments(2,0,0,0,lty=2)
  segments(1:1,df$lower[1:1],1:2,df$lower[1:1],lty=3)
  text(2, min(LL,low_eqbound)-0.2, paste("P-value",round(ptost, digits=3)), cex = .8)
  text(1, min(LL,low_eqbound)-0.2, paste("P-value",round(pttest, digits=3)), cex = .8)
  text(1.5, dif, paste("Mdif = ",round(dif, digits=3)), cex = .8)
  title(main=paste("Mdif = ",round(dif,digits=3),", 95% CI [",round(((m1-m2)-qnorm(0.975)*(sdpooled*sqrt(1/n1 + 1/n2))),digits=3),";",round(((m1-m2)+qnorm(0.975)*(sdpooled*sqrt(1/n1 + 1/n2))),digits=3),"]",", p = ",round(pttest,digits=3), sep=""), cex.main=1)
  testoutcome<-ifelse(pttest<0.05,"significant","non-significant")
  TOSToutcome<-ifelse(ptost<0.05,"significant","non-significant")
  message(cat("The NHST t-test was ",testoutcome,", t(",n1+n2-2,") = ",t,", p = ",pttest,sep=""))
  message(cat("The equivalence test was ",TOSToutcome,", t(",n1+n2-2,") = ",ttost,", p = ",ptost,sep=""))
  return(results)
}