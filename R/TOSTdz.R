#' TOST function
#' @param m1 mean of group 1
#' @param m2 mean of group 2
#' @param sd1 standard deviation of group 1
#' @param sd2 standard deviation of group 2
#' @param n sample size (pairs)
#' @param r12 correlation of dependent variable between group 1 and  group 2
#' @param low_eqbound_dz lower equivalence bounds (e.g., -0.5) expressed in standardized mean difference (Cohen's dz)
#' @param high_eqbound_dz upper equivalence bounds (e.g., 0.5) expressed in standardized mean difference (Cohen's dz)
#' @param alpha alpha level (default = 0.05)
#' @return Returns dataframe with t-value, degrees of freedom, p-value, and lower and upper limit of confidence interval
#' @examples
#' TOSTdz(n=65,m1=5.830769,m2=5.7538462,sd1=1.1668498,sd2=1.2994081,r12=0.744992,low_eqbound_dz=-0.337341706,high_eqbound_dz=0.337341706, alpha=0.05)
#' @export
#' 

TOSTdz<-function(n,m1,m2,sd1,sd2,r12,low_eqbound_dz, high_eqbound_dz, alpha){
  sdif<-sqrt(sd1^2+sd2^2-2*r12*sd1*sd2)
  se<-sdif/sqrt(n)
  t<-(m2-m1)/se
  degree_f<-n-1
  pttest<-2*pt(abs(t), degree_f, lower=FALSE)
  t1<-((m2-m1)+(low_eqbound_dz*sdif))/se
  p1<-1-pt(t1, degree_f, lower=FALSE)
  t2<-((m2-m1)+(high_eqbound_dz*sdif))/se
  p2<-pt(t2, degree_f, lower=FALSE)
  ttost<-ifelse(abs(t1) < abs(t2), t1, t2)
  LL_90<-((m2-m1)-qt(1-alpha, degree_f)*se)/sdif
  UL_90<-((m2-m1)+qt(1-alpha, degree_f)*se)/sdif
  ptost<-max(p1,p2)
  results<-data.frame(ttost,degree_f,ptost,LL_90,UL_90)
  LL_95<-((m2-m1)-qt(0.975, degree_f)*se)/sdif
  UL_95<-((m2-m1)+qt(0.975, degree_f)*se)/sdif
  df = data.frame(labels=c(paste(100*(1-2*alpha),"% CI dz", sep=""),"Equivalence Range"), mean=c((m2-m1)/sdif,0), lower=c(LL_90,low_eqbound_dz), upper = c(UL_90,high_eqbound_dz))
  plot(NA, xlim=c(.5,2.5), ylim=c(min(LL_90,low_eqbound_dz)-0.2, max(UL_90,high_eqbound_dz)+0.2), bty="l", xaxt="n", xlab="",ylab="Cohen's dz")
  points(df$mean[1:2], pch=19)
  points(1,UL_95,pch=10)
  points(1,LL_95,pch=10)
  axis(1, 1:2, df$labels)
  segments(1:2,df$lower[1:2],1:2,df$upper[1:2])
  segments(1:1,df$upper[1:1],1:2,df$upper[1:1],lty=3)
  segments(2,0,0,0,lty=2)
  segments(1:1,df$lower[1:1],1:2,df$lower[1:1],lty=3)
  text(2, min(LL_90,low_eqbound_dz)-0.15, paste("P-value",round(ptost, digits=3)), cex = .8)
  text(1, min(LL_90,low_eqbound_dz)-0.15, paste("P-value",round(pttest, digits=3)), cex = .8)
  text(1.5, (m2-m1), paste("dz = ",round((m2-m1)/sdif, digits=3)), cex = .8)
  title(main=paste("Mdif = ",round((m2-m1),digits=3),", 95% CI [",round((LL_95*sdif),digits=3),";",round((UL_95*sdif),digits=3),"]",", p = ",round(pttest,digits=3),sep=""))
  testoutcome<-ifelse(pttest<0.05,"significant","non-significant")
  TOSToutcome<-ifelse(ptost<0.05,"significant","non-significant")
  message(cat("The NHST t-test was ",testoutcome,", t(",degree_f,") = ",t,", p = ",pttest,sep=""))
  message(cat("The equivalence test was ",TOSToutcome,", t(",degree_f,") = ",ttost,", p = ",ptost,sep=""))
  return(results)
}