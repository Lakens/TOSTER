p_from_z = function(x,
                    alternative = "two.sided"){
  if(alternative == "two.sided"){
    2*pnorm(-abs(unlist(x)))
  } else  if (alternative  == "greater"){
    pnorm(x, lower.tail = FALSE)
  } else if (alternative  == "less"){
    pnorm(x, lower.tail = TRUE)
  } else{
    stop("alternative must be two.sided, greater, or less")
  }

}



se_dz = function(smd,n){
  sqrt( 1/n + (smd^2/(2*n)) )
}

se_ds = function(smd,n){
  if(length(n) == 1){
    n = c(n,n)
  }

  sqrt((n[1]+n[2])/(n[1]*n[2]) + smd^2/(2*(n[1]+n[2])))

}

ci_perc = fucntion(vec,
                   alternative = "two.sided",
                   alpha = 0.05){
  if(alternative == "two.sided"){
    alpha = alpha/2
    res = quantile(vec, c(alpha,1-alpha))
  } else {
    res = quantile(vec, c(alpha,1-alpha))
  }

}
