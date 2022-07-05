p_from_z = function(x,
                    alternative = "two.sided"){
  if(alternative == "two.sided"){
    2*pnorm(-abs(unlist(x)))
  } else (alternative  == "greater"){
    pnorm(x, lower.tail = FALSE)
  } else (alternative  == "less"){
    pnorm(x, lower.tail = TRUE)
  }

}
## compare paired samples studies (cohen's dz)
compare_dz = function(
    # cohen dz original study
  d_ori,
  # N original study
  n_ori,
  # cohen dz replication study
  d_rep,
  # N replication study
  n_rep) {

  se_ori = sqrt( 1/n_ori + (d_ori^2/(2*n_ori)) )
  se_rep = sqrt( 1/n_rep + (d_rep^2/(2*n_rep)) )
  se_diff = sqrt(se_ori^2 + se_rep^2)
  z = (d_ori - d_rep)/se_diff
  names(z) = "z"
  pval = p_from_z(z)
  alt_meth = "two-sided"
  null = 0
  par = c(n_ori,n_rep)
  names(par) = c("N original", "N replication")
  rval <- list(statistic = z, parameter = par, p.value = pval,
               #conf.int = cint,
               estimate = d_ori-d_rep,
               null.value = 0,
               alternative = "two.sided",
               method = "Difference in Cohen's dz",
               data.name = "Summary Statistics")
  class(rval) <- "htest"
  return(rval)
}
