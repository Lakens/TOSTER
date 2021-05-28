#' @title Distribution of smd's d for independent samples
#' #' Distribution of SMD with 2 independent groups.
#'
#' @inheritParams stats::dt
#' @param mult The location parameter of the distribution.
#'   If `ncp == 0` (or `NULL`), this is the median.
#' @param sigma The scale parameter of the distribution.
#'
#' @details
#'
#'   We recommend reading this documentation on
#'   <https://pkg.mitchelloharawild.com/distributional/>, where the math
#'   will render nicely.
#'
#'   In the following, let \eqn{X} be a **central** Students T random variable
#'   with `df` = \eqn{\nu}.
#'
#'   **Support**: \eqn{R}, the set of all real numbers
#'
#'   **Mean**: Undefined unless \eqn{\nu \ge 2}, in which case the mean is
#'     zero.
#'
#'   **Variance**:
#'
#'   \deqn{
#'     \frac{\nu}{\nu - 2}
#'   }{
#'     \nu / (\nu - 2)
#'   }
#'
#'   Undefined if \eqn{\nu < 1}, infinite when \eqn{1 < \nu \le 2}.
#'
#'   **Probability density function (p.d.f)**:
#'
#'   \deqn{
#'     f(x) = \frac{\Gamma(\frac{\nu + 1}{2})}{\sqrt{\nu \pi} \Gamma(\frac{\nu}{2})} (1 + \frac{x^2}{\nu} )^{- \frac{\nu + 1}{2}}
#'   }{
#'     f(x) = \Gamma((\nu + 1) / 2) / (\sqrt(\nu \pi) \Gamma(\nu / 2)) (1 + x^2 / \nu)^(- (\nu + 1) / 2)
#'   }
#'
#' @seealso [stats::TDist]
#'
#' @examples
#' dist <- dist_smd_ind(df = c(1,2,5), mult = c(0,1,2), sigma = c(1,2,3))
#'
#' dist
#' mean(dist)
#' variance(dist)
#'
#' generate(dist, 10)
#'
#' density(dist, 2)
#' density(dist, 2, log = TRUE)
#'
#' cdf(dist, 4)
#'
#' quantile(dist, 0.7)
#' @import vctrs
#' @name dist_smd_ind
#' @export

dist_smd_ind <- function(df, n1, n2, sigma, ncp){
  if(any(df <= 0)){
    abort("The degrees of freedom parameter of an SMD must be strictly positive.")
  }
  if(any(sigma[!is.na(sigma)] <= 0)){
    abort("The scale (sigma) parameter of an SMD distribution must be strictly positive.")
  }
  distributional::dist_wrap("dist_smd_ind",
                            package = "TOSTER",
                            df = df, n1 = n1, n2 = n2, sigma = sigma, ncp = ncp)
}
#dist_smd_ind(df= 12, n1= 20,n2= 20,sigma = .5, ncp =.25)
#' @export
print.dist_smd_ind <- function(x, ...){
  cat(format(x, ...))
}


#' @export
density.dist_smd_ind <- function(x, at, ...){
  ncp <- x[[1]]$ncp
  sigma <- x[[1]]$sigma
  df <- x[[1]]$df

  stats::dt((at - df)/sigma, df, ncp) / sigma
}


#' @export
quantile.dist_smd_ind <- function(x, p, ...){
  ncp <- x[[1]]$ncp
  df = x[[1]]$df
  n1 = x[[1]]$n1
  n2 = x[[1]]$n2

  stats::qt(p, df, ncp) * sqrt(1/n1+1/n2)
}


#' @export
cdf.dist_smd_ind <- function(x, q, ...){
  ncp <- x[[1]]$ncp
  sigma <- x[[1]]$sigma
  df <- x[[1]]$df

  # sqrt((n1+n2)/(n1*n2))*rt(1000000,df=df,ncp=ncp)
  # sqrt((x[["n1"]]+x[["n2"]])/(x[["n1"]]*x[["n2"]]))

  stats::pt((q)/sigma, df, ncp)
}

#' @export
generate.dist_smd_ind <- function(x, times, ...){
  ncp <- x[[1]]$ncp
  sigma <- x[[1]]$sigma
  df <- x[[1]]$df
  n1 = x[[1]]$n1
  n2 = x[[1]]$n2
  # sqrt((n1+n2)/(n1*n2))*rt(1000000,df=df,ncp=ncp)
  # sqrt((x[["n1"]]+x[["n2"]])/(x[["n1"]]*x[["n2"]]))

  stats::rt(times, df, ncp) * sqrt((n1+n2)/(n1*n2))
}

#' @export
mean.dist_smd_ind <- function(x, ...){
  ncp <- x[[1]]$ncp
  sigma <- x[[1]]$sigma
  df <- x[[1]]$df
  n1 = x[[1]]$n1
  n2 = x[[1]]$n2

  ncp
}

#' @export
variance.dist_smd_ind <- function(x, ...){
  ncp <- x[[1]]$ncp
  sigma <- x[[1]]$sigma
  df <- x[[1]]$df
  n1 = x[[1]]$n1
  n2 = x[[1]]$n2
  if(df <= 1) return(NA_real_)
  if(df <= 2) return(Inf)
  sigma^2
}
