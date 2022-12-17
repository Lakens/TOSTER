#' @title Bootstrapped correlation coefficients
#' @description A function for a bootstrap, percentile, method for correlation coefficients.
#' @inheritParams boot_t_TOST
#' @inheritParams z_cor_test
#' @details This function uses a percentile bootstrap methods for the confidence intervals.
#' The returned p-values are calculated from a re-sampled null distribution (similar to boot_t_TOST).
#'  @return A list with class "htest" containing the following components:
#' \describe{
#'   \item{\code{"statistic"}}{z-score}
#'   \item{\code{"p.value"}}{the p-value of the test.}
#'   \item{\code{"estimate"}}{the estimated measure of association, with name "cor", "tau", or "rho" corresponding to the method employed.}
#'   \item{\code{"null.value"}}{the value of the association measure under the null hypothesis.}
#'   \item{\code{"alternative"}}{character string indicating the alternative hypothesis (the value of the input argument alternative). Possible values are "greater", "less", or "two-sided".}
#'   \item{\code{"method"}}{a character string indicating how the association was measured.}
#'   \item{\code{"data.name"}}{a character string giving the names of the data.}
#'   \item{\code{"call"}}{the matched call.}
#' }
#' @references
#' TBA
#' @section References:
#'
#' Efron, B., & Tibshirani, R. J. (1994). An introduction to the bootstrap. CRC press.
#' @export


boot_cor_test <- function(x,
                          y,
                          alternative = c("two.sided", "less", "greater"),
                          method = c("pearson", "kendall", "spearman"),
                          alpha = 0.05,
                          null = 0,
                          TOST = FALSE,
                          R = 1999,
                          ...) {
  alternative = match.arg(alternative)
  method = match.arg(method)

  nulltest = z_cor_test(
    x=x,
    y=y,
    alternative = alternative,
    method = method,
    alpha = alpha,
    null = null,
    TOST = TOST
  )

  if(TOST && null <=0){
    stop("positive value for null must be supplied if using TOST.")
  }
  if(TOST){
    alternative = "less"
  }

  if(alternative != "two.sided"){
    ci = 1 - alpha*2
    intmult = c(1,1)
  } else {
    ci = 1 - alpha
    if(TOST){
      intmult = c(1,1)
    } else if(alternative == "less"){
      intmult = c(1,NA)
    } else {
      intmult = c(NA,1)
    }
  }
  r_xy = cor(x,y,
             method = method)
  df = data.frame(x=x,
                  y=y)
  df = na.omit(df)
  n_obs = nrow(df)

  z_xy = rho_to_z(r_xy)
  DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
  # get absolute value if TOST
  z_test = ifelse(TOST, abs(z_xy), z_xy)
  znull = rho_to_z(null)
  z_test = z_test-znull
  NVAL = null
  if (method == "pearson") {
    # Pearson # Fisher
    method2 <- "Pearson's product-moment correlation"
    names(NVAL) = "correlation"
    rfinal = c(cor = r_xy)
    z.se <- 1 / sqrt(n_obs - 3)
    cint = cor_to_ci(cor = r_xy, n = n_obs, ci = ci,
                     method = "pearson")
  }
  if (method == "spearman") {
    method2 <- "Spearman's rank correlation rho"
    #  # Fieller adjusted
    rfinal = c(rho = r_xy)
    names(NVAL) = "rho"
    z.se <- (1.06 / (n_obs - 3)) ^ 0.5
    cint = cor_to_ci(cor = r_xy, n = n_obs, ci = ci,
                     method = "spearman",
                     correction = "fieller")
  }
  if (method == "kendall") {
    method2 <- "Kendall's rank correlation tau"
    # # Fieller adjusted
    rfinal = c(tau = r_xy)
    names(NVAL) = "tau"
    z.se <- (0.437 / (n_obs - 4)) ^ 0.5

    cint = cor_to_ci(cor = r_xy, n = n_obs, ci = ci,
                     method = "kendall",
                     correction = "fieller")
  }
  z_test2 = z_test / z.se
  z_stat = z_xy / z.se
  cor.boot = rep(NA, times = length(R)) # corr difference vector
  zxy.boot <- rep(NA, times = length(R)) # z vector
  zstat.boot <- rep(NA, times = length(R)) # z vector
  for (i in 1:R) {
    idx <- sample.int(n_obs, n_obs, replace = TRUE)
    dat_run = df[idx, ]
    cor.boot[i] <- cor(dat_run,
                       method = method)[1,2]
    zxy.boot[i] <- rho_to_z(cor.boot[i])
    zstat.boot[i] <- zxy.boot[i]/z.se
  }

  znull = zstat.boot - z_stat
  #m_vec = append(m_vec, nullTOST$effsize$estimate[1])
  #d_vec = append(d_vec, nullTOST$effsize$estimate[2])
  if(alternative == "two.sided")
    boot.pval <- 2 * min(mean(znull <= z_test2), mean(znull > z_test2))
  if(alternative == "less"){
    boot.pval = mean(znull < z_test2)
  }
  if(alternative == "greater"){
    boot.pval = mean(znull > z_test2)
  }


  boot.se = sd(cor.boot, na.rm =TRUE)
  boot.cint <- quantile(cor.boot, c((1 - ci) / 2, 1 - (1 - ci) / 2))


  names(z_test2) = "z"
  attr(boot.cint, "conf.level") <- ci

  # Store as htest
  rval <- list(statistic = z_test2, p.value = boot.pval,
               conf.int = boot.cint,
               estimate = rfinal,
               null.value = NVAL,
               alternative = alternative,
               method = method2,
               data.name = DNAME,
               boot = list(r = cor.boot,
                           z = zxy.boot,
                           znull = znull),
               call = match.call())
  class(rval) <- "htest"
  return(rval)
}

