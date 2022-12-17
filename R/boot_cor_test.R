#' @title Bootstrapped TOST with t-tests
#' @description A function for a bootstrap method for TOST with all types of t-tests.
#' @param R number of bootstrap replicates
#' @inheritParams t_TOST
#' @details For details on the calculations in this function see vignette("robustTOST").
#' @return An S3 object of class
#'   \code{"TOSTt"} is returned containing the following slots:
#' \describe{
#'   \item{\code{"TOST"}}{A table of class \code{"data.frame"} containing two-tailed t-test and both one-tailed results.}
#'   \item{\code{"eqb"}}{A table of class \code{"data.frame"} containing equivalence bound settings.}
#'   \item{\code{"effsize"}}{ table of class \code{"data.frame"} containing effect size estimates}
#'   \item{\code{"hypothesis"}}{String stating the hypothesis being tested}
#'   \item{\code{"smd"}}{List containing the results of the standardized mean difference calculations (e.g., Cohen's d).Items include: d (estimate), dlow (lower CI bound), dhigh (upper CI bound), d_df (degrees of freedom for SMD), d_sigma (SE), d_lambda (non-centrality), J (bias correction), smd_label (type of SMD), d_denom (denominator calculation)}
#'   \item{\code{"alpha"}}{Alpha level set for the analysis.}
#'   \item{\code{"method"}}{Type of t-test.}
#'   \item{\code{"decision"}}{List included text regarding the decisions for statistical inference.}
#'   \item{\code{"boot"}}{List containing the bootstrap samples.}
#' }
#' @details The implemented test(s) corresponds to the proposal of Chapter 16 of Efron and Tibshirani (1994).
#'  Returns TOSTt class object with bootstrapped based results.
#'  Please note that the repeated measures "corrected" effect size is not available at this time.
#'
#' For details on the calculations in this function see vignette("robustTOST").
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
    } else if(alternative = "less"){
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
    method <- "Pearson's product-moment correlation"
    z.se <- 1 / sqrt(n_obs - 3)
  }
  if (method == "spearman") {
    method <- "Spearman's rank correlation rho"
    # Spearman # Fieller adjusted
    z.se <- (1.06 / (n_obs - 3)) ^ 0.5
  }
  if (method == "kendall") {
    method <- "Kendall's rank correlation tau"
    # Kendall # Fieller adjusted
    z.se <- (0.437 / (n_obs - 4)) ^ 0.5

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


  names(zstat) = "z"
  attr(cint, "conf.level") <- ci

  # Store as htest
  rval <- list(statistic = z_final, p.value = pvalue,
               conf.int = cint,
               estimate = rfinal,
               null.value = NVAL,
               alternative = alternative,
               method = method,
               data.name = DNAME,
               call = match.call())
  class(rval) <- "htest"
  return(rval)
}

