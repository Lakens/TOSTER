#' @title Bootstrapped t-test
#' @description
#' `r lifecycle::badge('stable')`
#'
#' Performs t-tests with bootstrapped p-values and confidence intervals. This function supports
#' standard hypothesis testing alternatives as well as equivalence and minimal effect testing,
#' all with the familiar `htest` output structure.
#'
#' @section Purpose:
#' Use this function when:
#'   * You need more robust inference than provided by standard t-tests
#'   * Your data don't meet the assumptions of normality or homogeneity
#'   * You want to perform equivalence or minimal effect testing with bootstrap methods
#'   * Sample sizes are small or standard parametric approaches may be unreliable
#'   * You prefer the standard `htest` output format for compatibility with other R functions
#'
#' @inheritParams simple_htest
#' @inheritParams boot_t_TOST
#' @param alternative the alternative hypothesis:
#'     * "two.sided": different from mu (default)
#'     * "less": less than mu
#'     * "greater": greater than mu
#'     * "equivalence": between specified bounds
#'     * "minimal.effect": outside specified bounds
#'
#' @param mu a number or vector specifying the null hypothesis value(s):
#'     * For standard alternatives: a single value (default = 0)
#'     * For equivalence/minimal.effect: two values representing the lower and upper bounds
#'
#' @details
#' This function performs bootstrapped t-tests, providing more robust inference than standard
#' parametric t-tests. It supports one-sample, two-sample (independent), and paired designs,
#' as well as five different alternative hypotheses.
#'
#' The bootstrap procedure follows these steps:
#'   * Calculate the test statistic from the original data
#'   * Generate R bootstrap samples by resampling with replacement
#'   * Calculate the test statistic for each bootstrap sample
#'   * Compute the p-value by comparing the original test statistic to the bootstrap distribution
#'   * Calculate confidence intervals using the specified bootstrap method
#'
#' Three bootstrap confidence interval methods are available:
#'   - *Studentized bootstrap ("stud")*: Accounts for the variability in standard error estimates
#'   - *Basic bootstrap ("basic")*: Uses the empirical distribution of bootstrap estimates
#'   - *Percentile bootstrap ("perc")*: Uses percentiles of the bootstrap distribution directly
#'
#' For different alternatives, the p-values are calculated as follows:
#'   * "two.sided": Proportion of bootstrap statistics at least as extreme as the observed statistic (in either direction), multiplied by 2
#'   * "less": Proportion of bootstrap statistics less than or equal to the observed statistic
#'   * "greater": Proportion of bootstrap statistics greater than or equal to the observed statistic
#'   * "equivalence": Maximum of two one-sided p-values (for lower and upper bounds)
#'   * "minimal.effect": Minimum of two one-sided p-values (for lower and upper bounds)
#'
#'
#' For two-sample tests, the test is of \eqn{\bar x - \bar y} (mean of x minus mean of y).
#' For paired samples, the test is of the difference scores (z),
#' wherein \eqn{z = x - y}, and the test is of \eqn{\bar z} (mean of the difference scores).
#' For one-sample tests, the test is of \eqn{\bar x} (mean of x).
#'
#' Unlike the `t_TOST` function, this function returns a standard `htest` object for
#' compatibility with other R functions, while still providing the benefits of bootstrapping.
#'
#' For detailed information on calculation methods, see `vignette("robustTOST")`.
#'
#' @return A list with class `"htest"` containing the following components:
#'
#'   - "p.value": the bootstrapped p-value for the test.
#'   - "stderr": the bootstrapped standard error.
#'   - "conf.int": a bootstrapped confidence interval for the mean appropriate to the specified alternative hypothesis.
#'   - "estimate": the estimated mean or difference in means.
#'   - "null.value": the specified hypothesized value(s) of the mean or mean difference.
#'   - "alternative": a character string describing the alternative hypothesis.
#'   - "method": a character string indicating what type of bootstrapped t-test was performed.
#'   - "boot": the bootstrap samples of the mean or mean difference.
#'   - "data.name": a character string giving the name(s) of the data.
#'   - "call": the matched call.
#'
#' @examples
#'
#' # Example 1: Basic two-sample test with formula notation
#' data(sleep)
#' result <- boot_t_test(extra ~ group, data = sleep)
#' result  # Standard htest output format
#'
#' # Example 2: One-sample bootstrapped t-test
#' set.seed(123)
#' x <- rnorm(20, mean = 0.5, sd = 1)
#' boot_t_test(x, mu = 0, R = 999) # Using fewer replicates for demonstration
#'
#' # Example 3: Paired samples test with percentile bootstrap CI
#' before <- c(5.1, 4.8, 6.2, 5.7, 6.0, 5.5, 4.9, 5.8)
#' after <- c(5.6, 5.2, 6.7, 6.1, 6.5, 5.8, 5.3, 6.2)
#' boot_t_test(x = before, y = after,
#'             paired = TRUE,
#'             alternative = "less",  # Testing if before < after
#'             boot_ci = "perc",
#'             R = 999)
#'
#' # Example 4: Equivalence testing with bootstrapped t-test
#' # Testing if the effect is within ±0.5 units
#' data(mtcars)
#' boot_t_test(mpg ~ am, data = mtcars,
#'             alternative = "equivalence",
#'             mu = c(-0.5, 0.5),
#'             boot_ci = "stud",
#'             R = 999)
#'
#' # Example 5: Minimal effect testing with bootstrapped t-test
#' # Testing if the effect is outside ±3 units
#' boot_t_test(mpg ~ am, data = mtcars,
#'             alternative = "minimal.effect",
#'             mu = c(-3, 3),
#'             R = 999)
#'
#' @references
#' Efron, B., & Tibshirani, R. J. (1994). An introduction to the bootstrap. CRC press.
#'
#' @family Robust tests
#' @name boot_t_test
#' @export boot_t_test

boot_t_test <- function(x, ...){
  UseMethod("boot_t_test")
}

#' @rdname boot_t_test
#' @method boot_t_test default
#' @export

boot_t_test.default <- function(x,
                                y = NULL,
                                var.equal = FALSE,
                                paired = FALSE,
                                alternative = c("two.sided",
                                                "less",
                                                "greater",
                                                "equivalence",
                                                "minimal.effect"),
                                mu = 0,
                                alpha = 0.05,
                                boot_ci = c("stud","basic","perc"),
                                R = 1999, ...){
  alternative = match.arg(alternative)
  boot_ci = match.arg(boot_ci)

  if(!missing(alpha) && (length(alpha) != 1 || !is.finite(alpha) ||
                              alpha < 0 || alpha > 1)) {
    stop("'alpha' must be a single number between 0 and 1")
  }


  if (!is.null(y)) {
    dname <- paste(deparse(substitute(x)), "and",
                   deparse(substitute(y)))
  }
  else {
    dname <- deparse(substitute(x))
  }

  null_test = simple_htest(x = x,
                           y = y,
                           test = "t.test",
                           var.equal = var.equal,
                           paired = paired,
                           alternative = alternative,
                           mu = mu,
                           alpha = 0.05)
  mu = null_test$null.value
  m_vec <- rep(NA, times=length(R)) # mean difference vector
  m_se_vec <- rep(NA, times=length(R)) # mean difference vector
  if(alternative %in% c("equivalence","minimal.effect")){
    conf.level = 1-alpha*2
  } else {
    conf.level = 1-alpha
  }


  if(!is.null(y)){
    dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
    if (paired) {
      i1 <- y
      i2 <- x
      data <- data.frame(i1 = i1, i2 = i2)
      data <- na.omit(data)
      y <- data$i1
      x <- data$i2
    }
    yok <- !is.na(y)
    xok <- !is.na(x)
    y <- y[yok]

  }else{
    dname <- deparse(substitute(x))
    #if (paired) {
    #  stop("'y' is missing for paired test")
    #}

    xok <- !is.na(x)
    yok <- NULL
  }
  x <- x[xok]
  if(paired && !is.null(y)){
    x <- x - y
    y <- NULL
  }
  nx <- length(x)
  mx <- mean(x)
  vx <- var(x)
  if (is.null(y)) {
    if (nx < 2)
      stop("not enough 'x' observations")
    df <- nx - 1
    stderr <- sqrt(vx/nx)
    if (stderr < 10 * .Machine$double.eps * abs(mx)){
      stop("data are essentially constant")
    }
    tstat <- (mx - mu)/stderr

    method <- if (paired) "Bootstrapped Paired t-test" else "Bootstrapped One Sample t-test"
    #estimate <- setNames(mx, if (paired) "mean of the differences" else "mean of x")
    #x.cent <- x - mx # remove to have an untransformed matrix
    X <- matrix(sample(x, size = nx*R, replace = TRUE), nrow = R)
    MX <- rowMeans(X - mx)
    VX <- rowSums((X - MX) ^ 2) / (nx - 1)
    MZ2 = NA
    VZ2 = NA
    for(i in 1:R){
      zi = X[i,]
      MZ2[i] = mean(zi - mx)
      VZ2[i] <- sum((zi - MZ2[i]) ^ 2) / (nx - 1) #rowSums((X - MX) ^ 2) / (nx - 1)
    }

    STDERR <- sqrt(VX/nx)
    TSTAT <- (MX)/STDERR
    #TSTAT_low <- (MX-low_eqbound)/STDERR
    #TSTAT_high <- (MX-high_eqbound)/STDERR
    EFF <- MX+mx

    for(i in 1:nrow(X)){
      dat = X[i,]

      m_vec[i] <- mean(dat, na.rm=TRUE) # mean difference vector
      m_se_vec[i] <- sd(dat, na.rm = TRUE)/sqrt(length(na.omit(dat)))

    }
  }

  if(!is.null(y) && !paired){
    ny <- length(y)
    if(nx < 1 || (!var.equal && nx < 2))
      stop("not enough 'x' observations")
    if(ny < 1 || (!var.equal && ny < 2))
      stop("not enough 'y' observations")
    if(var.equal && nx + ny < 3)
      stop("not enough observations")
    my <- mean(y)
    vy <- var(y)
    method <- paste("Bootstrapped", paste(if (!var.equal) "Welch", "Two Sample t-test"))
    estimate <- c(mx, my)
    names(estimate) <- c("mean of x", "mean of y")
    if(var.equal){
      df <- nx + ny - 2
      v <- 0
      if (nx > 1){
        v <- v + (nx - 1) * vx
      }

      if (ny > 1){
        v <- v + (ny - 1) * vy
      }

      v <- v/df
      stderr <- sqrt(v * (1/nx + 1/ny))
      z <- c(x, y)
      mz <- mean(z)
      #Z <- matrix(sample(z, size = (nx+ny)*R, replace = TRUE), nrow = R)
      X <- matrix(sample(x, size = nx*R, replace = TRUE), nrow = R)
      Y <- matrix(sample(y, size = ny*R, replace = TRUE), nrow = R)
      MX <- rowMeans(X - mx + mz)
      MY <- rowMeans(Y - my + mz)
      V <- (rowSums((X-MX)^2) + rowSums((Y-MY)^2))/df
      STDERR <- sqrt(V*(1/nx + 1/ny))
      EFF <- (MX+mx) - (MY+my)

      #d_vec <- rep(NA, times=length(R))
      for(i in 1:nrow(X)){
        #dat = Z[i,]
        dat_x = X[i,]#dat[1:nx]
        dat_y = Y[i,]#dat[(nx+1):(nx+ny)]

        m_vec[i] <- mean(dat_x, na.rm=TRUE) - mean(dat_y,na.rm=TRUE)  # mean difference vector
        m_se_vec[i] <- sqrt(sd(dat_x, na.rm=TRUE)^2/length(na.omit(dat_x)) + sd(dat_y, na.rm=TRUE)^2/length(na.omit(dat_y)))

      }
    }else{
      stderrx <- sqrt(vx/nx)
      stderry <- sqrt(vy/ny)
      stderr <- sqrt(stderrx^2 + stderry^2)
      df <- stderr^4/(stderrx^4/(nx - 1) + stderry^4/(ny - 1))
      z <- c(x, y)
      mz <- mean(z)
      x.cent <- x - mx + mz
      y.cent <- y - my + mz
      X <- matrix(sample(x, size = nx*R, replace = TRUE), nrow = R)
      Y <- matrix(sample(y, size = ny*R, replace = TRUE), nrow = R)
      MX <- rowMeans(X - mx + mz)
      MY <- rowMeans(Y - my + mz)
      VX <- rowSums((X-MX)^2)/(nx-1)
      VY <- rowSums((Y-MY)^2)/(ny-1)
      STDERR <- sqrt(VX/nx + VY/ny)
      EFF <- (MX+mx) - (MY+my)

      for(i in 1:nrow(X)){
        #dat = Z[i,]
        dat_x = X[i,]#dat[1:nx]
        dat_y = Y[i,]#dat[(nx+1):(nx+ny)]

        m_vec[i] <- mean(dat_x, na.rm=TRUE) - mean(dat_y,na.rm=TRUE)  # mean difference vector
        m_se_vec[i] <- sqrt(sd(dat_x, na.rm=TRUE)^2/length(na.omit(dat_x)) + sd(dat_y, na.rm=TRUE)^2/length(na.omit(dat_y)))

      }
    }
    if (stderr < 10 * .Machine$double.eps * max(abs(mx), abs(my))){
      stop("data are essentially constant")
    }

    tstat <- (mx - my - mu)/stderr
    # Remember tstat[which.max( abs(tstat) )]
    #TSTAT <- (MX - MY)/STDERR

    TSTAT <- (MX-MY)/STDERR
    #TSTAT_low <- (MX-low_eqbound)/STDERR
    #TSTAT_high <- (MX-high_eqbound)/STDERR
  }

  if(is.null(y)){
    diff = mx
  }else{
    diff = mx-my
  }

  if(alternative %in% c("equivalence", "minimal.effect")){


    tstat_l = (diff-min(mu))/stderr
    tstat_u = (diff-max(mu))/stderr
  }

  if (alternative == "less") {

    boot.pval <- mean(TSTAT < tstat)

    boot.cint <- switch(boot_ci,
                        "stud" = stud(m_vec,
                                      boots_se = m_se_vec,
                                      t0 = diff,
                                      se0 = null_test$stderr,
                                      alpha = alpha*2),
                        "basic" = basic(m_vec, t0 = diff, alpha*2),
                        "perc" = perc(m_vec, alpha*2))

  }

  if(alternative == "greater") {
    boot.pval <- mean(TSTAT > tstat)
    boot.cint <- switch(boot_ci,
                        "stud" = stud(m_vec,
                                      boots_se = m_se_vec,
                                      t0 = diff,
                                      se0 = null_test$stderr,
                                      alpha*2),
                        "basic" = basic(m_vec, t0 = diff, alpha*2),
                        "perc" = perc(m_vec, alpha*2))
  }

  if(alternative == "two.sided"){
    boot.pval <- 2*min(mean(TSTAT <= tstat), mean(TSTAT > tstat))
    boot.cint <- switch(boot_ci,
                        "stud" = stud(m_vec,
                                      boots_se = m_se_vec,
                                      t0 = diff,
                                      se0 = null_test$stderr,
                                      alpha),
                        "basic" = basic(m_vec, t0 = diff, alpha),
                        "perc" = perc(m_vec, alpha))
  }

  if(alternative == "equivalence") {
    p_l = mean(TSTAT > tstat_l)
    p_u = mean(TSTAT < tstat_u)
    boot.pval <- max(p_l, p_u)
    boot.cint <- switch(boot_ci,
                        "stud" = stud(m_vec,
                                      boots_se = m_se_vec,
                                      t0 = diff,
                                      se0 = null_test$stderr,
                                      alpha*2),
                        "basic" = basic(m_vec, t0 = diff, alpha*2),
                        "perc" = perc(m_vec, alpha*2))
  }

  if(alternative == "minimal.effect") {
    p_l = mean(TSTAT < tstat_l)
    p_u = mean(TSTAT > tstat_u)
    boot.pval <- min(p_l,p_u)
    boot.cint <- switch(boot_ci,
                        "stud" = stud(m_vec,
                                      boots_se = m_se_vec,
                                      t0 = diff,
                                      se0 = null_test$stderr,
                                      alpha*2),
                        "basic" = basic(m_vec, t0 = diff, alpha*2),
                        "perc" = perc(m_vec, alpha*2))
  }

  boot.se = sd(m_vec, na.rm = TRUE)
  attr(boot.cint, "conf.level") <- conf.level

  rval = list(
    p.value = boot.pval,
    stderr = boot.se,
    conf.int = boot.cint,
    estimate = null_test$estimate,
    null.value = null_test$null.value,
    alternative = alternative,
    method = method,
    boot = m_vec,
    data.name = null_test$data.name,
    call = match.call()
  )

  class(rval) = "htest"

  return(rval)
}

#' @rdname boot_t_test
#' @method boot_t_test formula
#' @export
#'
boot_t_test.formula <- function (formula, data, subset, na.action, ...){
  if(missing(formula)
     || (length(formula) != 3L)
     || (length(attr(terms(formula[-2L]), "term.labels")) != 1L))
    stop("'formula' missing or incorrect")
  
  # Check for paired argument in ...
  dots <- list(...)
  if(!is.null(dots$paired) && isTRUE(dots$paired))
    stop("cannot use 'paired' in formula method")
  
  m <- match.call(expand.dots = FALSE)
  if(is.matrix(eval(m$data, parent.frame())))
    m$data <- as.data.frame(data)
  ## need stats:: for non-standard evaluation
  m[[1L]] <- quote(stats::model.frame)
  m$... <- NULL
  mf <- eval(m, parent.frame())
  DNAME <- paste(names(mf), collapse = " by ")
  names(mf) <- NULL
  response <- attr(attr(mf, "terms"), "response")
  g <- factor(mf[[-response]])
  if(nlevels(g) != 2L)
    stop("grouping factor must have exactly 2 levels")
  DATA <- setNames(split(mf[[response]], g), c("x", "y"))
  y <- do.call("boot_t_test", c(DATA, list(...)))
  y$data.name <- DNAME
  y
}
