#' @title Bootstrapped TOST with t-tests
#' @description
#' `r lifecycle::badge('stable')`
#'
#' Performs equivalence testing using the Two One-Sided Tests (TOST) procedure with bootstrapped t-tests.
#' This provides a robust alternative to traditional TOST when data may not meet all parametric assumptions.
#'
#' @section Purpose:
#' Use this function when:
#'   * You want more robust confidence intervals for your effect sizes
#'   * Sample sizes are small and parametric assumptions may not hold
#'   * You want to avoid relying on asymptotic approximations
#'
#' @inheritParams t_TOST
#' @param glass Option to calculate Glass's delta instead of Cohen's d style SMD ('glass1' uses first group's SD, 'glass2' uses second group's SD).
#' @param mu a number indicating the true value of the mean for the two-tailed test (default = 0).
#' @param R number of bootstrap replications (default = 1999).
#' @param boot_ci method for bootstrap confidence interval calculation: "stud" (studentized, default), "basic" (basic bootstrap), or "perc" (percentile bootstrap).
#' @param subset an optional vector specifying a subset of observations to be used.
#' @param na.action a function indicating what should happen when the data contain NAs.
#' @param ... further arguments to be passed to or from methods.
#'
#' @details
#' The function implements a bootstrap method for TOST as described in Chapter 16 of Efron and Tibshirani (1994).
#' This approach provides a robust alternative to traditional parametric TOST when data distributions may not
#' meet standard assumptions.
#'
#' The bootstrap procedure follows these steps:
#'   * Resample with replacement from the original data to create R bootstrap samples
#'   * For each bootstrap sample, calculate test statistics and effect sizes
#'   * Use the distribution of bootstrap results to compute p-values and confidence intervals
#'   * Combine results using the specified bootstrap confidence interval method
#'
#' Three types of bootstrap confidence intervals are available:
#'   * Studentized ("stud"): Accounts for the variability in the standard error estimate
#'   * Basic/Empirical ("basic"): Uses the empirical distribution of bootstrap estimates
#'   * Percentile ("perc"): Uses percentiles of the bootstrap distribution
#'
#' For two-sample tests, the test is of \eqn{\bar x - \bar y} (mean of x minus mean of y).
#' For paired samples, the test is of the difference scores (z),
#' wherein \eqn{z = x - y}, and the test is of \eqn{\bar z} (mean of the difference scores).
#' For one-sample tests, the test is of \eqn{\bar x} (mean of x).
#'
#'
#' For details on the calculations in this function see `vignette("robustTOST")`.
#'
#' @return An S3 object of class `"TOSTt"` is returned containing the following slots:
#'
#'   - "TOST": A table of class "data.frame" containing two-tailed t-test and both one-tailed results.
#'   - "eqb": A table of class "data.frame" containing equivalence bound settings.
#'   - "effsize": Table of class "data.frame" containing effect size estimates.
#'   - "hypothesis": String stating the hypothesis being tested.
#'   - "smd": List containing the results of the standardized mean difference calculations (e.g., Cohen's d).
#'      - Items include: d (estimate), dlow (lower CI bound), dhigh (upper CI bound), d_df (degrees of freedom for SMD), d_sigma (SE), d_lambda (non-centrality), J (bias correction), smd_label (type of SMD), d_denom (denominator calculation).
#'   - "alpha": Alpha level set for the analysis.
#'   - "method": Type of t-test.
#'   - "decision": List included text regarding the decisions for statistical inference.
#'   - "boot": List containing the bootstrap samples for SMD and raw effect sizes.
#'
#' @examples
#' \dontrun{
#' # Example 1: Two-Sample Test with Symmetric Bounds
#' set.seed(1234)
#' group1 <- rnorm(30, mean = 5, sd = 2)
#' group2 <- rnorm(30, mean = 5.5, sd = 2.2)
#'
#' # Using symmetric bounds of ±1.5
#' result <- boot_t_TOST(x = group1,
#'                      y = group2,
#'                      eqb = 1.5,
#'                      R = 999)  # Using fewer replications for demonstration
#'
#' # Example 2: Paired Sample Test with Percentile Bootstrap
#' set.seed(5678)
#' pre <- rnorm(25, mean = 100, sd = 15)
#' post <- pre + rnorm(25, mean = 3, sd = 10)
#'
#' result <- boot_t_TOST(x = pre,
#'                      y = post,
#'                      paired = TRUE,
#'                      eqb = c(-5, 8),  # Asymmetric bounds
#'                      boot_ci = "perc")
#'
#' # Example 3: One Sample Test
#' set.seed(9101)
#' scores <- rnorm(40, mean = 0.3, sd = 1)
#'
#' # Testing if mean is equivalent to zero within ±0.5 units
#' result <- boot_t_TOST(x = scores,
#'                      eqb = 0.5,
#'                      boot_ci = "basic")
#'}
#' @references
#' Efron, B., & Tibshirani, R. J. (1994). An introduction to the bootstrap. CRC press.
#'
#' @family Robust tests
#' @family TOST
#' @importFrom stats var
#' @name boot_t_TOST
#' @export boot_t_TOST

boot_t_TOST <- function(x, ...){
  UseMethod("boot_t_TOST")
}

#' @rdname boot_t_TOST
#' @method boot_t_TOST default
#' @export

boot_t_TOST.default <- function(x,
                                y = NULL,
                                hypothesis = "EQU",
                                paired = FALSE,
                                var.equal = FALSE,
                                eqb,
                                low_eqbound,
                                high_eqbound,
                                eqbound_type = "raw",
                                alpha = 0.05,
                                bias_correction = TRUE,
                                rm_correction = FALSE,
                                glass = NULL,
                                mu = 0,
                                R = 1999,
                                boot_ci = c("stud", "basic", "perc"),
                                ...){
  boot_ci = match.arg(boot_ci)

  if(!missing(mu) && (length(mu) != 1 || is.na(mu))) {
    stop("'mu' must be a single number")
  }


  if(!missing(alpha) && (length(alpha) != 1 || !is.finite(alpha) ||
                              alpha < 0 || alpha > 1)) {
    stop("'alpha' must be a single number between 0 and 1")
  }

  if(!missing(eqb)){
    if(!is.numeric(eqb) || length(eqb) > 2){
      stop(
        "eqb must be a numeric of a length of 1 or 2"
      )
    }
    if(length(eqb) == 1){
      high_eqbound = abs(eqb)
      low_eqbound = -1*abs(eqb)
    } else {
      high_eqbound = max(eqb)
      low_eqbound = min(eqb)
    }


  }

  interval_no_zero = test_interval_no_zero(c(low_eqbound, high_eqbound))

  if(interval_no_zero){
    message("Equivalence interval does not include zero.")
  }

  if (!is.null(y)) {
    dname <- paste(deparse(substitute(x)), "and",
                   deparse(substitute(y)))
  }
  else {
    dname <- deparse(substitute(x))
  }

  if(!is.null(y)){
    nullTOST = t_TOST(x = x,
                      y = y,
                      hypothesis = hypothesis,
                      paired = paired,
                      var.equal = var.equal,
                      low_eqbound = low_eqbound,
                      high_eqbound = high_eqbound,
                      eqbound_type = eqbound_type,
                      alpha = alpha,
                      mu = mu,
                      bias_correction = bias_correction,
                      smd_ci = "z")
  } else{
    nullTOST = t_TOST(x = x,
                      hypothesis = hypothesis,
                      paired = FALSE,
                      var.equal = TRUE,
                      low_eqbound = low_eqbound,
                      high_eqbound = high_eqbound,
                      eqbound_type = eqbound_type,
                      alpha = alpha,
                      mu = mu,
                      bias_correction = bias_correction,
                      rm_correction = FALSE,
                      smd_ci = "z")
  }
  d_vec <- rep(NA, times=length(R)) # smd vector
  m_vec <- rep(NA, times=length(R)) # mean difference vector
  d_se_vec <- rep(NA, times=length(R)) # smd vector SE
  m_se_vec <- rep(NA, times=length(R)) # mean difference vector SE
  #t_vec <- rep(NA, times=length(R)) # t-test vector
  #tl_vec <- rep(NA, times=length(R)) # lower bound vector
  #tu_vec <- rep(NA, times=length(R)) # upper bound vector

  conf.level = 1-alpha*2

  if(bias_correction){
    smd_type = 'g'
  } else {
    smd_type = 'd'
  }

  if(!is.null(y)){
    dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
    if (paired) {
      i1 <- x
      i2 <- y
      data <- data.frame(i1 = i1, i2 = i2)
      data <- na.omit(data)
      x <- data$i1
      y <- data$i2
    }
    yok <- !is.na(y)
    xok <- !is.na(x)
    y <- y[yok]

  }else{
    # One sample ----
    dname <- deparse(substitute(x))
    #if (paired) {
    #  stop("'y' is missing for paired test")
    #}

    xok <- !is.na(x)
    yok <- NULL
  }
  x <- x[xok]
  #if(paired && !is.null(y)){
  #  x <- x - y
  #  y <- NULL
  #}
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
    #tstat <- (mx - mu)/stderr
    #tstat_low = (mx - low_eqbound)/stderr
    #tstat_high = (mx - high_eqbound)/stderr

    method <-  "Bootstrapped One Sample t-test"
    #estimate <- setNames(mx, if (paired) "mean of the differences" else "mean of x")
    #x.cent <- x - mx # remove to have an untransformed matrix
    X <- matrix(sample(x, size = nx*R, replace = TRUE), nrow = R)
    MX <- rowMeans(X - mx)
    VX <- rowSums((X - MX) ^ 2) / (nx - 1)
    STDERR <- sqrt(VX/nx)
    TSTAT <- (MX)/STDERR
    #TSTAT_low <- (MX-low_eqbound)/STDERR
    #TSTAT_high <- (MX-high_eqbound)/STDERR
    EFF <- MX+mx

    for(i in 1:nrow(X)){
      dat = X[i,]
      runTOST =  t_TOST(x = dat,
                         hypothesis = hypothesis,
                         paired = FALSE,
                         var.equal = TRUE,
                         low_eqbound = low_eqbound,
                         high_eqbound = high_eqbound,
                         eqbound_type = eqbound_type,
                         alpha = alpha,
                         mu = mu,
                         bias_correction = bias_correction,
                         rm_correction = FALSE,
                        smd_ci = "z")

      d_vec[i] <- runTOST$smd$d # smd vector
      m_vec[i] <- runTOST$effsize$estimate[1] # mean difference vector
      d_se_vec[i] <- runTOST$effsize$SE[2] # smd vector
      m_se_vec[i] <- runTOST$effsize$SE[1] # mean difference vector
      #t_vec[i] <- runTOST$TOST$t[1] - mx # t-test vector
      #tl_vec[i] <- runTOST$TOST$t[2] - mx # lower bound vector
      #tu_vec[i] <- runTOST$TOST$t[3] - mx # upper bound vector
    }
  }
  # paired -----
  if (paired){

    ny <- length(y)
    if(nx < 1 || (!var.equal && nx < 2))
      stop("not enough 'x' observations")
    if(ny < 1 || (!var.equal && ny < 2))
      stop("not enough 'y' observations")
    if(var.equal && nx + ny < 3)
      stop("not enough observations")
    my <- mean(y)
    vy <- var(y)

    z <- x - y

    nz <- length(z)
    mz <- mean(z)
    vz <- var(z)

    if (nz < 2)
      stop("not enough 'x' observations")
    df <- nz - 1
    stderr <- sqrt(vz/nz)
    if (stderr < 10 * .Machine$double.eps * abs(mz)){
      stop("data are essentially constant")
    }

    method <- "Bootstrapped Paired t-test"
    #estimate <- setNames(mx, if (paired) "mean of the differences" else "mean of x")
    #x.cent <- x - mx # remove to have an untransformed matrix
    #Z <- matrix(sample(z, size = nz*R, replace = TRUE), nrow = R)
    MZ <- rep(NA, times=length(R)) # Means
    VZ <- rep(NA, times=length(R)) # Variance
    STDERR <- rep(NA, times=length(R))
    TSTAT <- rep(NA, times=length(R))
    EFF <- rep(NA, times=length(R))
    #VZ <- rowSums((Z - MZ) ^ 2) / (nz - 1)
    #STDERR <- sqrt(VZ/nz)
    #TSTAT <- (MZ)/STDERR
    #TSTAT_low <- (MX-low_eqbound)/STDERR
    #TSTAT_high <- (MX-high_eqbound)/STDERR
    #EFF <- MZ+mz

    for(i in 1:R){
      sampler = sample(1:nrow(data), replace = TRUE)
      zi = data$i1[sampler]-data$i2[sampler]
      runTOST = t_TOST(x = data$i1[sampler],
                          y = data$i2[sampler],
                        hypothesis = hypothesis,
                        paired = TRUE,
                        var.equal = FALSE,
                        low_eqbound = low_eqbound,
                        high_eqbound = high_eqbound,
                        eqbound_type = eqbound_type,
                        alpha = alpha,
                        mu = mu,
                        bias_correction = bias_correction,
                       glass = glass,
                        rm_correction = rm_correction,
                       smd_ci = "z")
      MZ[i] = mean(zi - mz)
      VZ[i] <- sum((zi - MZ[i]) ^ 2) / (nz - 1) #rowSums((X - MX) ^ 2) / (nx - 1)
      STDERR[i] <- sqrt(VZ[i]/nz)
      TSTAT[i] <- MZ[i]/STDERR[i]
      EFF[i] <- MZ[i] + mz
      d_vec[i] <- runTOST$smd$d # smd vector
      m_vec[i] <- runTOST$effsize$estimate[1] # mean difference vector
      d_se_vec[i] <- runTOST$effsize$SE[2] # smd vector
      m_se_vec[i] <- runTOST$effsize$SE[1] # mean difference vector
    }



  }

  # two sample -----
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
      ## var equal true ----
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
        runTOST =  t_TOST(x = dat_x,
                          y = dat_y,
                          hypothesis = hypothesis,
                          paired = paired,
                          var.equal = var.equal,
                          low_eqbound = low_eqbound,
                          high_eqbound = high_eqbound,
                          eqbound_type = eqbound_type,
                          alpha = alpha,
                          mu = mu,
                          bias_correction = bias_correction,
                          rm_correction = FALSE,
                          smd_ci = "z")

        d_vec[i] <- runTOST$smd$d # smd vector
        m_vec[i] <- runTOST$effsize$estimate[1] # mean difference vector
        d_se_vec[i] <- runTOST$effsize$SE[2] # smd vector
        m_se_vec[i] <- runTOST$effsize$SE[1] # mean difference vector
        #t_vec[i] <- runTOST$TOST$t[1] # t-test vector
        #tl_vec[i] <- runTOST$TOST$t[2] # lower bound vector
        #tu_vec[i] <- runTOST$TOST$t[3] # upper bound vector
      }
    }else{
      ## welch -----
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
        runTOST =  t_TOST(x = dat_x,
                          y = dat_y,
                          hypothesis = hypothesis,
                          paired = paired,
                          var.equal = var.equal,
                          low_eqbound = low_eqbound,
                          high_eqbound = high_eqbound,
                          eqbound_type = eqbound_type,
                          alpha = alpha,
                          mu = mu,
                          bias_correction = bias_correction,
                          rm_correction = FALSE)

        d_vec[i] <- runTOST$smd$d # smd vector
        m_vec[i] <- runTOST$effsize$estimate[1] # mean difference vector
        d_se_vec[i] <- runTOST$effsize$SE[2] # smd vector
        m_se_vec[i] <- runTOST$effsize$SE[1] # mean difference vector
        #t_vec[i] <- runTOST$TOST$t[1] # t-test vector
        #tl_vec[i] <- runTOST$TOST$t[2] # lower bound vector
        #tu_vec[i] <- runTOST$TOST$t[3] # upper bound vector
      }
    }
    if (stderr < 10 * .Machine$double.eps * max(abs(mx), abs(my))){
      stop("data are essentially constant")
    }

    tstat <- (mx - my - mu)/stderr
    #TSTAT <- (MX - MY)/STDERR

    TSTAT <- (MX-MY)/STDERR
    #TSTAT_low <- (MX-low_eqbound)/STDERR
    #TSTAT_high <- (MX-high_eqbound)/STDERR
  }
  tstat = nullTOST$TOST$t[1]
  tstat_l = nullTOST$TOST$t[2]
  tstat_u = nullTOST$TOST$t[3]
  #m_vec = append(m_vec, nullTOST$effsize$estimate[1])
  #d_vec = append(d_vec, nullTOST$effsize$estimate[2])

  boot.pval <- 2 * min(mean(TSTAT <= tstat), mean(TSTAT > tstat))

  if(hypothesis == "EQU"){
    p_l = mean(TSTAT > tstat_l)
    p_u = mean(TSTAT < tstat_u)
  } else{
    p_l = mean(TSTAT < tstat_l)
    p_u = mean(TSTAT > tstat_u)
  }

  boot.se = sd(m_vec)
  boot.cint <- switch(boot_ci,
                      "stud" = stud(m_vec,
                                    boots_se = m_se_vec,
                                    t0 = nullTOST$effsize$estimate[1],
                                    se0 = nullTOST$effsize$SE[1],
                                    alpha*2),
                      "basic" = basic(m_vec, t0 = nullTOST$effsize$estimate[1], alpha*2),
                      "perc" = perc(m_vec, alpha*2))
  d.cint <- switch(boot_ci,
                   "stud" = stud(d_vec,
                                 boots_se = d_se_vec,
                                 t0 = nullTOST$effsize$estimate[2],
                                 se0 = nullTOST$effsize$SE[2],
                                 alpha*2),
                   "basic" = basic(d_vec,t0 = nullTOST$effsize$estimate[2], alpha*2),
                   "perc" = perc(d_vec, alpha*2))
  d.se = sd(d_vec)

  TOST = nullTOST$TOST
  TOST$p.value = c(boot.pval, p_l, p_u)
  effsize = nullTOST$effsize
  effsize$SE = c(boot.se,d.se)
  effsize$lower.ci = c(boot.cint[1],
                       d.cint[1])

  effsize$upper.ci = c(boot.cint[2],
                       d.cint[2])
  pTOST = max(p_l,p_u)
  TOSToutcome<-ifelse(pTOST<alpha,"significant","non-significant")
  testoutcome<-ifelse(boot.pval<alpha,"significant","non-significant")
  if(hypothesis == "EQU"){
    pTOST = max(p_l,
                p_u) # get highest p value for TOST result
    tTOST = ifelse(abs(tstat_l) < abs(tstat_u),
                   tstat_l,
                   tstat_u) #Get lowest t-value for summary TOST result
  } else {
    pTOST = min(p_l,
                p_u) # get highest p value for TOST result
    tTOST = ifelse(abs(tstat_l) > abs(tstat_u),
                   tstat_l,
                   tstat_u) #Get lowest t-value for summary TOST result

    if(!interval_no_zero){
      if(pTOST <= boot.pval){
        message("MET test may have higher error rates than a nil two-tailed test. Consider wider equivalence bounds.")
      }
    }
  }

  # Change text based on two tailed t test if mu is not zero
  if(mu == 0){
    mu_text = "zero"
  } else {
    mu_text = mu
  }

  if(hypothesis == "EQU"){
    #format(low_eqbound, digits = 3, nsmall = 3, scientific = FALSE)
    TOST_restext = paste0("The equivalence test was ",
                          TOSToutcome,", t(",round(df, digits=2),") = ",
                          format(tTOST, digits = 3,
                                 nsmall = 3, scientific = FALSE),", p = ",
                          format(pTOST, digits = 3,
                                 nsmall = 3, scientific = TRUE),sep="")
  } else {
    TOST_restext = paste0("The minimal effect test was ",
                          TOSToutcome,", t(",round(df, digits=2),") = ",
                          format(tTOST, digits = 3,
                                 nsmall = 3, scientific = FALSE),", p = ",
                          format(pTOST, digits = 3,
                                 nsmall = 3, scientific = TRUE),sep="")
  }

  ttest_restext = paste0("The null hypothesis test was ",
                         testoutcome,", t(",round(df, digits=2),") = ",
                         format(tstat, digits = 3,
                                nsmall = 3, scientific = FALSE),", p = ",
                         format(boot.pval, digits = 3,
                                nsmall = 3, scientific = TRUE),sep="")
  combined_outcome = tost_decision(hypothesis = hypothesis,
                                    alpha = alpha,
                                    pvalue = boot.pval,
                                    pTOST = pTOST,
                                    mu_text = mu_text)


  decision = list(
    TOST = TOST_restext,
    ttest = ttest_restext,
    combined = combined_outcome
  )
  call2 = match.call()
  call2$boot_ci = boot_ci
  rval = list(
    TOST = TOST,
    eqb = nullTOST$eqb,
    alpha = alpha,
    method = method,
    hypothesis = nullTOST$hypothesis,
    effsize = effsize,
    smd = nullTOST$smd,
    decision = decision,
    boot = list(SMD = d_vec,
                raw = m_vec),
    data.name = dname,
    call = call2
  )

  class(rval) = "TOSTt"

  return(rval)
}

#' @rdname boot_t_TOST
#' @method boot_t_TOST formula
#' @export
#'
boot_t_TOST.formula <- function (formula, data, subset, na.action, ...){
  if(missing(formula)
     || (length(formula) != 3L)
     || (length(attr(terms(formula[-2L]), "term.labels")) != 1L))
    stop("'formula' missing or incorrect")
  
  # Check for paired argument in ... and warn user
  dots <- list(...)
  if("paired" %in% names(dots)){
    if(isTRUE(dots$paired)){
      message("Using 'paired = TRUE' with the formula interface is not recommended. Please ensure your data is sorted appropriately to make the correct paired comparison.")
    }
  }
  
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
  y <- do.call("boot_t_TOST", c(DATA, list(...)))
  y$data.name <- DNAME
  y
}
