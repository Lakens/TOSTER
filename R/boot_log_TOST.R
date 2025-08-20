#' @title Bootstrapped TOST with Log Transformed t-tests
#' @description
#' `r lifecycle::badge('stable')`
#'
#' Performs equivalence testing using the Two One-Sided Tests (TOST) procedure with bootstrapped
#' log-transformed t-tests. This approach is particularly useful for ratio-scale data where the
#' equivalence bounds are expressed as ratios (e.g., bioequivalence studies).
#'
#' @section Purpose:
#' Use this function when:
#'   - Your data is on a ratio scale (all values must be positive)
#'   - You want to establish equivalence based on the ratio of means rather than their difference
#'   - Traditional parametric methods may not be appropriate due to skewed distributions
#'   - You need to analyze bioequivalence data where bounds are expressed as ratios
#'
#' @inheritParams boot_t_TOST
#' @inheritParams log_TOST
#' @param x a (non-empty) numeric vector of positive data values on a ratio scale.
#' @param y an optional (non-empty) numeric vector of positive data values on a ratio scale.
#' @param hypothesis 'EQU' for equivalence (default), or 'MET' for minimal effects test.
#' @param paired a logical indicating whether you want a paired t-test.
#' @param var.equal a logical variable indicating whether to treat the two variances as being equal.
#' @param eqb Equivalence bound expressed as a ratio. Can provide 1 value (e.g., 1.25 for bounds of 0.8 and 1.25)
#'   or 2 specific values that represent the lower and upper equivalence bounds (e.g., c(0.8, 1.25)).
#' @param alpha alpha level (default = 0.05).
#' @param null the ratio value under the null hypothesis (default = 1).
#' @param boot_ci method for bootstrap confidence interval calculation: "stud" (studentized, default),
#'   "basic" (basic bootstrap), or "perc" (percentile bootstrap).
#' @param R number of bootstrap replications (default = 1999).
#' @param ... further arguments to be passed to or from methods.
#'
#' @details
#' The function implements a bootstrap method for log-transformed TOST as recommended by
#' He et al. (2022) and corresponds to the proposal in Chapter 16 of Efron and Tibshirani (1994).
#' This is approximately equivalent to the percentile bootstrap method mentioned by He et al. (2014).
#'
#' For two-sample tests, the test is of \eqn{\bar log(x) - \bar log(y)} which corresponds to testing
#' the ratio of geometric means. For paired samples, the test is of difference scores on the log scale,
#' \eqn{z = log(x) - log(y) = log(x/y)}, which also corresponds to a ratio test.
#'
#'
#' The bootstrap procedure follows these steps:
#'   - Log-transform the data
#'   - Perform resampling with replacement to generate bootstrap samples
#'   - For each bootstrap sample, calculate test statistics and effect sizes
#'   - Use the distribution of bootstrap results to compute p-values and confidence intervals
#'   - Back-transform for the ratio of means
#'
#'
#' Note that all input data must be positive (ratio scale with a true zero) since log transformation
#' is applied. The function will stop with an error if any negative values are detected.
#'
#' For details on the calculations in this function see `vignette("robustTOST")`.
#'
#' @return An S3 object of class `"TOSTt"` is returned containing the following slots:
#'
#'   - "TOST": A table of class `"data.frame"` containing two-tailed t-test and both one-tailed results.
#'   - "eqb": A table of class `"data.frame"` containing equivalence bound settings.
#'   - "effsize": Table of class `"data.frame"` containing effect size estimates.
#'   - "hypothesis": String stating the hypothesis being tested.
#'   - "smd": List containing the results of the means ratio calculation.
#'      - Items include: d (means ratio estimate), dlow (lower CI bound), dhigh (upper CI bound), d_df (degrees of freedom for SMD), d_sigma (SE), d_lambda (non-centrality), J (bias correction), smd_label (type of SMD), d_denom (denominator calculation).
#'   - "alpha": Alpha level set for the analysis.
#'   - "method": Type of t-test.
#'   - "decision": List included text regarding the decisions for statistical inference.
#'   - "boot": List containing the bootstrap samples.
#'
#' @examples
#' # Example 1: Two-Sample Test for Bioequivalence
#' # Generate ratio scale data (e.g., drug concentrations)
#' test_group <- rlnorm(30, meanlog = 3.5, sdlog = 0.4)
#' ref_group <- rlnorm(30, meanlog = 3.6, sdlog = 0.4)
#'
#' # FDA standard bioequivalence bounds (80% to 125%)
#' result <- boot_log_TOST(x = test_group,
#'                         y = ref_group,
#'                         eqb = 1.25,  # Creates bounds of 0.8 and 1.25
#'                         R = 999)     # Reduce for demonstration
#'
#' # Example 2: Paired Sample Test
#' # Generate paired ratio scale data
#' n <- 20
#' baseline <- rlnorm(n, meanlog = 4, sdlog = 0.3)
#' followup <- baseline * rlnorm(n, meanlog = 0.05, sdlog = 0.2)
#'
#' # Test with asymmetric bounds
#' result <- boot_log_TOST(x = followup,
#'                         y = baseline,
#'                         paired = TRUE,
#'                         eqb = c(0.85, 1.20),
#'                         boot_ci = "perc")
#'
#' @references
#' Efron, B., & Tibshirani, R. J. (1994). An introduction to the bootstrap. CRC press.
#'
#' He, Y., Deng, Y., You, C., & Zhou, X. H. (2022). Equivalence tests for ratio of means in bioequivalence studies under crossover design. Statistical Methods in Medical Research, 09622802221093721.
#'
#' Food and Drug Administration (2014). Bioavailability and Bioequivalence Studies Submitted in NDAs or INDs — General Considerations.
#' Center for Drug Evaluation and Research. Docket: FDA-2014-D-0204.
#' https://www.fda.gov/regulatory-information/search-fda-guidance-documents/bioavailability-and-bioequivalence-studies-submitted-ndas-or-inds-general-considerations
#'
#' @name boot_log_TOST
#' @family Robust tests
#' @family TOST
#' @export boot_log_TOST

boot_log_TOST <- function(x, ...){
  UseMethod("boot_log_TOST")
}

#' @rdname boot_log_TOST
#' @importFrom stats var
#' @method boot_log_TOST default
#' @export

boot_log_TOST.default <- function(x,
                                y = NULL,
                                hypothesis = c("EQU","MET"),
                                paired = FALSE,
                                var.equal = FALSE,
                                eqb = 1.25,
                                alpha = 0.05,
                                null = 1,
                                boot_ci = c("stud","basic", "perc"),
                                R = 1999, ...){
  hypothesis = match.arg(hypothesis)
  boot_ci = match.arg(boot_ci)
  if(!missing(null) && (length(null) != 1 || is.na(null))) {
    stop("'null' must be a single number")
  }
  if(any(x < 0) || any(y < 0)){
    stop("Negative values detected. Values must be on ratio scale (true zero).")
  }

  if(!missing(alpha) && (length(alpha) != 1 || !is.finite(alpha) ||
                              alpha < 0 || alpha > 1)) {
    stop("'alpha' must be a single number between 0 and 1")
  }


    if(!is.numeric(eqb) || length(eqb) > 2){
      stop(
        "eqb must be a numeric of a length of 1 or 2"
      )
    }
    if(length(eqb) == 1){
      if(eqb > 1){
        high_eqbound = eqb
        low_eqbound = 1/eqb
      } else{
        high_eqbound = 1/eqb
        low_eqbound = eqb
      }

    } else {
      high_eqbound = max(eqb)
      low_eqbound = min(eqb)


    }
  interval_no_zero = test_interval_no_zero(c(log(low_eqbound), log(high_eqbound)))

  if(interval_no_zero){
    message("Equivalence interval does not include 1.")
  }

  if (!is.null(y)) {
    dname <- paste(deparse(substitute(x)), "and",
                   deparse(substitute(y)))
    x = log(x)
    y = log(y)
  }
  else {
    stop("One sample tests not supported at this time.")
  }


  #if (paired) {
  #  stop("'y' is missing for paired test")
  #}

  xok <- !is.na(x)
  yok <- NULL

x <- x[xok]
if(paired){
  x <- x - y
  y <- NULL
}
nx <- length(x)
mx <- mean(x)
vx <- var(x)

if(!paired){
  nullTOST = log_TOST(x = exp(x),
                      y = exp(y),
                      hypothesis = hypothesis,
                      paired = paired,
                      var.equal = var.equal,
                      eqb = eqb,
                      alpha = alpha,
                      null = null)
} else {
  nullTOST = log_pair(x = x,
                      hypothesis = hypothesis,
                      eqb = eqb,
                      alpha = alpha,
                      null = null)

}


  d_vec <- rep(NA, times=length(R)) # smd vector
  m_vec <- rep(NA, times=length(R)) # mean difference vector
  se_vec <- NA # Standard error vector
  #t_vec <- rep(NA, times=length(R)) # t-test vector
  #tl_vec <- rep(NA, times=length(R)) # lower bound vector
  #tu_vec <- rep(NA, times=length(R)) # upper bound vector

  conf.level = 1-alpha*2

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

  }
  x <- x[xok]

  nx <- length(x)
  mx <- mean(x)
  vx <- var(x)
  # Paired ----
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
    method <-  "Bootstrapped Log Paired t-test" # else "Bootstrapped One Sample t-test"
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
      runTOST =  log_pair(
        x = dat,
        hypothesis = hypothesis,
        eqb = eqb,
        alpha = alpha,
        null = null
      )

      d_vec[i] <- runTOST$smd$d # smd vector
      m_vec[i] <- runTOST$effsize$estimate[1] # mean difference vector
      se_vec[i] = runTOST$effsize$SE[1] # SE
      #t_vec[i] <- runTOST$TOST$t[1] - mx # t-test vector
      #tl_vec[i] <- runTOST$TOST$t[2] - mx # lower bound vector
      #tu_vec[i] <- runTOST$TOST$t[3] - mx # upper bound vector
    }
  }
# Two sample ----
  if(!is.null(y)){
    ny <- length(y)
    if(nx < 1 || (!var.equal && nx < 2))
      stop("not enough 'x' observations")
    if(ny < 1 || (!var.equal && ny < 2))
      stop("not enough 'y' observations")
    if(var.equal && nx + ny < 3)
      stop("not enough observations")
    my <- mean(y)
    vy <- var(y)
    method <- paste("Bootstrapped Log", paste(if (!var.equal) "Welch", "Two Sample t-test"))
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
        runTOST =  log_TOST(x = exp(dat_x),
                            y = exp(dat_y),
                            hypothesis = hypothesis,
                            paired = paired,
                            var.equal = var.equal,
                            eqb = eqb,
                            alpha = alpha,
                            null = null)

        d_vec[i] <- runTOST$smd$d # smd vector
        m_vec[i] <- runTOST$effsize$estimate[1] # mean difference vector
        se_vec[i] = runTOST$effsize$SE[1] # SE
        #t_vec[i] <- runTOST$TOST$t[1] # t-test vector
        #tl_vec[i] <- runTOST$TOST$t[2] # lower bound vector
        #tu_vec[i] <- runTOST$TOST$t[3] # upper bound vector
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
        runTOST =  log_TOST(x = exp(dat_x),
                          y = exp(dat_y),
                          hypothesis = hypothesis,
                          paired = paired,
                          var.equal = var.equal,
                          eqb = eqb,
                          alpha = alpha,
                          null = null)

        d_vec[i] <- runTOST$smd$d # smd vector
        m_vec[i] <- runTOST$effsize$estimate[1] # mean difference vector
        se_vec[i] = runTOST$effsize$SE[1] # SE
        #t_vec[i] <- runTOST$TOST$t[1] # t-test vector
        #tl_vec[i] <- runTOST$TOST$t[2] # lower bound vector
        #tu_vec[i] <- runTOST$TOST$t[3] # upper bound vector
      }
    }
    if (stderr < 10 * .Machine$double.eps * max(abs(mx), abs(my))){
      stop("data are essentially constant")
    }

    tstat <- (mx - my - null)/stderr
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
  d.se = sd(exp(m_vec))
  boot.cint <- switch(boot_ci,
                      "stud" = stud(boots_est = m_vec, boots_se = se_vec,
                                    se0=nullTOST$effsize$SE[1], t0 = nullTOST$effsize$estimate[1],
                                    alpha),
                      "basic" = basic(m_vec, t0 = nullTOST$effsize$estimate[1], alpha*2),
                      "perc" = perc(m_vec, alpha*2))
  d.cint = exp(boot.cint)
  #d.cint <- switch(boot_ci,
  #                 "basic" = basic(d_vec, t0 = nullTOST$effsize$estimate[2], alpha*2),
  #                 "perc" = perc(d_vec, alpha*2))
  #d.se = sd(d_vec)

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
  if(null == 1){
    mu_text = "1"
  } else {
    mu_text = null
  }

  if(hypothesis == "EQU"){
    #format(low_eqbound, digits = 3, nsmall = 3, scientific = FALSE)
    TOST_restext = paste0("The equivalence test was ",TOSToutcome,", t(",round(df, digits=2),") = ",format(tTOST, digits = 3, nsmall = 3, scientific = FALSE),", p = ",format(pTOST, digits = 3, nsmall = 3, scientific = TRUE),sep="")
  } else {
    TOST_restext = paste0("The minimal effect test was ",
                          TOSToutcome,", t(",round(df, digits=2),") = ",
                          format(tTOST, digits = 3, nsmall = 3,
                                 scientific = FALSE),", p = ",
                          format(pTOST, digits = 3, nsmall = 3,
                                 scientific = TRUE),
                          sep="")
  }

  ttest_restext = paste0("The null hypothesis test was ",
                         testoutcome,", t(",round(df, digits=2),") = ",
                         format(tstat, digits = 3, nsmall = 3,
                                scientific = FALSE),", p = ",
                         format(boot.pval, digits = 3, nsmall = 3,
                                scientific = TRUE),sep="")

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

#' @rdname boot_log_TOST
#' @method boot_log_TOST formula
#' @export
#'
boot_log_TOST.formula <- function (formula, data, subset, na.action, ...){
  if(missing(formula)
     || (length(formula) != 3L)
     || (length(attr(terms(formula[-2L]), "term.labels")) != 1L))
    stop("'formula' missing or incorrect")
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
  y <- do.call("boot_log_TOST", c(DATA, list(...)))
  y$data.name <- DNAME
  y
}
