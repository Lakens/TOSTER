#' @title Bootstrapped TOST with log t-tests
#' @description A function for a bootstrap method for TOST with all types of t-tests.
#' @inheritParams boot_t_TOST
#' @inheritParams log_TOST
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
#' @details The implemented test(s) corresponds to the proposal of Chapter 16 of Efron and Tibshirani (1993). Returns TOSTt class object with boostrapped based results. Please note that the repeated measures "corrected" effect size is not available at this time.
#' @section References:
#'
#' Efron, B., & Tibshirani, R. J. (1994). An introduction to the bootstrap. CRC press.
#'
#' He, Y., Deng, Y., You, C., & Zhou, X. H. (2022). Equivalence tests for ratio of means in bioequivalence studies under crossover design. Statistical Methods in Medical Research, 09622802221093721.
#' @importFrom stats var quantile
#' @name boot_log_TOST
#' @family compare studies
#' @export boot_log_TOST

boot_log_TOST <- function(x, ...){
  UseMethod("boot_log_TOST")
}

#' @rdname boot_log_TOST
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
                                R = 1999, ...){
  hypothesis = match.arg(hypothesis)
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
    method <- if (paired) "Bootstrapped Paired t-test" else "Bootstrapped One Sample t-test"
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
  boot.cint <- quantile(m_vec, c(alpha, 1 - alpha ))
  d.cint <- quantile(d_vec, c(alpha, 1 - alpha ))
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
    call = match.call()
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
