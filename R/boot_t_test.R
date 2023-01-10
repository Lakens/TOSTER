#' @title Bootstrapped t-test
#' @description A function for a bootstrap method for Tt-tests.
#' @param R number of bootstrap replicates
#' @inheritParams simple_htest
#' @details For details on the calculations in this function see vignette("robustTOST").
#' @return A list with class \code{"htest"} containing the following components:
#' \describe{
#'   \item{\code{statistic}}{the value of the t-statistic.}
#'   \item{\code{parameter}}{the degrees of freedom for the t-statistic.}
#'   \item{\code{p.value}}{the p-value for the test.}
#'   \item{\code{conf.int}}{a confidence interval for the mean appropriate to the specified alternative hypothesis.}
#'   \item{\code{estimate}}{the estimated mean or difference in means depending on whether it was a one-sample test or a two-sample test.}
#'   \item{\code{null.value}}{the specified hypothesized value of the mean or mean difference. May be 2 values.}
#'   \item{\code{stderr}}{the standard error of the mean (difference), used as denominator in the t-statistic formula.}
#'   \item{\code{alternative}}{a character string describing the alternative hypothesis.}
#'   \item{\code{method}}{a character string indicating what type of t-test was performed.}
#'   \item{\code{data.name}}{a character string giving the name(s) of the data..}
#' }
#' @details The implemented test(s) corresponds to the proposal of Chapter 16 of Efron and Tibshirani (1994).
#'
#' For details on the calculations in this function see vignette("robustTOST").
#' @section References:
#'
#' Efron, B., & Tibshirani, R. J. (1994). An introduction to the bootstrap. CRC press.
#' @family Robust TOST
#' @importFrom stats var quantile
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
                                paired = FALSE,
                                alternative = c("two.sided",
                                                "less",
                                                "greater",
                                                "equivalence",
                                                "minimal.effect"),
                                mu = 0,
                                alpha = 0.05,
                                R = 1999, ...){
  alternative = match.arg(alternative)

  if(!missing(mu) && (length(mu) != 1 || is.na(mu))) {
    stop("'mu' must be a single number")
  }


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
                           paired = paired,
                           alternative = alternative,
                           mu = 0,
                           alpha = 0.05)
  m_vec <- rep(NA, times=length(R)) # mean difference vector
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

      m_vec[i] <- mean(dat, na.rm=TRUE) # mean difference vector

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
  tstat = null_test
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

#' @rdname boot_t_test
#' @method boot_t_test formula
#' @export
#'
boot_t_test.formula <- function (formula, data, subset, na.action, ...){
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
  y <- do.call("boot_t_test", c(DATA, list(...)))
  y$data.name <- DNAME
  y
}
