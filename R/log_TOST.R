#' @title TOST with log transformed t-tests
#' @description A function for TOST with all types of t-tests.
#' @param eqb Equivalence bound; default is 1.25 (FDA guidelines). Can provide 1 value (reciprical value is taken as the lower bound) or 2 specific values that represent the upper and lower equivalence bounds.
#' @param null Null hypothesis value for two-tailed test (default is 1).
#' @inheritParams t_TOST
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
#' }
#' @name log_TOST
#' @export log_TOST

#log_TOST <- setClass("log_TOST")
log_TOST <- function(x, ...,
                     hypothesis = "EQU",
                     paired = FALSE,
                     var.equal = FALSE,
                     eqb = 1.25,
                     alpha = 0.05,
                     null = 1){
  UseMethod("log_TOST")
}

#' @rdname log_TOST
#' @importFrom stats sd cor na.omit setNames t.test terms nlm optim optimize
#' @method log_TOST default
#' @export

# @method log_TOST default
log_TOST.default = function(x,
                            y = NULL,
                            hypothesis = c("EQU","MET"),
                            eqb = 1.25,
                            alpha = 0.05,
                            null = 1,
                            ...) {
  hypothesis = match.arg(hypothesis)

  if(is.null(y)){
    stop("One sample tests not supported at this time.")
  } else if(paired == TRUE) {
    sample_type = "Paired Sample"
  } else {
    sample_type = "Two Sample"
  }


  if(hypothesis == "EQU"){
    alt_low = "greater"
    alt_high = "less"
    test_hypothesis = "Hypothesis Tested: Equivalence"

  } else if(hypothesis == "MET"){
    alt_low = "less"
    alt_high = "greater"
    test_hypothesis = "Hypothesis Tested: Minimal Effect"

  } else{
    stop("hypothesis must be set to EQU or MET")
  }

  if(any(eqb <= 0)){
    stop("Equivalence bounds must be a positive value")
  }

  if(!is.numeric(alpha) || alpha <=0 || alpha >=1){
    stop("The alpha must be a numeric value between 0 and 1")
  }

  if (!is.null(y)) {
    dname <- paste(deparse(substitute(x)), "and",
                   deparse(substitute(y)))
  }
  else {
    dname <- deparse(substitute(x))
  }

  if(any(x < 0) || any(y < 0)){
    stop("Negative values detected. Values must be on ratio scale (true zero).")
  }
  x = log(x)
  y = log(y)

  tresult = t.test(x = x,
                   y = y,
                   paired = paired,
                   var.equal = var.equal,
                   mu = log(null),
                   conf.level = 1 - alpha*2,
                   alternative = "two.sided")

  if(paired == TRUE && !missing(y)){
    i1 <- x
    i2 <- y
    data <- data.frame(i1 = i1, i2 = i2)
    data <- na.omit(data)
    colnames(data) = c("i1", "i2")
    data2 =  data
    data2$diff = data2$i2 - data2$i1

    n <- nrow(data)
    i1 <- data$i1
    i2 <- data$i2
    m1 <- mean(i1)
    m2 <- mean(i2)
    sd1  <- sd(i1)
    sd2  <- sd(i2)
    r12 <- cor(i1, i2)

    # Calculate log rom ------
    cohen_res = d_est_pair(
      n = n,
      m1 = m1,
      m2 = m2,
      sd1 = sd1,
      sd2 = sd2,
      r12 = r12,
      type = smd_type,
      denom = denom,
      alpha = alpha,
      smd_ci = smd_ci
    )

  } else {

    x1 = na.omit(x)
    y1 = na.omit(y)

    n1 = length(x1)
    n2 = length(y1)

    m1 = mean(x1)
    m2 = mean(y1)

    sd1 = sd(x1)
    sd2 = sd(y1)

    cohen_res = d_est_ind(
      n1 = n1,
      n2 = n2,
      m1 = m1,
      m2 = m2,
      sd1 = sd1,
      sd2 = sd2,
      type = smd_type,
      var.equal = var.equal,
      alpha = alpha,
      denom = denom,
      smd_ci = smd_ci
    )

  }

  if(!missing(eqb)){
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
  }

  if(hypothesis == "EQU"){
    null_hyp = paste0(round(low_eqbound,2),
                      " >= (Mean1/Mean2) or (Mean1/Mean2) >= ",
                      round(high_eqbound,2))
    alt_hyp = paste0(round(low_eqbound,2),
                     " < (Mean1/Mean2) < ",
                     round(high_eqbound,2))
  } else if(hypothesis == "MET"){
    null_hyp = paste0(round(low_eqbound,2),
                      " <= (Mean1/Mean2)  <= ",
                      round(high_eqbound,2))
    alt_hyp = paste0(round(low_eqbound,2),
                     " > (Mean1/Mean2) or (Mean1/Mean2)  > ",
                     round(high_eqbound,2))
  }

  low_ttest <- t.test(
    y = y,
    x = x,
    paired = paired,
    var.equal = var.equal,
    alternative = alt_low,
    mu = log(low_eqbound),
    conf.level = 1-alpha*2
  )

  high_ttest <- t.test(
    y = y,
    x = x,
    paired = paired,
    var.equal = var.equal,
    alternative = alt_high,
    mu = log(high_eqbound),
    conf.level = 1-alpha*2
  )

  if(hypothesis == "EQU"){
    pTOST = max(low_ttest$p.value,
                high_ttest$p.value) # get highest p value for TOST result
    tTOST = ifelse(abs(low_ttest$statistic) < abs(high_ttest$statistic),
                   low_ttest$statistic,
                   high_ttest$statistic) #Get lowest t-value for summary TOST result
  } else {
    pTOST = min(low_ttest$p.value,
                high_ttest$p.value) # get highest p value for TOST result
    tTOST = ifelse(abs(low_ttest$statistic) > abs(high_ttest$statistic),
                   low_ttest$statistic,
                   high_ttest$statistic) #Get lowest t-value for summary TOST result
  }


  TOST = data.frame(
    t = c(tresult$statistic,
          low_ttest$statistic,
          high_ttest$statistic),
    SE = c(tresult$stderr,
           low_ttest$stderr,
           high_ttest$stderr),
    df = c(tresult$parameter,
           low_ttest$parameter,
           high_ttest$parameter),
    p.value = c(tresult$p.value,
                low_ttest$p.value,
                high_ttest$p.value),
    row.names = c("t-test","TOST Lower","TOST Upper")
  )

  eqb = data.frame(
    type = c("log(Means Ratio)","Means Ratio"),
    low_eq = c(log(low_eqbound),low_eqbound_d),
    high_eq = c(log(high_eqbound),high_eqbound_d)
  )

  effsize = data.frame(
    estimate = c(tresult$statistic * tresult$stderr,
                 exp(tresult$statistic * tresult$stderr)),
    SE = c(tresult$stderr,NA),
    lower.ci = c(tresult$conf.int[1], exp(tresult$conf.int[1])),
    upper.ci = c(tresult$conf.int[2], exp(tresult$conf.int[2]),
    conf.level = c((1-alpha*2),(1-alpha*2)),
    row.names = c("log(Means Ratio)","Means Ratio")
  )
  TOSToutcome<-ifelse(pTOST<alpha,"significant","non-significant")
  testoutcome<-ifelse(tresult$p.value<alpha,"significant","non-significant")

  # Change text based on two tailed t test if mu is not zero
  if(null == 1){
    mu_text = "one"
  } else {
    mu_text = null
  }

  if(hypothesis == "EQU"){
    #format(low_eqbound, digits = 3, nsmall = 3, scientific = FALSE)
    TOST_restext = paste0("The equivalence test was ",TOSToutcome,", t(",round(tresult$parameter, digits=2),") = ",format(tTOST, digits = 3, nsmall = 3, scientific = FALSE),", p = ",format(pTOST, digits = 3, nsmall = 3, scientific = TRUE),sep="")
  } else {
    TOST_restext = paste0("The minimal effect test was ",TOSToutcome,", t(",round(tresult$parameter, digits=2),") = ",format(tTOST, digits = 3, nsmall = 3, scientific = FALSE),", p = ",format(pTOST, digits = 3, nsmall = 3, scientific = TRUE),sep="")
  }

  ttest_restext = paste0("The null hypothesis test was ",testoutcome,", t(",round(tresult$parameter, digits=2),") = ",format(tresult$statistic, digits = 3, nsmall = 3, scientific = FALSE),", p = ",format(tresult$p.value, digits = 3, nsmall = 3, scientific = TRUE),sep="")
  if (hypothesis == "EQU"){
    if(tresult$p.value <= alpha && pTOST <= alpha){
      combined_outcome <- paste0("NHST: reject null significance hypothesis that the effect is equal to ",mu_text," \n",
                                 "TOST: reject null equivalence hypothesis")
    }
    if(tresult$p.value < alpha && pTOST > alpha){
      combined_outcome <- paste0("NHST: reject null significance hypothesis that the effect is equal to ",mu_text," \n ",
                                 "TOST: don't reject null equivalence hypothesis")
      # paste0("statistically different from ",mu_text," and not statistically equivalent")
    }
    if(tresult$p.value > alpha && pTOST <= alpha){
      combined_outcome <- paste0("NHST: don't reject null significance hypothesis that the effect is equal to ",mu_text," \n ",
                                 "TOST: reject null equivalence hypothesis")
      #paste0("statistically not different from ",mu_text," and statistically equivalent")
    }
    if(tresult$p.value > alpha && pTOST > alpha){
      combined_outcome <- paste0("NHST: don't reject null significance hypothesis that the effect is equal to ",mu_text," \n ",
                                 "TOST: don't reject null equivalence hypothesis")
      #paste0("statistically not different from ",mu_text," and not statistically equivalent")
    }
  } else {
    if(tresult$p.value <= alpha && pTOST <= alpha){
      combined_outcome <- paste0("NHST: reject null significance hypothesis that the effect is equal to ",mu_text," \n ",
                                 "TOST: reject null MET hypothesis")
      #paste0("statistically different from ",mu_text," and statistically greater than the minimal effect threshold")
    }
    if(tresult$p.value < alpha && pTOST > alpha){
      combined_outcome <- paste0("NHST: reject null significance hypothesis that the effect is equal to ",mu_text," \n ",
                                 "TOST: don't reject null MET hypothesis")
      #paste0("statistically different from ",mu_text," but not statistically greater than the minimal effect threshold")
    }
    if(tresult$p.value > alpha && pTOST <= alpha){
      combined_outcome <- paste0("NHST: don't reject null significance hypothesis that the effect is equal to ",mu_text," \n ",
                                 "TOST: reject null MET hypothesis")
      #paste0("statistically not different from ",mu_text," and statistically greater than the minimal effect threshold")
    }
    if(tresult$p.value > alpha && pTOST > alpha){
      combined_outcome <- paste0("NHST: don't reject null significance hypothesis that the effect is equal to ",mu_text," \n ",
                                 "TOST: don't reject null MET hypothesis")
      #paste0("statistically not different from ",mu_text," and not statistically greater than the minimal effect threshold")
    }
  }

  decision = list(
    TOST = TOST_restext,
    ttest = ttest_restext,
    combined = combined_outcome
  )

  #message(cat("Based on the equivalence test and the null-hypothesis test combined, we can conclude that the observed effect is ",combined_outcome,".",sep=""))


  rval = list(
    TOST = TOST,
    eqb = eqb,
    alpha = alpha,
    method = tresult$method,
    hypothesis = test_hypothesis,
    effsize = effsize,
    smd = cohen_res,
    decision = decision,
    data.name = dname,
    call = match.call()
  )

  class(rval) = "TOSTt"

  return(rval)

}

#' @rdname log_TOST
#' @method log_TOST formula
#' @export

log_TOST.formula = function(formula,
                            data,
                            subset,
                            na.action, ...) {

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
  y <- do.call("log_TOST", c(DATA, list(...)))
  y$data.name <- DNAME
  y

}

#' @importFrom stats sd
#' @keywords internal

logrom_calc = function(paired = FALSE,
                       bias_c = TRUE,
                       vtype = "LS",
                       m1i,
                       sd1i,
                       n1i,
                       m2i,
                       sd2i,
                       n2i,
                       ri = NULL) {
  if (!paired) {
    yi <- log(m1i / m2i)


    ### sample size weighted average of the coefficient of variation in group 1
    mn1wcvi <- .wmean(sd1i / m1i,
                      n1i,
                      na.rm = TRUE)
    ### sample size weighted average of the coefficient of variation in group 2
    mn2wcvi <- .wmean(sd2i / m2i,
                      n2i,
                      na.rm = TRUE)
    ### sample size weighted average of the two CV values

    mnwcvi  <-
      (sum(n1i * (sd1i / m1i)) + sum(n2i * (sd2i / m2i))) / sum((n1i +
                                                                   n2i))

    ### large sample approximation to the sampling variance (does not assume homoscedasticity)
    if (vtype == "LS") {
      vi <-
        sd1i ^ 2 / (n1i * m1i ^ 2) + sd2i ^ 2 / (n2i * m2i ^
                                                   2)
    }
    ### estimator assuming homoscedasticity
    if (vtype == "HO") {
      mi   <- n1i+n2i - 2
      sdpi <- sqrt(((n1i-1)*sd1i^2 + (n2i-1)*sd2i^2)/mi)
      vi <-
        sdpi ^ 2 / (n1i * m1i ^ 2) + sdpi ^ 2 / (n2i * m2i ^
                                                   2)
    }
    ### estimator using the weighted averages of the CV values
    if (vtype == "AV") {
      vi <- mn1wcvi ^ 2 / n1i + mn2wcvi ^ 2 / n2i
    }
    ### estimator using the weighted average of two weighted averages of the CV values
    if (vtype == "AVHO"){
      vi <- mnwcvi ^ 2 * (1 / n1i + 1 / n2i)
    }

  }

  if (paired) {
    yi <- log(m1i / m2i)
    vi <-
      sd1i ^ 2 / (n1i * m1i ^ 2) + sd2i ^ 2 / (n1i * m2i ^ 2) - 2 * ri * sd1i *
      sd2i / (m1i * m2i * n1i)

  }

  if(bias_c){
    J = 0.5 * (sd1i^2 / (n1i * m1i^2) - sd2i^2 / (n2i * m2i^2))
    yi = yi + J

    Jvar = 0.5 * (sd1i^4 / (n1i^2 * m1i^4) - sd2i^4 / (n2i^2 * m2i^4))
    vi = vi + Jvar
  }


  rval = list(
    log_rom = yi,
    var_rom = vi
  )
  return(rval)
}
