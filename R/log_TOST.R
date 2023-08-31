#' @title TOST with log transformed t-tests
#' @description
#' `r lifecycle::badge('stable')`
#'
#'A function for TOST on the log-transformed data using parametric t-tests.
#' @param eqb Equivalence bound; default is 1.25 (FDA guidelines). Can provide 1 value (reciprocal value is taken as the lower bound) or 2 specific values that represent the upper and lower equivalence bounds.
#' @param null Null hypothesis value for a two-tailed test (default is 1).
#' @inheritParams t_TOST
#' @details
#' For details on the calculations in this function see `vignette("robustTOST")`.
#'
#' For two-sample tests, the test is of \eqn{\bar log(x) - \bar log(y)} (mean of x minus mean of y).
#' For paired samples, the test is of the difference scores (z),
#' wherein \eqn{z =  log(x) - log(y) = log(x)/log(y)}, and the test is of \eqn{\bar z} (mean of the difference/ratio scores).
#'
#' @return An S3 object of class
#'   `"TOSTt`" is returned containing the following slots:
#'
#'   - "TOST": A table of class \code{"data.frame"} containing two-tailed t-test and both one-tailed results.
#'   - "eqb": A table of class \code{"data.frame"} containing equivalence bound settings.
#'   - "effsize":  table of class \code{"data.frame"} containing effect size estimates.
#'   - "hypothesis": String stating the hypothesis being tested
#'   - "smd": List containing the results of the means ratio calculation.
#'      - Items include: d (means ratio estimate), dlow (lower CI bound), dhigh (upper CI bound), d_df (degrees of freedom for SMD), d_sigma (SE), d_lambda (non-centrality), J (bias correction), smd_label (type of SMD), d_denom (denominator calculation)
#'   - "alpha": Alpha level set for the analysis.
#'   - "method": Type of t-test.
#'   - "decision": List included text regarding the decisions for statistical inference.
#'
#' @references
#' He, Y., Deng, Y., You, C., & Zhou, X. H. (2022). Equivalence tests for ratio of means in bioequivalence studies under crossover design. Statistical Methods in Medical Research, 09622802221093721.
#'
#' Food and Drug Administration (2014). Bioavailability and Bioequivalence Studies Submitted in NDAs or INDs — General Considerations.
#' Center for Drug Evaluation and Research. Docket: FDA-2014-D-0204.
#' https://www.fda.gov/regulatory-information/search-fda-guidance-documents/bioavailability-and-bioequivalence-studies-submitted-ndas-or-inds-general-considerations
#' @examples
#' data(mtcars)
#' # Default FDA bioequivalence bounds
#' log_TOST(mpg ~ am,
#' data = mtcars)
#' @family Robust tests
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
                            var.equal = FALSE,
                            paired = FALSE,
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
  logrom = (tresult$statistic) * tresult$stderr + log(null)
  logSE = tresult$stderr

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
    d_df = n-1

    # Calculate log rom ------
    #log_romres = logrom_calc(
    #  paired = TRUE,
    #  bias_c = TRUE,
    #  vtype = "LS",
    #  m1i = m1,
    #  sd1i = sd1,
    #  n1i = n1,
    #  m2i = m2,
    #  sd2i = sd2,
    #  n2i = n2,
    #  ri = r12
    #)

    rom_res = list(
      d = exp(logrom),
      d_df = d_df,
      dlow = exp(tresult$conf.int[1]),
      dhigh = exp(tresult$conf.int[2]),
      d_sigma = logSE,
      d_lambda = NULL,
      #hn = hn,
      smd_label = "Means Ratio",
      J = 0.5 * (sd1^2 / (n * m1^2) - sd2^2 / (n * m2^2)),
      d_denom = 1,
      ntilde = 1,
      t_stat = NULL,
      smd_ci = "t"
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
    d_df = n1 + n2 - 2

    rom_res = list(
      d = exp(logrom),
      d_df = d_df,
      dlow = exp(tresult$conf.int[1]),
      dhigh = exp(tresult$conf.int[2]),
      d_sigma = logSE,
      d_lambda = NULL,
      #hn = hn,
      smd_label = "Means Ratio",
      J = 0.5 * (sd1^2 / (n1 * m1^2) - sd2^2 / (n2 * m2^2)),
      d_denom = 1,
      ntilde = 1,
      t_stat = NULL,
      smd_ci = "t"
    )

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
    low_eq = c(log(low_eqbound),low_eqbound),
    high_eq = c(log(high_eqbound),high_eqbound)
  )

  effsize = data.frame(
    estimate = c(logrom,
                 exp(logrom)),
    SE = c(logSE,NA),
    lower.ci = c(tresult$conf.int[1], exp(tresult$conf.int[1])),
    upper.ci = c(tresult$conf.int[2], exp(tresult$conf.int[2])),
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

  combined_outcome = tost_decision(hypothesis = hypothesis,
                                    alpha = alpha,
                                    pvalue = tresult$p.value,
                                    pTOST = pTOST,
                                    mu_text = mu_text)


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
    method = paste0("Log-transformed ",tresult$method),
    hypothesis = test_hypothesis,
    effsize = effsize,
    smd = rom_res,
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


