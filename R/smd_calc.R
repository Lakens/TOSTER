#' @title SMD Calculation
#' @description
#' `r lifecycle::badge('stable')`
#'
#' A function to only calculate standardized mean differences.
#' @inheritParams t_TOST
#' @inheritParams boot_t_TOST
#' @param mu Null value. Deviating from zero will give the x-y-mu.
#' @details For details on the calculations in this function see vignette("SMD_calcs").
#' @return A data frame containing the SMD estimates.
#' @examples
#' \dontrun{
#' smd_calc(formula = extra ~ group,data = sleep, paired = TRUE, smd_ci = "nct")
#' }
#' @name smd_calc
#' @family effect sizes
#' @export smd_calc

#smd_calc <- setClass("smd_calc")
smd_calc <- function(x, ...,
                     paired = FALSE,
                     var.equal = FALSE,
                     alpha = 0.05,
                     bias_correction = TRUE,
                     rm_correction = FALSE,
                     glass = NULL,
                     smd_ci = c("nct", "goulet", "t", "z")){
  UseMethod("smd_calc")
}

#' @rdname smd_calc
#' @importFrom stats sd cor na.omit setNames t.test terms nlm optim optimize
#' @method smd_calc default
#' @export

# @method smd_calc default
smd_calc.default = function(x,
                            y = NULL,
                            paired = FALSE,
                            var.equal = FALSE,
                            alpha = 0.05,
                            mu = 0,
                            bias_correction = TRUE,
                            rm_correction = FALSE,
                            glass = NULL,
                            smd_ci = c("nct", "goulet", "t", "z"),
                            ...) {

  if(is.null(glass)){
    glass = "no"
  }
  smd_ci = match.arg(smd_ci)

  if(bias_correction){
    smd_type = 'g'
  } else {
    smd_type = 'd'
  }

  if(is.null(y)){
    sample_type = "One Sample"
  } else if(paired == TRUE) {
    sample_type = "Paired Sample"
  } else {
    sample_type = "Two Sample"
  }

  if(glass == "glass1" || glass == "glass2"){
    if(glass == "glass1"){
      denom = "glass1"
    }

    if(glass == "glass2"){
      denom = "glass2"
    }
  } else{
    if(sample_type != "Two Sample" ){
      if(rm_correction){
        denom = "rm"
      } else {
        denom = "z"
      }
    } else{
      denom = "d"
    }
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
    m2 <- mean(i2) + mu
    sd1  <- sd(i1)
    sd2  <- sd(i2)
    r12 <- cor(i1, i2)

    # Calculate Cohens d
    cohen_res = d_est_pair(
      n = n,
      m1 = m1,
      m2 = m2,
      sd1 = sd1,
      sd2 = sd2,
      r12 = r12,
      type = smd_type,
      denom = denom,
      alpha = alpha/2,
      smd_ci = smd_ci
    )

  } else if(!missing(y)){

    x1 = na.omit(x)
    y1 = na.omit(y)
    n1 = length(x1)
    n2 = length(y1)
    m1 = mean(x1)
    m2 = mean(y1) + mu
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
      alpha = alpha/2,
      denom = denom,
      smd_ci = smd_ci
    )

  } else {

    x1 = na.omit(x)
    n1 = length(x1)
    m1 = mean(x1) + mu
    sd1 = sd(x1)

    cohen_res = d_est_one(
      n = n1,
      mu = m1,
      sd = sd1,
      type = smd_type,
      testValue = 0,
      alpha = alpha/2,
      smd_ci = smd_ci
    )

  }

  effsize = data.frame(
    estimate = c(cohen_res$d),
    SE = c(cohen_res$d_sigma),
    lower.ci = c(cohen_res$dlow),
    upper.ci = c( cohen_res$dhigh),
    conf.level = c((1-alpha)),
    row.names = c(cohen_res$smd_label)
  )


  return(effsize)

}

#' @rdname smd_calc
#' @method smd_calc formula
#' @export

smd_calc.formula = function(formula,
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
  y <- do.call("smd_calc", c(DATA, list(...)))
  #y$data.name <- DNAME
  y

}

# Bootstrap -------

#' @rdname smd_calc

#smd_calc <- setClass("smd_calc")
boot_smd_calc <- function(x, ...,
                     paired = FALSE,
                     var.equal = FALSE,
                     alpha = 0.05,
                     bias_correction = TRUE,
                     rm_correction = FALSE,
                     glass = NULL,
                     boot_ci = c("stud","basic","perc"),
                     R = 1999){
  UseMethod("boot_smd_calc")
}



# @method smd_calc default
boot_smd_calc.default = function(x,
                                 y = NULL,
                                 paired = FALSE,
                                 var.equal = FALSE,
                                 alpha = 0.05,
                                 mu = 0,
                                 bias_correction = TRUE,
                                 rm_correction = FALSE,
                                 glass = NULL,
                                 boot_ci = c("stud","basic","perc"),
                                 R = 1999,
                                 ...) {
  boot_ci = match.arg(boot_ci)
  if(paired == TRUE && !missing(y)){
    i1 <- x
    i2 <- y
    data <- data.frame(x = i1, y = i2)
    data <- na.omit(data)
    raw_smd = smd_calc(x = data$x,
                       y = data$y,
                       paired = paired,
                       var.equal = var.equal,
                       alpha = alpha,
                       mu = mu,
                       bias_correction = bias_correction,
                       rm_correction = rm_correction,
                       glass = glass,
                       smd_ci = "z")

    boots = c()
    boots_se = c()

    for(i in 1:R){
      sampler = sample(1:nrow(data), replace = TRUE)
      res_boot = smd_calc(x = data$x[sampler],
                         y = data$y[sampler],
                         paired = paired,
                         var.equal = var.equal,
                         alpha = alpha,
                         mu = mu,
                         bias_correction = bias_correction,
                         rm_correction = rm_correction,
                         glass = glass,
                         smd_ci = "z")
      boots = c(boots, res_boot$estimate)
      boots_se = c(boots_se, res_boot$SE)
    }


  } else if(!missing(y)){

    i1 <- na.omit(x)
    i2 <- na.omit(y)
    data <- data.frame(values = c(i1,i2),
                       group = c(rep("x",length(i1)),
                                      rep(i2)))

    raw_smd = smd_calc(x = i1,
                       y = i2,
                       paired = paired,
                       var.equal = var.equal,
                       alpha = alpha,
                       mu = mu,
                       bias_correction = bias_correction,
                       rm_correction = rm_correction,
                       glass = glass,
                       smd_ci = "z")

    boots = c()
    boots_se = c()
    for(i in 1:R){
      sampler = sample(1:nrow(data), replace = TRUE)
      boot_dat = data[sampler,]
      x_boot = subset(boot_dat,
                      group == "x")
      y_boot = subset(boot_dat,
                      group == "y")
      res_boot = smd_calc(x = x_boot,
                          y = y_boot,
                          paired = paired,
                          var.equal = var.equal,
                          alpha = alpha,
                          mu = mu,
                          bias_correction = bias_correction,
                          rm_correction = rm_correction,
                          glass = glass,
                          smd_ci = "z")
      boots = c(boots, res_boot$estimate)
      boots_se = c(boots_se, res_boot$SE)
    }

  } else {

    x1 = na.omit(x)
    n1 = length(x1)
    raw_smd = smd_calc(x = x1,
                       paired = paired,
                       var.equal = var.equal,
                       alpha = alpha,
                       mu = mu,
                       bias_correction = bias_correction,
                       rm_correction = rm_correction,
                       glass = glass,
                       smd_ci = "z")

    boots = c()
    boots_se = c()
    for(i in 1:R){
      sampler = sample(1:nrow(x1), replace = TRUE)
      x_boot = x1[sampler]

      res_boot = smd_calc(x = x_boot,
                          paired = paired,
                          var.equal = var.equal,
                          alpha = alpha,
                          mu = mu,
                          bias_correction = bias_correction,
                          rm_correction = rm_correction,
                          glass = glass,
                          smd_ci = "z")
      boots = c(boots, res_boot$estimate)
      boots_se = c(boots_se, res_boot$SE)
    }

  }

  ci = switch(boot_ci,
              "stud" = stud(boots, boots_se,
                            se0=raw_smd$SE, t0 = raw_smd$estimate,
                            alpha),
              "perc" = perc(boots, alpha),
              "basic" = basic(boots, t0 = raw_smd$estimate, alpha))

  effsize = data.frame(
    estimate = raw_smd$estimate,
    bias = raw_smd$estimate - median(boots),
    SE = sd(boots),
    lower.ci = ci[1],
    upper.ci = ci[2],
    conf.level = c((1-alpha)),
    boot_ci = boot_ci,
    row.names = c(raw_smd$smd_label)
  )


  return(effsize)

}

#' @rdname smd_calc
#' @method smd_calc formula
#' @export

boot_smd_calc.formula = function(formula,
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
  y <- do.call("boot_smd_calc", c(DATA, list(...)))
  #y$data.name <- DNAME
  y

}
