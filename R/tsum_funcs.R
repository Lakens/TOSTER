tsum_test = function (m1, sd1 = NULL, n1 = NULL,
                      m2 = NULL, sd2 = NULL, n2 = NULL,
                      r12 = NULL,
                      paired = FALSE,
                      alternative = "two.sided",
                      mu = 0, var.equal = FALSE,
                      conf.level = 0.95) {
  alt.expanded <- if (!missing(alternative))
    char.expand(alternative, c("two.sided", "greater", "less"),
                stop("argument 'alternative' must match one of \"greater\", \"less\", \"two.sided\"."))
  else alternative
  if (!missing(mu))
    if ((length(mu) != 1) || !is.finite(mu))
      stop("argument 'mu' must be a single finite numeric value.")
  if (!missing(conf.level))
    if ((length(conf.level) != 1) || !is.finite(conf.level) ||
        (conf.level <= 0) || (conf.level >= 1))
      stop("argument 'conf.level' must be a single number greater than zero and less than one \n.")
  if (!is.null(m1) && is.null(m2) && is.null(n1) &&
      is.null(sd1))
    stop("You must enter the value for both sd1 and n1")
  if (is.null(n1) && !is.null(m1) && !is.null(sd1) &&
      is.null(m2))
    stop("You must enter the value for n1")
  if (is.null(sd1) && !is.null(m1) && !is.null(n1) &&
      is.null(m2))
    stop("You must enter the value for sd1")
  if (is.null(n2) && !is.null(m1) && !is.null(m2) &&
      !is.null(sd2) && !is.null(sd1) && !is.null(n1))
    stop("You must enter the value for n2")
  if (is.null(n2) && is.null(n1) && !is.null(m1) &&
      !is.null(m2) && !is.null(sd2) && !is.null(sd1))
    stop("You must enter the value for both n1 and n2")
  if (is.null(sd1) && is.null(sd2) && !is.null(m1) &&
      !is.null(m2) && !is.null(n1) && !is.null(n2))
    stop("You must enter the value for both sd1 and sd2")
  if (!is.null(sd1) && is.null(sd2) && !is.null(m1) &&
      !is.null(m2) && !is.null(n1) && !is.null(n2))
    stop("You must enter the value for sd2")
  if (is.null(n2) && is.null(sd2) && !is.null(m1) &&
      !is.null(m2) && !is.null(sd1) && !is.null(n1))
    stop("You must enter the value for both sd2 and n2")
  alpha <- 1 - conf.level
  if (is.null(m2)) {

    conf.int.xbar <- m1
    conf.int.s <- sqrt(sd1^2/n1)
    ret.val <- list(statistic = (conf.int.xbar - mu)/conf.int.s,
                    parameter = n1 - 1, estimate = conf.int.xbar,
                    stderr = conf.int.s,
                    null.value = mu, alternative = alt.expanded, method = "One-sample t-Test",
                    data.name = c("Summarized x"))
    names(ret.val$estimate) <- "mean of x"
    names(ret.val$null.value) <- "mean"
  } else if(paired == TRUE){
    if(n1 != n2){
      warning("Unequal number of pairs; using smallest n")
      n1 = min(n1,n2)
    }

    conf.int.xbar = m1-m2
    sd_diff = sqrt(sd1^2+sd2^2-2*r12*sd1*sd2)
    conf.int.s <- sqrt(sd_diff^2/n1)

    ret.val <- list(statistic = (conf.int.xbar - mu)/conf.int.s,
                    parameter = n1 - 1, estimate = conf.int.xbar,
                    stderr = conf.int.s,
                    null.value = mu, alternative = alt.expanded, method = "Paired t-test",
                    data.name = c("Summarized mean differences"))
    names(ret.val$estimate) <- "mean of differences"
    names(ret.val$null.value) <- "difference in means"


  } else {
    m1 <- m1
    m2 <- m2
    conf.int.xbar <- m1 - m2
    var.x <- sd1^2
    var.y <- sd2^2
    conf.int.s <- if (var.equal)
      sqrt((((n1 - 1) * var.x + (n2 - 1) * var.y) *
              (1/n1 + 1/n2))/(n1 + n2 - 2))
    else sqrt((var.x/n1) + (var.y/n2))
    ret.val <- c(if (var.equal) list(method = "Standard Two-Sample t-Test",
                                     parameter = n1 + n2 - 2) else list(method = "Welch Modified Two-Sample t-Test",
                                                                           parameter = {
                                                                             const <- 1/(1 + (n1 * var.y)/(n2 * var.x))
                                                                             1/((const^2)/(n1 - 1) + ((1 - const)^2)/(n2 -
                                                                                                                         1))
                                                                           }), list(statistic = (conf.int.xbar - mu)/conf.int.s,
                                                                                    estimate = c(m1, m2), stderr = conf.int.s, null.value = mu, alternative = alt.expanded,
                                                                                    data.name = paste("Summarized ", deparse(substitute(x)),
                                                                                                      " and ", deparse(substitute(y)), sep = "")))
    names(ret.val$estimate) <- c("mean of x", "mean of y")
    names(ret.val$null.value) <- "difference in means"
  }
  ret.val <- c(ret.val, switch(alt.expanded, two.sided = {
    conf.int.hw <- qt((1 - alpha/2), ret.val$parameter) *
      conf.int.s
    list(p.value = 2 * pt(-abs(ret.val$statistic), ret.val$parameter),
         conf.int = c(conf.int.xbar - conf.int.hw, conf.int.xbar +
                        conf.int.hw))
  }, greater = {
    list(p.value = 1 - pt(ret.val$statistic, ret.val$parameter),
         conf.int = c(conf.int.xbar - qt((1 - alpha), ret.val$parameter) *
                        conf.int.s, NA))
  }, less = {
    list(p.value = pt(ret.val$statistic, ret.val$parameter),
         conf.int = c(NA, conf.int.xbar + qt((1 - alpha),
                                             ret.val$parameter) * conf.int.s))
  }))
  names(ret.val$statistic) <- "t"
  names(ret.val$parameter) <- "df"
  attr(ret.val$conf.int, "conf.level") <- conf.level
  ret.val <- ret.val[c("statistic", "parameter", "p.value", "stderr",
                       "conf.int", "estimate", "null.value", "alternative",
                       "method", "data.name")]
  oldClass(ret.val) <- "htest"
  return(ret.val)
}
