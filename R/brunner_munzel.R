#' @title Brunner-Munzel Test
#' @description
#' `r lifecycle::badge("maturing")`
#'
#' This is a generic function that performs a generalized asymptotic Brunner-Munzel test in a fashion similar to [t.test].
#' @param x a (non-empty) numeric vector of data values.
#' @param y an optional (non-empty) numeric vector of data values.
#' @param formula a formula of the form lhs ~ rhs where lhs is a numeric variable giving the data values and rhs either 1 for a one-sample or paired test or a factor with two levels giving the corresponding groups. If lhs is of class "Pair" and rhs is 1, a paired test is done.
#' @param data an optional matrix or data frame (or similar: see model.frame) containing the variables in the formula formula. By default the variables are taken from environment(formula).
#' @param paired a logical indicating whether you want a paired test.
#' @param mu 	a number specifying an optional parameter used to form the null hypothesis (Default = 0.5). This can be thought of as the null in terms of the relative effect, p = P (X < Y ) + 0.5 * P (X = Y); See ‘Details’.
#' @param perm a logical indicating whether or not to perform a permutation test over approximate t-distribution based test (default is FALSE). Highly recommend to set perm = TRUE when sample size per condition is less than 15.
#' @param max_n_perm the maximum number of permutations (default is 10000).
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less". You can specify just the initial letter.
#' @param subset an optional vector specifying a subset of observations to be used.
#' @param na.action a function which indicates what should happen when the data contain NAs. Defaults to getOption("na.action").
#' @inheritParams t_TOST
#' @param ...  further arguments to be passed to or from methods.
#' @details
#'
#' This function is made to provide a test of stochastic equality between two samples (paired or independent), and is referred to as the Brunner-Munzel test.
#'
#' This tests the hypothesis that the relative effect, discussed below, is equal to the null value (default is mu = 0.5).
#'
#' The estimate of the relative effect, which can be considered as value similar to the probability of superiority, refers to the following:
#'
#'  \deqn{\hat p = p(X<Y) + \frac{1}{2} \cdot P(X=Y)}
#'
#'  Note, for paired samples, this does *not* refer to the probability of an increase/decrease in paired sample but rather the probability that a randomly sampled value of X.
#'  This is also referred to as the "relative" effect in the literature. Therefore, the results will differ from the concordance probability provided by the ses_calc function.
#'
#'  The brunner_munzel function is based on the npar.t.test and npar.t.test.paired functions within the nparcomp package (Konietschke et al. 2015).
#'
#' @return A list with class `"htest"` containing the following components:
#'
#'   - "statistic": the value of the test statistic.
#'   - "parameter": the degrees of freedom for the test statistic.
#'   - "p.value": the p-value for the test.
#'   - "conf.int": a confidence interval for the relative effect appropriate to the specified alternative hypothesis.
#'   - "estimate": the estimated relative effect.
#'   - "null.value": the specified hypothesized value of the relative effect.
#'   - "stderr": the standard error of the relative effect.
#'   - "alternative": a character string describing the alternative hypothesis.
#'   - "method": a character string indicating what type of test was performed.
#'   - "data.name": a character string giving the name(s) of the data.
#'
#' @examples
#' data(mtcars)
#' brunner_munzel(mpg ~ am, data = mtcars)
#' @references
#' Brunner, E., Munzel, U. (2000). The Nonparametric Behrens-Fisher Problem: Asymptotic Theory and a Small Sample Approximation. Biometrical Journal 42, 17 -25.
#'
#' Neubert, K., Brunner, E., (2006). A Studentized Permutation Test for the Nonparametric Behrens-Fisher Problem. Computational Statistics and Data Analysis.
#'
#' Munzel, U., Brunner, E. (2002). An Exact Paired Rank Test. Biometrical Journal 44, 584-593.
#'
#' Konietschke, F., Placzek, M., Schaarschmidt, F., & Hothorn, L. A. (2015). nparcomp: an R software package for nonparametric multiple comparisons and simultaneous confidence intervals. Journal of Statistical Software 64 (2015), Nr. 9, 64(9), 1-17. http://www.jstatsoft.org/v64/i09/
#' @name brunner_munzel
#' @family Robust tests
#' @export brunner_munzel

#brunner_munzel <- setClass("brunner_munzel")
brunner_munzel <- function(x,
                           ...,
                           paired = FALSE,
                           alternative = c("two.sided",
                                           "less",
                                           "greater"),
                           mu = 0.5,
                           alpha = 0.05,
                           perm = FALSE,
                           max_n_perm = 10000) {

  UseMethod("brunner_munzel")
}

#' @rdname brunner_munzel
#' @importFrom stats sd cor na.omit setNames t.test terms nlm optim optimize
#' @method brunner_munzel default
#' @export

# @method brunner_munzel default
brunner_munzel.default = function(x,
                                  y,
                                  paired = FALSE,
                                  alternative = c("two.sided",
                                                  "less",
                                                  "greater"),
                                  mu = 0.5,
                                  alpha = 0.05,
                                  perm = FALSE,
                                  max_n_perm = 10000,
                                  ...) {
  alternative = match.arg(alternative)
  if(!missing(mu) &&
     ((length(mu) > 1L) || !is.finite(mu)) ||
    ( mu < 0 || mu > 1))
    stop("'mu' must be a single number between 0 and 1.")

  if(alternative == "two.sided"){
    if(mu == 0 | mu == 1){
      stop("'mu' cannot be 0 or 1 when alternative is 'two.sided'")
    }
  }
  if(!is.numeric(x)) stop("'x' must be numeric")
  if(!missing(y)) {
    if(!is.numeric(y)) stop("'y' must be numeric")
    DNAME <- paste(deparse(substitute(x)), "and",
                   deparse(substitute(y)))
    if(paired) {
      if(length(x) != length(y))
        stop("'x' and 'y' must have the same length")
      ok <- complete.cases(x,y)
      x = x[ok]
      y = y[ok]
    }
    else {
      x <- x[is.finite(x)]
      y <- y[is.finite(y)]
    }

    if(min(length(x),length(y)) < 15 && !(perm)){
      message("Sample size in at least one group is small. Permutation test (perm = TRUE) is highly recommended.")
    }
    if(min(length(x),length(y)) > 250 &&
       sum(length(x),length(y)) > 500 &&
       (perm)){
      message("Sample size is fairly large. Use of a permutation test is probably unnessary.")
    }
  } else {

    stop("'y' is missing. One sample tests currently not supported.")

  }


  # Paired -----
  if(paired){

    n = length(x)
    n1 = n+1
    df.sw = n - 1
    all_data <- c(y, x)
    N = length(all_data)

    xinverse <- c(x, y)
    x1 <- y
    x2 <- x
    rx <- rank(all_data)
    rxinverse <- rank(xinverse)
    rx1 <- rx[1:n]
    rx2 <- rx[(n+1):N]
    rix1 <- rank(x1)
    rix2 <- rank(x2)
    BM1 <- 1 / n * (rx1 - rix1)
    BM2 <- 1 / n * (rx2 - rix2)
    BM3 <- BM1 - BM2
    BM4 <- 1 / (2 * n) * (rx1 - rx2)
    pd <- mean(BM2)

    m <- mean(BM3)
    v <- (sum(BM3 ^ 2) - n * m ^ 2) / (n - 1)
    v0 <- (v == 0)
    v[v0] <- 1 / n
    test_stat <- sqrt(n) * (pd - mu) / sqrt(v)
    std_err = sqrt(v)

    if(perm == TRUE){
      METHOD = "Paired Brunner-Munzel permutation test"
      # Directly from nparcomp
      if(n<=13){
        max_n_perm=2^n
        p<-0
        for (i in 1:n){
          a<-rep(c(rep(c(i,i+n),max_n_perm/(2^i)),rep(c(i+n,i),max_n_perm/(2^i))),2^(i-1))
          p<-rbind(p,a)
        }
        p<-p[2:(n+1),]
        P<-matrix(p,ncol=max_n_perm)

        xperm<-matrix(all_data[P],nrow=N,ncol=max_n_perm)
        rxperm<-matrix(rx[P],nrow=N,ncol=max_n_perm)
      }
      else{
        P<-matrix(nrow=n,ncol=10000)
        permu<-function(all_data){
          n<-length(all_data)
          result<-sample(c(0,1),size=n,replace=TRUE)
          return(result)
        }
        P1<-apply(P,2,permu)
        P2<-rbind(P1,P1)
        xperm<-all_data*P2+xinverse*(1-P2)
        rxperm<-rx*P2+rxinverse*(1-P2)
      }
      xperm1<-xperm[1:n,]
      xperm2<-xperm[n1:N,]
      rperm1<-rxperm[1:n,]
      rperm2<-rxperm[n1:N,]
      riperm1<-apply(xperm1,2,rank)
      riperm2<-apply(xperm2,2,rank)
      BMperm2<-1/n*(rperm2-riperm2)
      BMperm3<-1/n*(rperm1-riperm1)-BMperm2
      pdperm<-colMeans(BMperm2)
      mperm3<-colMeans(BMperm3)
      vperm3<-(colSums(BMperm3^2)-n*mperm3^2)/(n-1)
      vperm30<-(vperm3==0)
      vperm3[vperm30]<-1/n
      #print(vperm3)
      Tperm<-sqrt(n)*(pdperm-mu)/sqrt(vperm3)
      p1perm<-mean(Tperm<=test_stat)
      pq1<-sort(Tperm)[(floor((1-alpha/2)*max_n_perm)+1)]
      pq2<-sort(Tperm)[(floor((1-alpha)*max_n_perm)+1)]

      p.value = switch(alternative,
                       "two.sided" = min(2*p1perm,2*(1-p1perm)),
                       "less" = p1perm,
                       "greater" =  1-p1perm)

      pd.lower = switch(alternative,
                        "two.sided" = pd-pq1*sqrt(v/n),
                        "less" = 0,
                        "greater" =  pd-pq2*sqrt(v/n))

      pd.upper = switch(alternative,
                        "two.sided" = pd+pq1*sqrt(v/n),
                        "less" = pd+pq2*sqrt(v/n),
                        "greater" =  1)
    } else {

      METHOD = "exact paired Brunner-Munzel test"
      p.value = switch(alternative,
                       "two.sided" = (2*min(pt(test_stat,n-1),
                                            1-pt(test_stat,n-1))),
                       "less" = pt(test_stat,
                                   df=df.sw),
                       "greater" =  1-pt(test_stat,
                                         df=df.sw))

      pd.lower = switch(alternative,
                        "two.sided" = pd-qt(1-alpha/2,df.sw)*sqrt(v/n),
                        "less" = 0,
                        "greater" =  pd-qt(1-alpha,df.sw)*sqrt(v/n))

      pd.upper = switch(alternative,
                        "two.sided" = pd+qt(1-alpha/2,df.sw)*sqrt(v/n),
                        "less" = pd+qt(1-alpha,df.sw)*sqrt(v/n),
                        "greater" =  1)
    }



  } else{

  # Two-sample ------
    rxy <- rank(c(x, y))
    rx <- rank(x)
    ry <- rank(y)
    n.x <- as.double(length(x))
    n.y <- as.double(length(y))
    N = n.x + n.y

    pl2 <- 1/n.y*(rxy[1:n.x]-rx)
    pl1 <- 1/n.x*(rxy[(n.x+1):N]-ry)
    pd <- mean(pl2)
    pd1 <- (pd == 1)
    pd0 <- (pd == 0)
    pd[pd1] <- 0.9999
    pd[pd0] <- 0.0001
    s1 <- var(pl2)/n.x
    s2 <- var(pl1)/n.y

    V <- N*(s1+s2)
    singular.bf <- (V == 0)
    V[singular.bf] <- N/(2 * n.x * n.y)
    std_err = sqrt(V)


    test_stat <- sqrt(N)*(pd - mu)/sqrt(V)
    df.sw <- (s1 + s2)^2/(s1^2/(n.x - 1) + s2^2/(n.y - 1))
    df.sw[is.nan(df.sw)] <- 1000

    if(perm){

      ## permutation -----
      METHOD = "two-sample Brunner-Munzel permutation test"

      Tprob<-qnorm(pd)*exp(-0.5*qnorm(pd)^2)*sqrt(N/(V*2*pi))
      P<-apply(matrix(rep(1:N,max_n_perm),ncol=max_n_perm),2,sample)
      Px<-matrix(c(x,y)[P],ncol=max_n_perm)
      Tperm<-t(apply(perm_loop(x=Px[1:n.x,],y=Px[(n.x+1):N,],
                               n.x=n.x,n.y=n.y,max_n_perm),1,sort))
      p.PERM1<-mean((test_stat <= Tperm[1,]))
      if(alternative == "two.sided"){
        c1<-0.5*(Tperm[1,floor((1-alpha/2)*max_n_perm)]+Tperm[1,ceiling((1-alpha/2)*max_n_perm)])
        c2<-0.5*(Tperm[1,floor(alpha/2*max_n_perm)]+Tperm[1,ceiling(alpha/2*max_n_perm)])
      } else {
        c1<-0.5*(Tperm[1,
                       floor((1-alpha)*max_n_perm)]+Tperm[1,
                                                          ceiling((1-alpha)*max_n_perm)])
        c2<-0.5*(Tperm[1,
                       floor(alpha*max_n_perm)]+Tperm[1,
                                                      ceiling(alpha*max_n_perm)])
      }

      lower_ci = pd-sqrt(V/N)*c1
      upper_ci = pd-sqrt(V/N)*c2

      p.value = switch(alternative,
                       "two.sided" = min(2-2*p.PERM1,
                                         2*p.PERM1),
                       "less" =  p.PERM1,
                       "greater" =  1-p.PERM1)

      pd.lower = switch(alternative,
                        "two.sided" = lower_ci,
                        "less" = 0,
                        "greater" =  lower_ci)
      pd.lower = ifelse(pd.lower < 0, 0, pd.lower)

      pd.upper = switch(alternative,
                        "two.sided" = upper_ci,
                        "less" = upper_ci,
                        "greater" =  1)
      pd.upper = ifelse(pd.upper > 1, 1, pd.upper)

      } else{

        ## t-approx ----
    METHOD = "two-sample Brunner-Munzel test"
    p.value = switch(alternative,
                     "two.sided" = min(c(2 - 2 * pt(test_stat,
                                                  df=df.sw),
                                       2 * pt(test_stat,
                                              df=df.sw))),
                     "less" = pt(test_stat, df=df.sw),
                     "greater" =  1-pt(test_stat, df=df.sw))

    pd.lower = switch(alternative,
                      "two.sided" = pd - qt(1-alpha/2, df=df.sw)/sqrt(N)*sqrt(V),
                      "less" = 0,
                      "greater" =  pd - qt(1-alpha, df=df.sw)/sqrt(N)*sqrt(V))
    pd.lower = ifelse(pd.lower < 0, 0, pd.lower)

    pd.upper = switch(alternative,
                      "two.sided" = pd + qt(1-alpha/2, df=df.sw)/sqrt(N)*sqrt(V),
                      "less" = pd + qt(1-alpha, df=df.sw)/sqrt(N)*sqrt(V),
                      "greater" =  1)
    pd.upper = ifelse(pd.upper > 1, 1, pd.upper)
    }
  }

  names(mu) <- "relative effect"
  names(test_stat) = "t"
  names(df.sw) = "df"
  cint = c(pd.lower,pd.upper)
  attr(cint,"conf.level") = ifelse(alternative == "two.sided",
                                   1-alpha,
                                   1-2*alpha)
  estimate = pd
  names(estimate) = "p(X>Y) + .5*P(X=Y)"

  rval <- list(statistic = test_stat,
               parameter = df.sw,
               p.value = as.numeric(p.value),
               estimate = estimate,
               stderr = std_err,
               conf.int = cint,
               null.value = mu,
               alternative = alternative,
               method = METHOD,
               data.name = DNAME)
  class(rval) <- "htest"
  return(rval)

}

#' @rdname brunner_munzel
#' @method brunner_munzel formula
#' @export

brunner_munzel.formula = function(formula,
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
  y <- do.call("brunner_munzel", c(DATA, list(...)))
  y$data.name <- DNAME
  y

}


perm_loop <-function(x,y,n.x,n.y,max_n_perm){

  pl1P<-matrix(0,nrow=n.x,ncol=max_n_perm)
  pl2P<-matrix(0,nrow=n.y,ncol=max_n_perm)

  for(h1 in 1:n.x){
    help1<-matrix(t(x[h1,]),
                  ncol=max_n_perm,
                  nrow=n.y,byrow=TRUE)
    pl1P[h1,]<-1/n.y*(colSums((y<help1)+1/2*(y==help1)))
  }
  for(h2 in 1:n.y){
    help2<-matrix(t(y[h2,]),
                  ncol=max_n_perm,
                  nrow=n.x,byrow=TRUE)
    pl2P[h2,]<-1/n.x*(colSums((x<help2)+1/2*(x==help2)))
  }

  pdP<-colMeans(pl2P)
  pd2P<-colMeans(pl1P)

  v1P<-(colSums(pl1P^2)-n.x*pd2P^2)/(n.x-1)
  v2P<-(colSums(pl2P^2)-n.y*pdP^2)/(n.y-1)
  vP<-v1P/n.x + v2P/n.y

  v0P<-(vP==0)
  vP[v0P]<-0.5/(n.x*n.y)^2

  res1<-matrix(rep(0,max_n_perm*3),nrow=3)

  res1[1,]<-(pdP-1/2)/sqrt(vP)

  pdP0<-(pdP==0)
  pdP1<-(pdP==1)
  pdP[pdP0]<-0.01
  pdP[pdP1]<-0.99

  res1[2,]<-log(pdP/(1-pdP))*pdP*(1-pdP)/sqrt(vP)
  res1[3,]<-qnorm(pdP)*exp(-0.5*qnorm(pdP)^2)/(sqrt(2*pi*vP))

  res1
}

