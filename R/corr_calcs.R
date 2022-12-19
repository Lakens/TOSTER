cor_to_ci <- function(cor, n, ci = 0.95,
                      method = "pearson",
                      correction = "fieller", ...) {
  method <- match.arg(tolower(method),
                      c("pearson", "kendall", "spearman"),
                      several.ok = FALSE)

  if (method == "kendall") {
    out <- .cor_to_ci_kendall(cor, n,
                              ci = ci,
                              correction = correction, ...)
  } else if (method == "spearman") {
    out <- .cor_to_ci_spearman(cor, n,
                               ci = ci,
                               correction = correction, ...)
  } else {
    out <- .cor_to_ci_pearson(cor, n, ci = ci, ...)
  }

  out
}

# Kendall -----------------------------------------------------------------
#' @importFrom stats qnorm
.cor_to_ci_kendall <- function(cor, n,
                               ci = 0.95,
                               correction = "fieller", ...) {
  # by @tsbaguley (https://rpubs.com/seriousstats/616206)

  if (correction == "fieller") {
    tau.se <- (0.437 / (n - 4))^0.5
  } else {
    tau.se <- 1 / (n - 3)^0.5
  }

  moe <- stats::qnorm(1 - (1 - ci) / 2) * tau.se
  zu <- atanh(cor) + moe
  zl <- atanh(cor) - moe

  # Convert back to r
  ci_low <- tanh(zl)
  ci_high <- tanh(zu)

  c(ci_low, ci_high)
}


# Spearman -----------------------------------------------------------------
.cor_to_ci_spearman <- function(cor, n,
                                ci = 0.95,
                                correction = "fieller", ...) {
  # by @tsbaguley (https://rpubs.com/seriousstats/616206)

  if (correction == "fieller") {
    zrs.se <- (1.06 / (n - 3))^0.5
  } else if (correction == "bw") {
    zrs.se <- ((1 + (cor^2) / 2) / (n - 3))^0.5
  } else {
    zrs.se <- 1 / (n - 3)^0.5
  }

  moe <- stats::qnorm(1 - (1 - ci) / 2) * zrs.se

  zu <- atanh(cor) + moe
  zl <- atanh(cor) - moe

  # Convert back to r
  ci_low <- tanh(zl)
  ci_high <- tanh(zu)

  c(ci_low, ci_high)
}


# Pearson -----------------------------------------------------------------
.cor_to_ci_pearson <- function(cor, n, ci = 0.95, ...) {
  z <- atanh(cor)
  se <- 1 / sqrt(n - 3) # Sample standard error

  # CI
  alpha <- 1 - (1 - ci) / 2
  ci_low <- z - se * stats::qnorm(alpha)
  ci_high <- z + se * stats::qnorm(alpha)

  # Convert back to r
  ci_low <- tanh(ci_low)
  ci_high <- tanh(ci_high)

  c(ci_low, ci_high)
}
# TOSTER:::cor_to_ci(.5,18,"pearson",ci=.95)

corr_curv = function (r,
                      n,
                      type = "pearson",
                      steps = 5000) {
  intrvls <- (0:steps)/steps
  intrvls = subset(intrvls,intrvls>0 & intrvls<1)


  results <-
    suppressWarnings({
      lapply(
        intrvls,
        FUN = function(i)
          cor_to_ci(
            cor = r, n = n,
            method = type,
            correction = "fieller",
            ci = i
          )
      )
    })

  df <- data.frame(do.call(rbind, results))
  intrvl.limit <- c("lower.limit", "upper.limit")
  colnames(df) <- intrvl.limit
  df$intrvl.width <- (abs((df$upper.limit) - (df$lower.limit)))
  df$intrvl.level <- intrvls
  df$cdf <- (abs(df$intrvl.level/2)) + 0.5
  df$pvalue <- 1 - intrvls
  df$svalue <- -log2(df$pvalue)
  df <- head(df, -1)
  class(df) <- c("data.frame", "concurve")
  densdf <- data.frame(c(df$lower.limit, df$upper.limit))
  colnames(densdf) <- "x"
  densdf <- head(densdf, -1)
  class(densdf) <- c("data.frame", "concurve")

  return(list(df, densdf))

}

# plot_cor(.5,18)

.corboot <- function(isub, x, y, method){
  res <- cor(x[isub], y[isub], method = method)
  res
}

# Other corrs -----


# winsorized

wincor <- function(x, y, tr=0.2){

  df <- cbind(x,y)
  df <- df[complete.cases(df), ]
  n <- nrow(df)
  x <- df[,1]
  y <- df[,2]
  g <- floor(tr*n)
  xvec <- winval(x,tr)
  yvec <- winval(y,tr)
  corv <- cor(xvec,yvec)

  return(corv)
}

winval <- function(x, tr=0.2){
  y <- sort(x)
  n <- length(x)
  ibot <- floor(tr*n)+1
  itop <- length(x)-ibot+1
  xbot <- y[ibot]
  xtop <- y[itop]
  winval <- ifelse(x<=xbot,xbot,x)
  winval <- ifelse(winval>=xtop,xtop,winval)
  winval
}

# percentage bend measure of location
pbos <- function(x, beta=.2){
  n <- length(x)
  temp <- sort(abs(x-median(x)))
  omhatx <- temp[floor((1-beta)*n)]
  psi <- (x-median(x))/omhatx
  i1 <- length(psi[psi<(-1)])
  i2 <- length(psi[psi>1])
  sx <- ifelse(psi<(-1),0,x)
  sx <- ifelse(psi>1,0,sx)
  pbos <- (sum(sx) + omhatx*(i2-i1)) / (n-i1-i2)
  pbos
}

pbcor <- function(x, y, beta=.2){

  df <- cbind(x,y)
  df <- df[complete.cases(df), ]
  n <- nrow(df)
  x <- df[,1]
  y <- df[,2]
  temp <- sort(abs(x-median(x)))
  omhatx <- temp[floor((1-beta)*n)]
  temp <- sort(abs(y-median(y)))
  omhaty <- temp[floor((1-beta)*n)]
  a <- (x-pbos(x,beta))/omhatx
  b <- (y-pbos(y,beta))/omhaty
  a <- ifelse(a<=-1,-1,a)
  a <- ifelse(a>=1,1,a)
  b <- ifelse(b<=-1,-1,b)
  b <- ifelse(b>=1,1,b)
  corv <- sum(a*b)/sqrt(sum(a^2)*sum(b^2))

  return(corv)
}

.corboot_wincor <- function(isub, x, y, ...){
  res <- wincor(x[isub], y[isub], ...)
  res
}

.corboot_pbcor <- function(isub, x, y, ...){
  res <- pbcor(x[isub], y[isub], ...)
  res
}
