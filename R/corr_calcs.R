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

.corboot <- function(isub, x, y, ...){
  res <- cor(x[isub], y[isub], ...)
  res
}
