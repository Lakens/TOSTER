
#' @import ggplot2
#' @import ggdist
#' @import distributional
#' @importFrom psych harmonic.mean
dataTOSTpairedClass <- R6::R6Class(
  "dataTOSTpairedClass",
  inherit = dataTOSTpairedBase,
  private = list(
    .init = function() {

      ci <- 100 - as.integer(self$options$alpha * 200)

      tt <- self$results$tost
      eqb <- self$results$eqb
      desc <- self$results$desc
      plots <- self$results$plots

      for (pair in self$options$pairs) {
        tt$setRow(rowKey=pair,  list(i1=pair[[1]], i2=pair[[2]]))
        eqb$setRow(rowKey=pair, list(i1=pair[[1]], i2=pair[[2]]))
        desc$setRow(rowKey=pair, list(`name[1]`=pair[[1]], `name[2]`=pair[[2]]))
        plots$get(key=pair)$setTitle(paste(pair[[1]], '-', pair[[2]]))
      }

      eqb$getColumn('cil[cohen]')$setSuperTitle(jmvcore::format('{}% Confidence interval', ci))
      eqb$getColumn('ciu[cohen]')$setSuperTitle(jmvcore::format('{}% Confidence interval', ci))
      eqb$getColumn('cil[raw]')$setSuperTitle(jmvcore::format('{}% Confidence interval', ci))
      eqb$getColumn('ciu[raw]')$setSuperTitle(jmvcore::format('{}% Confidence interval', ci))
    },
    .run = function() {

      tt <- self$results$tost
      eqb <- self$results$eqb
      desc <- self$results$desc
      plots <- self$results$plots

      for (pair in self$options$pairs) {

        if (is.null(pair[[1]]))
          next()
        if (is.null(pair[[2]]))
          next()

        i1 <- jmvcore::toNumeric(self$data[[pair[[1]] ]])
        i2 <- jmvcore::toNumeric(self$data[[pair[[2]] ]])
        data <- data.frame(i1=i1, i2=i2)
        data <- na.omit(data)
        n <- nrow(data)
        i1 <- data$i1
        i2 <- data$i2
        m1 <- base::mean(i1)
        m2 <- base::mean(i2)
        med1 <- stats::median(i1)
        med2 <- stats::median(i2)
        sd1  <- stats::sd(i1)
        sd2  <- stats::sd(i2)
        se1  <- sd1/sqrt(n)
        se2  <- sd2/sqrt(n)

        res <- t.test(i1, i2, paired=TRUE)
        t <- unname(res$statistic)
        p <- unname(res$p.value)
        df <- unname(res$parameter)

        alpha <- self$options$alpha
        low_eqbound    <- self$options$low_eqbound
        high_eqbound   <- self$options$high_eqbound

        low_eqbound_dz <- self$options$low_eqbound_dz  # deprecated
        high_eqbound_dz <- self$options$high_eqbound_dz

        r12 <- stats::cor(i1, i2)
        sdif<-sqrt(sd1^2+sd2^2-2*r12*sd1*sd2)

        cohend = abs(m1-m2) / sdif

        J <- gamma(df/2)/(sqrt(df/2)*gamma((df-1)/2))



        if(self$options$smd_type == 'g') {
          cohend <-  cohend * J
          smd_label = "Hedges' g(z)"
        } else {
          smd_label = "Cohen's d(z)"
        }

        d_lambda <- cohend * sqrt(n/(2(1-r12)))
        d_sigma = sqrt((df + 1)/(df - 1)*(2/n)*(1 + cohend^2/8))

        tlow <- qt(1 / 2 - (1-alpha*2) / 2, df = df, ncp = d_lambda)
        thigh <- qt(1 / 2 + (1-alpha*2) / 2, df = df, ncp = d_lambda)
        if(m1 == m2) {
          dlow <- (tlow*(2*n)) / (sqrt(df) * sqrt(2*n))
          dhigh <- (thigh*(2*n)) / (sqrt(df) * sqrt(2*n))
        } else {
          dlow <- tlow / d_lambda * cohend
          dhigh <- thigh / d_lambda * cohend
        }

        if (low_eqbound_dz != -999999999 && low_eqbound_dz != -999999999) {
          # low_eqbound_dz and high_eqbound_dz options are deprecated
          low_eqbound  <- low_eqbound_d * sdif
          high_eqbound <- high_eqbound_d * sdif
        }
        else if (self$options$eqbound_type == 'd') {
          low_eqbound_dz <- low_eqbound
          high_eqbound_dz <- high_eqbound
          low_eqbound  <- low_eqbound * sdif
          high_eqbound <- high_eqbound * sdif
        } else {
          low_eqbound_dz <- low_eqbound / sdif
          high_eqbound_dz <- high_eqbound / sdif
        }

        low_ttest <- t.test(dep ~ group,
                            dataTTest,
                            paired = FALSE,
                            alternative = alt_low,
                            mu = low_eqbound)
        high_ttest <- t.test(dep ~ group,
                             dataTTest,
                             paired = FALSE,
                             alternative = alt_high,
                             mu = high_eqbound)

        t1 = low_ttest$statistic
        p1 = low_ttest$p.value
        t2 = high_ttest$statistic
        p2 = high_ttest$p.value

        degree_f = res$parameter
        pttest = res$p.value

        se <- res$stderr
        SE_val = res$stderr
        t <- res$statistic

        pttest <- res$p.value

        ttost<-ifelse(abs(t1) < abs(t2), t1, t2)
        LL90<-((m1-m2)-qt(1-alpha, degree_f)*se)
        UL90<-((m1-m2)+qt(1-alpha, degree_f)*se)
        ptost<-max(p1,p2)
        dif<-(m1-m2)
        LL95<-((m1-m2)-qt(1-(alpha/2), degree_f)*se)
        UL95<-((m1-m2)+qt(1-(alpha/2), degree_f)*se)

        tt$setRow(rowKey=pair, list(
          `t[0]`=t,  `df[0]`=df,       `p[0]`=p,
          `t[1]`=t2, `df[1]`=degree_f, `p[1]`=p2,
          `t[2]`=t1, `df[2]`=degree_f, `p[2]`=p1))

        eqb$setRow(rowKey=pair, list(
          `low[raw]`=low_eqbound, `high[raw]`=high_eqbound, `cil[raw]`=LL90, `ciu[raw]`=UL90,
          `low[cohen]`=low_eqbound_dz, `high[cohen]`=high_eqbound_dz))

        desc$setRow(rowKey=pair, list(
          `n[1]`=n, `m[1]`=m1, `med[1]`=med1, `sd[1]`=sd1, `se[1]`=se1,
          `n[2]`=n, `m[2]`=m2, `med[2]`=med2, `sd[2]`=sd2, `se[2]`=se2))

        plot <- plots$get(key=pair)
        points <- data.frame(
          m=dif,
          degree_f = unname(degree_f),
          SE= unname(SE_val),
          cil=LL90,
          ciu=UL90,
          low=low_eqbound,
          high=high_eqbound,
          stringsAsFactors=FALSE)
        plot$setState(points)
      }
    },
    .plot=function(image, ggtheme, theme, ...) {

      if (is.null(image$state))
        return(FALSE)

      tostplot(image, ggtheme, theme)

      return(TRUE)
    })
)
