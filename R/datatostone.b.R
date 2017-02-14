
#' @import ggplot2
dataTOSToneClass <- R6::R6Class(
  "dataTOSToneClass",
  inherit = dataTOSToneBase,
  private = list(
    .init = function() {

      ci <- 100 - as.integer(self$options$alpha * 200)

      tt <- self$results$tost
      eqb <- self$results$eqb

      eqb$getColumn('cil[cohen]')$setSuperTitle(jmvcore::format('{}% Confidence interval', ci))
      eqb$getColumn('ciu[cohen]')$setSuperTitle(jmvcore::format('{}% Confidence interval', ci))
      eqb$getColumn('cil[raw]')$setSuperTitle(jmvcore::format('{}% Confidence interval', ci))
      eqb$getColumn('ciu[raw]')$setSuperTitle(jmvcore::format('{}% Confidence interval', ci))

      for (key in eqb$rowKeys)
        eqb$setRow(rowKey=key, values=list(`cil[cohen]`='', `ciu[cohen]`=''))
    },
    .run = function() {

      tt <- self$results$tost
      eqb <- self$results$eqb
      desc <- self$results$desc
      plots <- self$results$plots

      alpha <- self$options$alpha
      mu  <- self$options$mu

      for (name in self$options$vars) {
        var <- jmvcore::toNumeric(self$data[[name]])
        var <- na.omit(var)
        n   <- length(var)
        m   <- base::mean(var)
        med <- stats::median(var)
        sd  <- stats::sd(var)
        se  <- sd/sqrt(n)
        res <- t.test(var)
        t   <- unname(res$statistic)
        p   <- unname(res$p.value)

         low_eqbound_d <- self$options$lowEqBD
        high_eqbound_d <- self$options$highEqBD
         low_eqbound <-  low_eqbound_d * sd
        high_eqbound <- high_eqbound_d * sd

        degree_f<-n-1
        t1<-(m-mu-low_eqbound)/(sd/sqrt(n))# t-test
        p1<-pt(t1, degree_f, lower.tail=FALSE)
        t2<-(m-mu-high_eqbound)/(sd/sqrt(n)) #t-test
        p2<-pt(t2, degree_f, lower.tail=TRUE)
        t<-(m-mu)/(sd/sqrt(n))
        pttest<-2*pt(-abs(t), df=degree_f)
        LL90<-m-mu-qt(1-alpha, degree_f)*(sd/sqrt(n))
        UL90<-m-mu+qt(1-alpha, degree_f)*(sd/sqrt(n))
        LL95<-m-mu-qt(1-(alpha/2), degree_f)*(sd/sqrt(n))
        UL95<-m-mu+qt(1-(alpha/2), degree_f)*(sd/sqrt(n))
        ptost<-max(p1,p2) #Get highest p-value for summary TOST result
        ttost<-ifelse(abs(t1) < abs(t2), t1, t2) #Get lowest t-value for summary TOST result
        dif<-(m-mu)

        tt$setRow(rowKey=name, list(
          `t[0]`=t,  `df[0]`=degree_f, `p[0]`=p,
          `t[1]`=t1, `df[1]`=degree_f, `p[1]`=p1,
          `t[2]`=t2, `df[2]`=degree_f, `p[2]`=p2))

        eqb$setRow(rowKey=name, list(
          `low[raw]`=low_eqbound, `high[raw]`=high_eqbound, `cil[raw]`=LL90, `ciu[raw]`=UL90,
          `low[cohen]`=low_eqbound_d, `high[cohen]`=high_eqbound_d))

        desc$setRow(rowKey=name, list(n=n, m=m, med=med, sd=sd, se=se))

        plot <- plots$get(key=name)
        points <- data.frame(
          m=dif,
          cil=LL90,
          ciu=UL90,
          low=low_eqbound,
          high=high_eqbound,
          stringsAsFactors=FALSE)
        plot$setState(points)
      }

    },
    .plot=function(image, ...) {

      if (is.null(image$state))
        return(FALSE)

      tostplot(image)

      return(TRUE)
    })
)
