
#' @import ggplot2
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
        low_eqbound_dz <- self$options$low_eqbound_dz
        high_eqbound_dz <- self$options$high_eqbound_dz
        r12 <- stats::cor(i1, i2)

        sdif<-sqrt(sd1^2+sd2^2-2*r12*sd1*sd2)
        low_eqbound<-low_eqbound_dz*sdif
        high_eqbound<-high_eqbound_dz*sdif
        se<-sdif/sqrt(n)
        t<-(m1-m2)/se
        degree_f<-n-1
        pttest<-2*pt(abs(t), degree_f, lower.tail=FALSE)
        t1<-((m1-m2)+(low_eqbound_dz*sdif))/se
        p1<-1-pt(t1, degree_f, lower.tail=FALSE)
        t2<-((m1-m2)+(high_eqbound_dz*sdif))/se
        p2<-pt(t2, degree_f, lower.tail=FALSE)
        ttost<-ifelse(abs(t1) < abs(t2), t1, t2)
        LL90<-((m1-m2)-qt(1-alpha, degree_f)*se)
        UL90<-((m1-m2)+qt(1-alpha, degree_f)*se)
        ptost<-max(p1,p2)
        dif<-(m1-m2)
        LL95<-((m1-m2)-qt(1-(alpha/2), degree_f)*se)
        UL95<-((m1-m2)+qt(1-(alpha/2), degree_f)*se)

        tt$setRow(rowKey=pair, list(
          `t[0]`=t,  `df[0]`=df,       `p[0]`=p,
          `t[1]`=t1, `df[1]`=degree_f, `p[1]`=p1,
          `t[2]`=t2, `df[2]`=degree_f, `p[2]`=p2))

        eqb$setRow(rowKey=pair, list(
          `low[raw]`=low_eqbound, `high[raw]`=high_eqbound, `cil[raw]`=LL90, `ciu[raw]`=UL90,
          `low[cohen]`=low_eqbound_dz, `high[cohen]`=high_eqbound_dz))

        desc$setRow(rowKey=pair, list(
          `n[1]`=n, `m[1]`=m1, `med[1]`=med1, `sd[1]`=sd1, `se[1]`=se1,
          `n[2]`=n, `m[2]`=m2, `med[2]`=med2, `sd[2]`=sd2, `se[2]`=se2))

        plot <- plots$get(key=pair)
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
