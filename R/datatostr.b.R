
#' @importFrom stats cor
dataTOSTrClass <- R6::R6Class(
  "dataTOSTrClass",
  inherit = dataTOSTrBase,
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

      eqb$getColumn('cil')$setSuperTitle(jmvcore::format('{}% Confidence interval', ci))
      eqb$getColumn('ciu')$setSuperTitle(jmvcore::format('{}% Confidence interval', ci))

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

        i1 <- jmvcore::toNumeric(self$data[[ pair[[1]] ]])
        i2 <- jmvcore::toNumeric(self$data[[ pair[[2]] ]])
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

        alpha <- self$options$alpha
        low_eqbound_r <- self$options$low_eqbound_r
        high_eqbound_r <- self$options$high_eqbound_r
        r <- stats::cor(i1, i2)
        res <- stats::cor.test(i1, i2)

        z1<-((log((1+abs(r))/(1-abs(r)))/2)-(log((1+low_eqbound_r)/(1-low_eqbound_r))/2))/(sqrt(1/(n-3)))
        z2<-((log((1+abs(r))/(1-abs(r)))/2)-(log((1+high_eqbound_r)/(1-high_eqbound_r))/2))/(sqrt(1/(n-3)))
        p1<-ifelse(low_eqbound_r<r,pnorm(-abs(z1)),1-pnorm(-abs(z1)))
        p2<-ifelse(high_eqbound_r>r,pnorm(-abs(z2)),1-pnorm(-abs(z2)))
        ptost<-max(p1,p2)
        pttest<-2*(1-pt(abs(r)*sqrt(n-2)/sqrt(1-abs(r)^2),n-2))
        zLL90<-(log((1+r)/(1-r))/2)-qnorm(1-alpha)*sqrt(1/(n-3))
        zUL90<-(log((1+r)/(1-r))/2)+qnorm(1-alpha)*sqrt(1/(n-3))
        LL90<-(exp(1)^(2*zLL90)-1)/(exp(1)^(2*zLL90)+1)
        UL90<-(exp(1)^(2*zUL90)-1)/(exp(1)^(2*zUL90)+1)

        tt$setRow(rowKey=pair, list(
          `r[0]`=r, `p[0]`=unname(res$p.value),
          `r[1]`=r, `p[1]`=p1,
          `r[2]`=r, `p[2]`=p2))

        eqb$setRow(rowKey=pair, list(
          `low`=low_eqbound_r, `high`=high_eqbound_r, `cil`=LL90, `ciu`=UL90))

        desc$setRow(rowKey=pair, list(
          `n[1]`=n, `m[1]`=m1, `med[1]`=med1, `sd[1]`=sd1, `se[1]`=se1,
          `n[2]`=n, `m[2]`=m2, `med[2]`=med2, `sd[2]`=sd2, `se[2]`=se2))

        plot <- plots$get(key=pair)
        points <- data.frame(
          m=r,
          cil=LL90,
          ciu=UL90,
          low=low_eqbound_r,
          high=high_eqbound_r,
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
