
#' @import ggplot2
dataTOSToneClass <- R6::R6Class(
  "dataTOSToneClass",
  inherit = dataTOSToneBase,
  private = list(
    .init = function() {

      ci <- 100 - as.integer(self$options$alpha * 200)

      tt <- self$results$tost
      eqb <- self$results$eqb
      effsize <- self$results$effsize

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
      effsize <- self$results$effsize
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
        res <- t.test(var-mu)
        t   <- unname(res$statistic)
        pttest   <- unname(res$p.value)

        low_eqbound    <- self$options$low_eqbound
        high_eqbound   <- self$options$high_eqbound
        low_eqbound_d  <- self$options$low_eqbound_d  # deprecated
        high_eqbound_d <- self$options$high_eqbound_d

        if (low_eqbound_d != -999999999 && low_eqbound_d != -999999999) {
          # low_eqbound_d and high_eqbound_d options are deprecated
          low_eqbound  <- low_eqbound_d * sd
          high_eqbound <- high_eqbound_d * sd
        }
        else if (self$options$eqbound_type == 'd') {
          low_eqbound_d <- low_eqbound
          high_eqbound_d <- high_eqbound
          low_eqbound  <- low_eqbound * sd
          high_eqbound <- high_eqbound * sd
        } else {
          low_eqbound_d <- low_eqbound / sd
          high_eqbound_d <- high_eqbound / sd
        }

        if(self$options$hypothesis == "EQU"){
          alt_low = "greater"
          alt_high = "less"
        } else if(self$options$hypothesis == "MET"){
          alt_low = "less"
          alt_high = "greater"
        }

        degree_f <- res$parameter
        SE_val <- res$stderr

        low_ttest <- t.test(x = var,
                            alternative = alt_low,
                            mu = low_eqbound)

        t1 = low_ttest$statistic
        p1 = low_ttest$p.value

        high_ttest <- t.test(x = var,
                             alternative = alt_high,
                             mu = high_eqbound)
        t2 = high_ttest$statistic
        p2 = high_ttest$p.value

        #t1<-(m-mu-low_eqbound)/(sd/sqrt(n))# t-test
        #p1<-pt(t1, degree_f, lower.tail=FALSE)
        #t2<-(m-mu-high_eqbound)/(sd/sqrt(n)) #t-test
        #p2<-pt(t2, degree_f, lower.tail=TRUE)
        #t<-(m-mu)/(sd/sqrt(n))
        #pttest<-2*pt(-abs(t), df=degree_f)
        LL90<-m-mu-qt(1-alpha, degree_f)*res$stderr
        UL90<-m-mu+qt(1-alpha, degree_f)*res$stderr
        LL95<-m-mu-qt(1-(alpha/2), degree_f)*res$stderr
        UL95<-m-mu+qt(1-(alpha/2), degree_f)*res$stderr
        ptost<-max(p1,p2) #Get highest p-value for summary TOST result
        ttost<-ifelse(abs(t1) < abs(t2), t1, t2) #Get lowest t-value for summary TOST result
        dif<-(m-mu)

        cohend <- abs(dif)/sd # Cohen's d
        df <- res$parameter
        # Compute unbiased Hedges' g
        # Use the lgamma function, and update to what Goulet-Pelletier & cousineau used; works with larger inputs
        J <- gamma(df/2)/(sqrt(df/2)*gamma((df-1)/2))

        if(self$options$smd_type == 'g') {
          cohend <-  cohend * J
          smd_label = "Hedges' g(z)"
        } else {
          smd_label = "Cohen's d(z)"
        }

        d_lambda <- cohend * sqrt(n/2)
        d_sigma = sqrt((df + 1)/(df - 1)*(2/n)*(1 + cohend^2/8))

        # Confidence interval of the SMD from Goulet-Pelletier & Cousineau
        tlow <- qt(1 / 2 - (1-alpha*2) / 2, df = df, ncp = d_lambda)
        thigh <- qt(1 / 2 + (1-alpha*2) / 2, df = df, ncp = d_lambda)
        if(dif == 0) {
          dlow <- (tlow*(2*n)) / (sqrt(df) * sqrt(2*n))
          dhigh <- (thigh*(2*n)) / (sqrt(df) * sqrt(2*n))
        } else {
          dlow <- tlow / d_lambda * cohend
          dhigh <- thigh / d_lambda * cohend
        }

        if ((dif) < 0) {
          cohend <- cohend * -1
          tdlow <- dlow
          dlow <- dhigh * -1
          dhigh <- tdlow * -1
        }


        tt$setRow(rowKey=name, list(
          `t[0]`=t,  `df[0]`=degree_f, `p[0]`=p,
          `t[1]`=t2, `df[1]`=degree_f, `p[1]`=p2,
          `t[2]`=t1, `df[2]`=degree_f, `p[2]`=p1))

        eqb$setRow(rowKey=name, list(
          `low[raw]` = low_eqbound,
          `high[raw]` = high_eqbound,
          `cil[raw]` = LL90,
          `ciu[raw]` = UL90,
          `low[cohen]` = low_eqbound_d,
          `high[cohen]` = high_eqbound_d))

        effsize$setRow(
          rowKey = pair,
          list(
            `stat[cohen]` = smd_label,
            `est[cohen]` = cohend,
            `cil[cohen]` = dlow,
            `ciu[cohen]` = dhigh,
            `est[raw]` = dif,
            `cil[raw]` = LL90,
            `ciu[raw]` = UL90
          )
        )

        desc$setRow(rowKey=name, list(n=n, m=m, med=med, sd=sd, se=se))

        plot <- plots$get(key=name)
        points <- data.frame(
          type = c("Mean Difference", smd_label),
          mu = c(dif,0),
          param = c(round(unname(res$parameter),0),round(unname(res$parameter),0)),
          sigma = c(unname(res$stderr), d_sigma),
          lambda = c(0, d_lambda),
          low=c(low_eqbound, low_eqbound_d),
          high=c(high_eqbound, high_eqbound_d),
          alpha = c(alpha, alpha),
          stringsAsFactors=FALSE)
        plot$setState(points)
      }

    },
    .plot=function(image, ggtheme, theme, ...) {

      if (is.null(image$state))
        return(FALSE)

      points <- image$state
      c1 = 1-points$alpha[1]
      c2 = 1-points$alpha[1]*2
      if(c1 < .999 && c2 > .5){
        sets = c(.5,c2,c1,.999)
      } else if(c2 <=.5 && c1 < .999) {
        sets = c(c2,c1,.999)
      } else {
        sets = c(.5,c2,c1)
      }
      plot = ggplot(data = points,
                    aes_string(y = 0)) +
        stat_dist_halfeye(aes(
          dist = dist_student_t(
            mu = mu,
            df = param,
            sigma = sigma,
            ncp = lambda
          ),
          fill = stat(cut_cdf_qi(p=cdf,
                                 .width = sets))
        ),
        .width = c(c2, c1)) +
        scale_fill_brewer(direction = -1,
                          na.translate = FALSE) +
        labs(x = '', y = '',
             fill = "Confidence Interval") +
        geom_vline(aes(xintercept = low),
                   linetype="dashed") +
        geom_vline(aes(xintercept = high),
                   linetype="dashed") +
        geom_text(aes(y=1.5, x=low,
                      vjust=-.9, hjust=1),
                  angle = 90,
                  label='Lower Bound') +
        geom_text(aes(y=1.5, x=high, vjust=1.5, hjust=1),
                  angle = 90,
                  label='Upper Bound') +
        theme_tidybayes() +
        theme(legend.position="top",
              strip.text = element_text(face="bold", size=10),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank()) +
        facet_wrap(~type,
                   ncol = 1,
                   scales = "free")
      print(plot)

      return(TRUE)
    })
)
