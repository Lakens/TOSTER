
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
      effsize <- self$results$effsize
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
      effsize <- self$results$effsize
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

        res <- t.test(x = i1,
                      y = i2,
                      paired=TRUE,
                      mu = 0,
                      alternative = "two.sided")

        t <- unname(res$statistic)
        p <- unname(res$p.value)
        df <- unname(res$parameter)

        alpha <- self$options$alpha
        low_eqbound    <- self$options$low_eqbound
        high_eqbound   <- self$options$high_eqbound

        low_eqbound_dz <- self$options$low_eqbound_dz  # deprecated
        high_eqbound_dz <- self$options$high_eqbound_dz

        r12 <- stats::cor(i1, i2)
        sdif <- sqrt(sd1^2+sd2^2-2*r12*sd1*sd2)

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

        # Get cohend confidence interval
        tlow <- qt(1 / 2 - (1-alpha*2) / 2, df = df, ncp = d_lambda)
        thigh <- qt(1 / 2 + (1-alpha*2) / 2, df = df, ncp = d_lambda)
        if(m1 == m2) {
          dlow <- (tlow*(2*n)) / (sqrt(df) * sqrt(2*n))
          dhigh <- (thigh*(2*n)) / (sqrt(df) * sqrt(2*n))
        } else {
          dlow <- tlow / d_lambda * cohend
          dhigh <- thigh / d_lambda * cohend
        }

        if ((m[1]-m[2]) < 0) {
          cohend <- cohend * -1
          tdlow <- dlow
          dlow <- dhigh * -1
          dhigh <- tdlow * -1
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

        eqb$setRow(
          rowKey = pair,
          list(
            `stat[cohen]` = smd_label,
            `low[cohen]` = low_eqbound_d,
            `high[cohen]` = high_eqbound_d,
            `stat[raw]` = "Raw",
            `low[raw]` = low_eqbound,
            `high[raw]` = high_eqbound
          )
        )

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

        plot <- plots$get(key=pair)
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
    },
    .diffplot = function(image, ggtheme, theme, ...) {

      if (is.null(image$state))
        return(FALSE)

      dep <- self$data[[depName]]
      dep <- jmvcore::toNumeric(dep)
      dataTTest <- data.frame(dep=dep, group=group)
      dataTTest <- na.omit(dataTTest)

      pos_1 <- position_jitterdodge(
        jitter.width  = 0.25,
        jitter.height = 0,
        dodge.width   = 0.9
      )

      data_summary <- function(x) {
        m <- mean(x)
        ymin <- m-sd(x)
        ymax <- m+sd(x)
        return(c(y=m,ymin=ymin,ymax=ymax))
      }

      p = ggplot(dataTTest,
                 aes(x=group,
                     y=dep,
                     color = group)) +
        geom_jitter(alpha = 0.5, position = pos_1) +
        stat_slab(aes(x=as.numeric(group)-.2),
                  fill="lightgrey",
                  side = "left") +
        stat_summary(aes(x=as.numeric(group)-.2),
                     fun.data=data_summary) +
        theme_tidybayes() +
        labs(x="Group",
             y="",
             color = "Group")  +
        scale_colour_manual(values=c("red2","dodgerblue"))

      return(p)
    },
    .indplot = function(image, ggtheme, theme, ...) {

      if (is.null(image$state))
        return(FALSE)

      i1 <- jmvcore::toNumeric(self$data[[pair[[1]] ]])
      i2 <- jmvcore::toNumeric(self$data[[pair[[2]] ]])
      data <- data.frame(i1=i1, i2=i2)
      data <- na.omit(data)
      dat_i1 = data$i1
      colnames(dat_i1) = c("val")
      dat_i1$pair = pair$key[1]
      dat_i1$id = row.names(dat_i1$id)
      dat_i2 = data$i2
      colnames(dat_i2) = c("val")
      dat_i2$pair = pair$key[2]
      dat_i2$id = row.names(dat_i2$id)

      dat = rbind(dat_i1,dat_i2)

      data_summary <- function(x) {
        m <- mean(x)
        ymin <- m-sd(x)
        ymax <- m+sd(x)
        return(c(y=m, ymin=ymin, ymax=ymax))
      }

      dat$xj <- jitter(dat$pair, amount=.03)

      ggplot(data=dat, aes(y=val)) +
        geom_boxplot(aes(x=pair,
                         group=pair),
                     width=0.2,
                     outlier.shape = NA) +
        geom_point(aes(x=xj)) +
        geom_line(aes(x=xj,
                      group=id)) +
        xlab("Condition") + ylab("Value") +
        #scale_x_continuous(breaks=c(1,2), labels=c("Before", "After"), limits=c(0.5, 2.5)) +
        theme_bw()

      return(p)
    })
)
