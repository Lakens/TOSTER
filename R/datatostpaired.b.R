
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
      indplot <- self$results$indplot
      diffplot <- self$results$diffplot

      pair = paste0(self$options$pair1, " - ",self$options$pair2)
      plots$setTitle(pair)

      effsize$getColumn('cil')$setSuperTitle(jmvcore::format('{}% Confidence Interval', ci))
      effsize$getColumn('ciu')$setSuperTitle(jmvcore::format('{}% Confidence Interval', ci))

    },
    .run = function() {

      if (is.null(self$options$pair1) || length(self$options$pair2) == 0)
        return()

      tt <- self$results$tost
      eqb <- self$results$eqb
      effsize <- self$results$effsize
      desc <- self$results$desc
      plots <- self$results$plots
      indplot <- self$results$indplot
      diffplot <- self$results$diffplot

      pair = paste0(self$options$pair1, " - ",self$options$pair2)

      pair1Name <- paste(self$options$pair1)
      pair2Name <- paste(self$options$pair2)

      i1 <- jmvcore::toNumeric(self$data[pair1Name])
      i2 <- jmvcore::toNumeric(self$data[pair2Name])
      data <- data.frame(i1 = i1, i2 = i2)
      data <- na.omit(data)
      colnames(data) = c("i1","i2")
      data2 =  data
      data2$diff = data2$i2-data2$i1

      n <- nrow(data)
      i1 <- data$i1
      i2 <- data$i2
      m1 <- base::mean(i1)
      m2 <- base::mean(i2)
      med1 <- stats::median(i1)
      med2 <- stats::median(i2)
      sd1  <- stats::sd(i1)
      sd2  <- stats::sd(i2)
      se1  <- sd1 / sqrt(n)
      se2  <- sd2 / sqrt(n)

      res <- t.test(
        x = data$i1,
        y = data$i2,
        paired = TRUE,
        mu = 0,
        alternative = "two.sided"
      )

      t <- unname(res$statistic)
      p <- unname(res$p.value)
      df <- unname(res$parameter)

      alpha <- self$options$alpha
      low_eqbound    <- self$options$low_eqbound
      high_eqbound   <- self$options$high_eqbound

      low_eqbound_dz <- self$options$low_eqbound_dz  # deprecated
      high_eqbound_dz <- self$options$high_eqbound_dz


      r12 <- stats::cor(i1, i2)
      sdif <- sqrt(sd1 ^ 2 + sd2 ^ 2 - 2 * r12 * sd1 * sd2)

      cohend = abs(unname(res$estimate)) / sdif
      cohen_df = 2*(n-1)
      J <- gamma(df / 2) / (sqrt(df / 2) * gamma((df - 1) / 2))

      if (self$options$smd_type == 'g') {
        cohend <-  cohend * J
        smd_label = "Hedges' g(z)"
      } else {
        smd_label = "Cohen's d(z)"
      }

      d_lambda <- cohend * sqrt(n / (2*(1 - r12)))

      # Equation 4b Goulet-Pelletier and Cousineau, 2018
      # ((2*(1-r12))/n)
      d_sigma = sqrt((df + 1) / (df - 1) * ((2*(1-r12))/n) * (1 + cohend ^ 2 /8))

      # Get cohend confidence interval
      tlow <- suppressWarnings(qt(1 / 2 - (1 - alpha * 2) / 2, df = df, ncp = d_lambda))
      thigh <- suppressWarnings(qt(1 / 2 + (1 - alpha * 2) / 2, df = df, ncp = d_lambda))
      if (m1 == m2) {
        dlow <- (tlow * (2 * n)) / (sqrt(df) * sqrt(2 * n))
        dhigh <- (thigh * (2 * n)) / (sqrt(df) * sqrt(2 * n))
      } else {
        dlow <- tlow / d_lambda * cohend
        dhigh <- thigh / d_lambda * cohend
      }

      if ((m1 - m2) < 0) {
        cohend <- cohend * -1
        tdlow <- dlow
        dlow <- dhigh * -1
        dhigh <- tdlow * -1
        d_lambda <- cohend * sqrt(n / (2*(1 - r12)))
      }

      if (low_eqbound_dz != -999999999 &&
          low_eqbound_dz != -999999999) {
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

      if (self$options$hypothesis == "EQU") {
        alt_low = "greater"
        alt_high = "less"
      } else if (self$options$hypothesis == "MET") {
        alt_low = "less"
        alt_high = "greater"
      }

      low_ttest <- t.test(
        x = data$i1,
        y = data$i2,
        paired = TRUE,
        alternative = alt_low,
        mu = low_eqbound
      )

      high_ttest <- t.test(
        x = data$i1,
        y = data$i2,
        paired = TRUE,
        alternative = alt_high,
        mu = high_eqbound
      )

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

      ttost <- ifelse(abs(t1) < abs(t2), t1, t2)
      LL90 <- ((m1 - m2) - qt(1 - alpha, degree_f) * se)
      UL90 <- ((m1 - m2) + qt(1 - alpha, degree_f) * se)
      ptost <- max(p1, p2)
      dif <- (m1 - m2)
      LL95 <- ((m1 - m2) - qt(1 - (alpha / 2), degree_f) * se)
      UL95 <- ((m1 - m2) + qt(1 - (alpha / 2), degree_f) * se)

      tt$setRow(
        rowNo=1,
        values = list(
          `i1` = pair1Name,
          `i2` = pair2Name,
          `b` = "t-test",
          `t` = t,
          `df` = df,
          `p` = p
        )
      )

      tt$setRow(
        rowNo=2,
        values = list(
          `i1` = "",
          `i2` = "",
          `b` = "TOST Lower",
          `t` = t2,
          `df` = degree_f,
          `p` = p2
        )
      )

      tt$setRow(
        rowNo=3,
        values = list(
          `i1` = "",
          `i2` = "",
          `b` = "TOST Upper",
          `t` = t1,
          `df` = degree_f,
          `p` = p1
        )
      )

      eqb$setRow(
        rowNo=1,
        values = list(
          `stat` = smd_label,
          `low` = low_eqbound_dz,
          `high` = high_eqbound_dz
        )
      )

      eqb$setRow(
        rowNo=2,
        values = list(
          `stat` = "Raw",
          `low` = low_eqbound,
          `high` = high_eqbound
        )
      )

      effsize$setRow(
        rowNo = 1,
        values = list(
          `stat` = smd_label,
          `est` = cohend,
          `cil` = dlow,
          `ciu` = dhigh
        )
      )

      effsize$setRow(
        rowNo = 2,
        values = list(
          `stat` = "Raw",
          `est` = dif,
          `cil` = LL90,
          `ciu` = UL90
        )
      )

      desc$setRow(
        rowNo = 1,
        values = list(
          `name` = pair1Name,
          `n` = n,
          `m` = m1,
          `med` = med1,
          `sd` = sd1,
          `se` = se1
        )
      )

      desc$setRow(
        rowNo = 2,
        values = list(
          `name` = pair2Name,
          `n` = n,
          `m` = m2,
          `med` = med2,
          `sd` = sd2,
          `se` = se2
        )
      )

      #plot <- plots$get(key = 1)
      points <- data.frame(
        type = c("Mean Difference", smd_label),
        mu = c(dif, 0),
        param = c(round(unname(res$parameter), 0), round(unname(res$parameter), 0)),
        sigma = c(unname(res$stderr), d_sigma),
        lambda = c(0, d_lambda),
        low = c(low_eqbound, low_eqbound_dz),
        high = c(high_eqbound, high_eqbound_dz),
        alpha = c(alpha, alpha),
        stringsAsFactors = FALSE
      )

      #print(points)
      plots$setState(points)

      #indplot <- indplot$get(key = 1)
      colnames(data) = c(pair1Name,pair2Name)

      indplot$setState(data)

      diffplot$setState(data2)


      #}
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

      data2 = image$state

      pos_1 <- position_jitterdodge(
        jitter.width  = 0,
        jitter.height = 0.25,
        dodge.width   = 0.9
      )

      data_summary <- function(y) {
        m <- mean(y)
        ymin <- m-sd(y)
        ymax <- m+sd(y)
        return(c(x=m,xmin=xmin,xmax=xmax))
      }

      p2 = ggplot(data2,
                  aes(x = 1,
                      y = diff)) +
        geom_jitter(aes(x = .5),
                    alpha = 0.5,
                    height = 0,
                    width = 0.1) +
        stat_slab(
          fill = "steelblue2",
          side = "top",
          alpha = 0.5,
          color = "black"
        ) +
        stat_summary(aes(x = 0.9),
                     fun.data = data_summary,
                     size = 1.5) +
        theme_tidybayes() +
        labs(y = "Difference (Mean ± SD)",
             x = "") +
        theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank()) +
        coord_flip()

      print(p2)
      return(TRUE)
    },
    .indplot = function(image, ggtheme, theme, ...) {

      if (is.null(image$state))
        return(FALSE)

      data = image$state
      pairnames = colnames(data)
      dat_i1 = data[1]
      colnames(dat_i1) = c("val")
      dat_i1$pair = pairnames[1]
      dat_i1$id = row.names(dat_i1)
      dat_i2 = data[2]
      colnames(dat_i2) = c("val")
      dat_i2$pair = pairnames[2]
      dat_i2$id = row.names(dat_i2)

      dat = rbind(dat_i1,dat_i2)
      dat$pair = as.factor(dat$pair)


      data_summary <- function(x) {
        m <- mean(x)
        ymin <- m-sd(x)
        ymax <- m+sd(x)
        return(c(y=m, ymin=ymin, ymax=ymax))
      }


      dat$xj <- jitter(as.numeric(dat$pair), amount=.05)

      p = ggplot(data = dat, aes(y = val)) +
        geom_boxplot(aes(x = pair,
                         group = pair),
                     width = 0.2,
                     size = 1.5,
                     fatten = 1.5,
                     colour="steelblue2",
                     alpha = .4,
                     outlier.shape = NA) +
        geom_point(aes(x = xj),
                   alpha = .5) +
        geom_line(aes(x = xj,
                      group = id),
                  alpha = .5) +
        xlab("Condition") + ylab("") +
        #scale_x_continuous(breaks=c(1,2), labels=c("Before", "After"), limits=c(0.5, 2.5)) +
        theme_tidybayes()

      print(p)
      return(TRUE)
    })
)
