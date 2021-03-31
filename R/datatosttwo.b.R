
#' @import ggplot2
#' @import ggdist
#' @import distributional
#' @importFrom psych harmonic.mean
dataTOSTtwoClass <- R6::R6Class(
  "dataTOSTtwoClass",
  inherit = dataTOSTtwoBase,
  private = list(
    .init = function() {

      ci <- 100 - as.integer(self$options$alpha * 200)

      tt <- self$results$tost
      eqb <- self$results$eqb
      effsize <- self$results$effsize
      desc <- self$results$desc

      effsize$getColumn('cil[cohen]')$setSuperTitle(jmvcore::format('{}% Confidence interval', ci))
      effsize$getColumn('ciu[cohen]')$setSuperTitle(jmvcore::format('{}% Confidence interval', ci))
      effsize$getColumn('cil[raw]')$setSuperTitle(jmvcore::format('{}% Confidence interval', ci))
      effsize$getColumn('ciu[raw]')$setSuperTitle(jmvcore::format('{}% Confidence interval', ci))

      groupName <- self$options$group
      groups <- NULL
      if ( ! is.null(groupName))
        groups <- base::levels(self$data[[groupName]])
      if (length(groups) != 2)
        groups <- c('Group 1', 'Group 2')

      desc <- self$results$desc
      for (key in desc$rowKeys) {
        desc$setRow(rowKey=key, values=list(
          `name[1]`=groups[1],
          `name[2]`=groups[2]))
      }

      if ( ! self$options$var_equal)
        tt$setNote('var_equal', "Welch's t-test")
    },
    .run = function() {

      if (is.null(self$options$group) || length(self$options$deps) == 0)
        return()

      tt <- self$results$tost
      eqb <- self$results$eqb
      effsize <- self$results$effsize
      desc <- self$results$desc
      plots <- self$results$plots
      old = FALSE

      groupName <- self$options$group
      group <- self$data[[groupName]]
      group <- as.factor(group)
      group <- droplevels(group)

      groupLevels <- base::levels(group)
      if (length(groupLevels) != 2)
        jmvcore::reject("Grouping variable must have exactly 2 levels", code="grouping_var_must_have_2_levels")

      for (depName in self$options$deps) {

        dep <- self$data[[depName]]
        dep <- jmvcore::toNumeric(dep)
        dataTTest <- data.frame(dep=dep, group=group)
        dataTTest <- na.omit(dataTTest)

        v <- tapply(dataTTest$dep, dataTTest$group, function(x) jmvcore::tryNaN(var(x)))
        n <- tapply(dataTTest$dep, dataTTest$group, length)
        m <- tapply(dataTTest$dep, dataTTest$group, function(x) jmvcore::tryNaN(mean(x)))
        med <- tapply(dataTTest$dep, dataTTest$group, function(x) jmvcore::tryNaN(median(x)))
        se <- sqrt(v/n)
        sd <- sqrt(v)

        sediff <- jmvcore::tryNaN(sqrt((v[1]/n[1])+(v[2]/n[2])))

        n[is.na(n)] <- 0
        m[is.na(m)] <- NaN
        med[is.na(med)] <- NaN
        se[is.na(se)] <- NaN
        sd[is.na(sd)] <- NaN
        sediff[is.na(sediff)] <- NaN



        var.equal <- self$options$var_equal
        alpha <- self$options$alpha

        low_eqbound    <- self$options$low_eqbound
        high_eqbound   <- self$options$high_eqbound
        low_eqbound_d  <- self$options$low_eqbound_d  # deprecated
        high_eqbound_d <- self$options$high_eqbound_d

        tresult <- t.test(dep ~ group,
                          dataTTest,
                          paired = FALSE,
                          alternative = "two.sided",
                          mu = 0,
                          var.equal=var.equal)

        t <- unname(tresult$statistic)
        p <- unname(tresult$p.value)
        df <- unname(tresult$parameter)

        if (var.equal) {
          denomSD <- sqrt((((n[1] - 1)*(sd[1]^2)) + (n[2] - 1)*(sd[2]^2))/((n[1]+n[2])-2)) #calculate sd pooled
        } else {
          denomSD <- sqrt((sd[1]^2 + sd[2]^2)/2) #calculate sd root mean squared for Welch's t-test
        }

        denomSD[is.na(denomSD)] <- NaN

        #denomSD <- jmvcore::tryNaN(sqrt(((n[1]-1)*v[1]+(n[2]-1)*v[2])/(n[1]+n[2]-2)))
        d <- abs(m[1]-m[2])/denomSD # Cohen's d

        d[is.na(d)] <- NaN

        cohend = d
        ntilde <- harmonic.mean(c(n[1],n[2]))

        # Compute unbiased Hedges' g
        # Use the lgamma function, and update to what Goulet-Pelletier & cousineau used; works with larger inputs
        J <- gamma(df/2)/(sqrt(df/2)*gamma((df-1)/2))

        if(self$options$smd_type == 'g') {
          cohend <-  cohend * J
          smd_label = "Hedges' g"
        } else {
          smd_label = "Cohen's d"
        }

        # add options for cohend here
        d_lambda <- cohend * sqrt(ntilde/2)
        d_sigma = sqrt((n[1]+n[2])/(n[1]*n[2])+(cohend^2/(2*(n[1]+n[2]))))

        # Confidence interval of the SMD from Goulet-Pelletier & Cousineau
        tlow <- qt(1 / 2 - (1-alpha*2) / 2, df = df, ncp = d_lambda)
        thigh <- qt(1 / 2 + (1-alpha*2) / 2, df = df, ncp = d_lambda)
        if(m[1] == m[2]) {
          dlow <- (tlow*(n[1]+n[2])) / (sqrt(df) * sqrt(n[1]*n[2]))
          dhigh <- (thigh*(n[1]+n[2])) / (sqrt(df) * sqrt(n[1]*n[2]))
        } else {
          dlow <- tlow / d_lambda * cohend
          dhigh <- thigh / d_lambda * cohend
        }

        # The function provided by Goulet-Pelletier & Cousineau works with +mdiff.
        # Now we fix the signs and directions for negative differences
        if ((m[1]-m[2]) < 0) {
          cohend <- cohend * -1
          tdlow <- dlow
          dlow <- dhigh * -1
          dhigh <- tdlow * -1
        }



        if (low_eqbound_d != -999999999 && low_eqbound_d != -999999999) {
          # low_eqbound_d and high_eqbound_d options are deprecated
          low_eqbound  <- low_eqbound_d * denomSD
          high_eqbound <- high_eqbound_d * denomSD
        }
        else if (self$options$eqbound_type == 'd') {
          low_eqbound_d <- low_eqbound
          high_eqbound_d <- high_eqbound
          low_eqbound  <- low_eqbound * denomSD
          high_eqbound <- high_eqbound * denomSD
        } else {
          low_eqbound_d <- low_eqbound / denomSD
          high_eqbound_d <- high_eqbound / denomSD
        }

        if(self$options$hypothesis == "EQU"){
          alt_low = "greater"
          alt_high = "less"
        } else if(self$options$hypothesis == "MET"){
          alt_low = "less"
          alt_high = "greater"
        }

        low_ttest <- t.test(dep ~ group,
                            dataTTest,
                            paired = FALSE,
                            var.equal = var.equal,
                            alternative = alt_low,
                            mu = low_eqbound)

        high_ttest <- t.test(dep ~ group,
                             dataTTest,
                             paired = FALSE,
                             var.equal = var.equal,
                             alternative = alt_high,
                             mu = high_eqbound)

        t1 = low_ttest$statistic
        p1 = low_ttest$p.value
        t2 = high_ttest$statistic
        p2 = high_ttest$p.value

        degree_f = tresult$parameter
        pttest = tresult$p.value

        LL90 <- (m[1]-m[2])-qt(1-alpha, degree_f)*tresult$stderr
        UL90 <- (m[1]-m[2])+qt(1-alpha, degree_f)*tresult$stderr
        LL95 <- (m[1]-m[2])-qt(1-(alpha/2), degree_f)*tresult$stderr
        UL95 <- (m[1]-m[2])+qt(1-(alpha/2), degree_f)*tresult$stderr

        dif<-(m[1]-m[2])

        tt$setRow(rowKey=depName, list(
          `t[0]`=t,  `df[0]`=df, `p[0]`=p,
          `t[1]`=t2, `df[1]`=degree_f, `p[1]`=p2,
          `t[2]`=t1, `df[2]`=degree_f, `p[2]`=p1))

        eqb$setRow(
          rowKey = depName,
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
          rowKey = depName,
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

        desc$setRow(
          rowKey = depName,
          list(
            `n[1]` = n[1],
            `m[1]` = m[1],
            `med[1]` = med[1],
            `sd[1]` = sd[1],
            `se[1]` = se[1],
            `n[2]` = n[2],
            `m[2]` = m[2],
            `med[2]` = med[2],
            `sd[2]` = sd[2],
            `se[2]` = se[2]
          )
        )

        # Get SE value
        SE_val = tresult$stderr

        plot <- plots$get(key=depName)
        points <- data.frame(
          type = c("Mean Difference", smd_label),
          mu = c(dif,0),
          #param = c(round(unname(tresult$parameter),round(n[1]+n[2]-2,0))),
          param = c(round(unname(tresult$parameter),0),round((n[1]+n[2]-2),0)),
          #param = c(60,90),
          sigma = c(unname(tresult$stderr), d_sigma),
          lambda = c(0, d_lambda),
          #cil=LL90,
          #ciu=UL90,
          low=c(low_eqbound, low_eqbound_d),
          high=c(high_eqbound, high_eqbound_d),
          alpha = c(alpha, alpha),
          stringsAsFactors=FALSE)
        plot$setState(points)
        #print(points)
      }
    },
    .plot=function(image, ggtheme, theme, ...) {

      if (is.null(image$state))
        return(FALSE)

      cut_cdf_qi = function(p, .width = c(.66, .95, 1), labels = NULL) {
        .width = sort(.width)

        if (is.function(labels)) {
          labels = labels(.width)
        } else if (is.null(labels)) {
          labels = .width
        }

        cut(abs(1 - p*2), labels = labels, breaks = c(0, .width), include.lowest = TRUE, ordered_result = TRUE)
      }

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
