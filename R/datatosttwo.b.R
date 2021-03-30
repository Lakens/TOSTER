
#' @import ggplot2
#' @import ggdist
#' @import distributional
#'
dataTOSTtwoClass <- R6::R6Class(
  "dataTOSTtwoClass",
  inherit = dataTOSTtwoBase,
  private = list(
    .init = function() {

      ci <- 100 - as.integer(self$options$alpha * 200)

      tt <- self$results$tost
      eqb <- self$results$eqb
      desc <- self$results$desc

      eqb$getColumn('cil[cohen]')$setSuperTitle(jmvcore::format('{}% Confidence interval', ci))
      eqb$getColumn('ciu[cohen]')$setSuperTitle(jmvcore::format('{}% Confidence interval', ci))
      eqb$getColumn('cil[raw]')$setSuperTitle(jmvcore::format('{}% Confidence interval', ci))
      eqb$getColumn('ciu[raw]')$setSuperTitle(jmvcore::format('{}% Confidence interval', ci))

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

        pooledSD <- jmvcore::tryNaN(sqrt(((n[1]-1)*v[1]+(n[2]-1)*v[2])/(n[1]+n[2]-2)))
        d <- (m[1]-m[2])/pooledSD # Cohen's d

        n[is.na(n)] <- 0
        m[is.na(m)] <- NaN
        med[is.na(med)] <- NaN
        se[is.na(se)] <- NaN
        sd[is.na(sd)] <- NaN
        sediff[is.na(sediff)] <- NaN
        pooledSD[is.na(pooledSD)] <- NaN
        d[is.na(d)] <- NaN


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
          sdpooled <- sqrt((((n[1] - 1)*(sd[1]^2)) + (n[2] - 1)*(sd[2]^2))/((n[1]+n[2])-2)) #calculate sd pooled
        } else {
          sdpooled <- sqrt((sd[1]^2 + sd[2]^2)/2) #calculate sd root mean squared for Welch's t-test
        }

        if (low_eqbound_d != -999999999 && low_eqbound_d != -999999999) {
          # low_eqbound_d and high_eqbound_d options are deprecated
          low_eqbound  <- low_eqbound_d * sdpooled
          high_eqbound <- high_eqbound_d * sdpooled
        }
        else if (self$options$eqbound_type == 'd') {
          low_eqbound_d <- low_eqbound
          high_eqbound_d <- high_eqbound
          low_eqbound  <- low_eqbound * sdpooled
          high_eqbound <- high_eqbound * sdpooled
        } else {
          low_eqbound_d <- low_eqbound / sdpooled
          high_eqbound_d <- high_eqbound / sdpooled
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

        eqb$setRow(rowKey=depName, list(
          `low[raw]`=low_eqbound, `high[raw]`=high_eqbound, `cil[raw]`=LL90, `ciu[raw]`=UL90,
          `low[cohen]`=low_eqbound_d, `high[cohen]`=high_eqbound_d))

        desc$setRow(rowKey=depName, list(
          `n[1]`=n[1], `m[1]`=m[1], `med[1]`=med[1], `sd[1]`=sd[1], `se[1]`=se[1],
          `n[2]`=n[2], `m[2]`=m[2], `med[2]`=med[2], `sd[2]`=sd[2], `se[2]`=se[2]))

        # Get SE value
        SE_val = tresult$stderr

        plot <- plots$get(key=depName)
        points <- data.frame(
          mu = dif,
          param = round(unname(tresult$parameter),0),
          sigma = unname(tresult$stderr),
          cil=LL90,
          ciu=UL90,
          low=low_eqbound,
          high=high_eqbound,
          alpha = alpha,
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
      c1 = 1-points$alpha
      c2 = 1-points$alpha*2
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
            mu = points$mu,
            df = points$param,
            sigma = points$sigma,
          ),
          fill = stat(cut_cdf_qi(p=cdf,
                                 .width = sets))
        ),
        .width = c(c2, c1)) +
        scale_fill_brewer(direction = -1,
                          na.translate = FALSE) +
        labs(x = 'Mean Difference', y = 'Density',
             fill = "Confidence Interval") +
        geom_vline(xintercept = low,linetype="dashed") +
        geom_vline(xintercept = high,linetype="dashed") +
        geom_text(aes(y=1.55, x=low,  vjust=-.9, hjust=1),
                  angle = 90,
                  label='Lower Bound') +
        geom_text(aes(y=1.55, x=high, vjust=1.5, hjust=1),
                  angle = 90,
                  label='Upper Bound') +
        theme_bw()
      print(plot)

      return(TRUE)
    })
)
