
#' @import ggplot2
#' @import ggdist
#' @import distributional

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

      effsize$getColumn('cil[cohen]')$setSuperTitle(jmvcore::format('{}% Confidence Interval', ci))
      effsize$getColumn('ciu[cohen]')$setSuperTitle(jmvcore::format('{}% Confidence Interval', ci))
      effsize$getColumn('cil[raw]')$setSuperTitle(jmvcore::format('{}% Confidence Interval', ci))
      effsize$getColumn('ciu[raw]')$setSuperTitle(jmvcore::format('{}% Confidence Interval', ci))

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

      if ( ! self$options$var_equal){
        tt$setNote('var_equal', "Welch's t-test")
      }


      if ( ! self$options$var_equal){
        effsize$setNote('var_equal', "Denominator set to the average SD")
      } else {
        effsize$setNote('var_equal', "Denominator set to the pooled SD")
      }

    },
    .run = function() {

      if (is.null(self$options$group) || length(self$options$deps) == 0)
        return()

      tt <- self$results$tost
      eqb <- self$results$eqb
      effsize <- self$results$effsize
      desc <- self$results$desc
      plots <- self$results$plots
      descplot <- self$results$descplots
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

        if(self$options$smd_type == "g"){
          bias_c = TRUE
        } else {
          bias_c = FALSE
        }

        low_eqbound    <- self$options$low_eqbound
        high_eqbound   <- self$options$high_eqbound

        TOSTres = t_TOST(formula = dep ~ group,
                         data = dataTTest,
                         paired = FALSE,
                         eqbound_type = self$options$eqbound_type,
                         hypothesis = self$options$hypothesis,
                         var.equal = var.equal,
                         alpha = alpha,
                         bias_correction = bias_c,
                         low_eqbound = low_eqbound,
                         high_eqbound = high_eqbound)

        if (self$options$eqbound_type == 'SMD') {

          pr_l_eqb = low_eqbound * TOSTres$smd$d_denom
          pr_h_eqb = high_eqbound * TOSTres$smd$d_denom
        } else if(self$options$eqbound_type == 'raw') {

          pr_l_eqb = low_eqbound
          pr_h_eqb = high_eqbound
        }


        if(self$options$hypothesis == "EQU"){
          alt_low = "greater"
          alt_high = "less"
          test_hypothesis = "Hypothesis Tested: Equivalence"
          null_hyp = paste0(round(pr_l_eqb,2),
                            " &ge; (Mean<sub>1</sub> - Mean<sub>2</sub>) or (Mean<sub>1</sub> - Mean<sub>2</sub>) &ge; ",
                            round(pr_h_eqb,2))
          alt_hyp = paste0(round(pr_l_eqb,2),
                           " < (Mean<sub>1</sub> - Mean<sub>2</sub>) < ",
                           round(pr_h_eqb,2))
        } else if(self$options$hypothesis == "MET"){
          alt_low = "less"
          alt_high = "greater"
          test_hypothesis = "Hypothesis Tested: Minimal Effect"
          null_hyp = paste0(round( pr_l_eqb,2),
                            " &le; (Mean<sub>1</sub> - Mean<sub>2</sub>)  &le; ",
                            round(pr_h_eqb,2))
          alt_hyp = paste0(round( pr_l_eqb,2),
                           " > (Mean<sub>1</sub> - Mean<sub>2</sub>) or (Mean<sub>1</sub> - Mean<sub>2</sub>)  > ",
                           round(pr_h_eqb,2))
        }


        tt$setRow(
          rowKey = depName,
          list(
            `t[0]` = TOSTres$TOST$t[1],
            `df[0]` = TOSTres$TOST$df[1],
            `p[0]` = TOSTres$TOST$p.value[1],
            `t[1]` = TOSTres$TOST$t[2],
            `df[1]` = TOSTres$TOST$df[2],
            `p[1]` = TOSTres$TOST$p.value[2],
            `t[2]` = TOSTres$TOST$t[3],
            `df[2]` = TOSTres$TOST$df[3],
            `p[2]` = TOSTres$TOST$p.value[3]
          )
        )

        eqb$setRow(
          rowKey = depName,
          list(
            `stat[cohen]` = TOSTres$smd$smd_label,
            `low[cohen]` = TOSTres$eqb$low_eq[2],
            `high[cohen]` = TOSTres$eqb$high_eq[2],
            `stat[raw]` = "Raw",
            `low[raw]` = TOSTres$eqb$low_eq[1],
            `high[raw]` = TOSTres$eqb$high_eq[1]
          )
        )

        effsize$setRow(
          rowKey = depName,
          list(
            `stat[cohen]` = TOSTres$smd$smd_label,
            `est[cohen]` = TOSTres$smd$d,
            `cil[cohen]` = TOSTres$smd$dlow,
            `ciu[cohen]` = TOSTres$smd$dhigh,
            `est[raw]` = TOSTres$effsize$estimate[1],
            `cil[raw]` = TOSTres$effsize$lower.ci[1],
            `ciu[raw]` = TOSTres$effsize$upper.ci[1]
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

        if(grepl(TOSTres$decision$ttest, pattern="non")){
          nhst_text = "&#10060 NHST: don't reject null significance hypothesis that the effect is equal to zero"
        } else{
          nhst_text = "&#9989 NHST: reject null significance hypothesis that the effect is equal to zero"

        }

        if(self$options$hypothesis == "EQU"){
          tost_hypt = "equivalence"
        } else{
          tost_hypt = "MET"
        }

        if(grepl(TOSTres$decision$TOST, pattern="non")){
          TOST_text = paste0("&#10060 TOST: don't reject null ", tost_hypt," hypothesis")
        } else{
          TOST_text = paste0("&#9989 TOST: reject null ", tost_hypt," hypothesis")

        }

        text_res = paste0(test_hypothesis,
                          "<br> <br>",
                          "Null Hypothesis: ", null_hyp,"<br>",
                          "Alternative: ", alt_hyp,"<br>",
                          nhst_text, "<br>",
                          TOST_text, "<br>",
                          "<br> Note: SMD confidence intervals are an approximation. See our <a href=\"https://aaroncaldwell.us/TOSTERpkg/articles/SMD_calcs.html\">documentation</a>. <br>",
                          ifelse(self$options$eqbound_type == 'SMD',
                                 "<br> &#x1f6a8; Warning: standardized bounds produce biased results. Consider setting bounds in raw units.", ""))
        self$results$text$setContent(text_res)

        plot <- plots$get(key=depName)

        plot$setState(TOSTres)

        descplot <- descplot$get(key=depName)

        descplot$setState(dataTTest)
        #print(points)
      }
    },
    .plot=function(image, ggtheme, theme, ...) {

      if (is.null(image$state))
        return(FALSE)

      TOSTres <- image$state

      plotTOSTr = plot_tost_jam(TOSTres,
                                ggtheme = ggtheme)
      print(plotTOSTr)

      return(TRUE)
    },
    .descplot = function(image, ggtheme, theme, ...) {

      if (is.null(image$state))
        return(FALSE)

      dataTTest = image$state


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
             aes(x = group,
                 y = dep,
                 color = group)) +
        geom_jitter(alpha = 0.5, position = pos_1) +
        stat_slab(aes(x=as.numeric(group)+.2),
                  fill="lightgrey",
                  side = "right") +
        stat_summary(aes(x=as.numeric(group)+.2),
                     fun.data=data_summary) +
        # theme_tidybayes() +
        labs(x="Group",
             y="",
             color = "Group")  +
        #scale_colour_manual(values=c("red2","dodgerblue")) +
        ggtheme +
        theme(
          legend.position = "top"
        )

      print(p)

      return(TRUE)
    })
)
