
#' @import ggplot2
#' @import ggdist
#' @import distributional

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

      alpha <- self$options$alpha
      low_eqbound    <- self$options$low_eqbound
      high_eqbound   <- self$options$high_eqbound
      eqbound_type = self$options$eqbound_type

      if(self$options$smd_type == 'g'){
        bias_c = TRUE
      } else {
        bias_c = FALSE
      }

      TOSTres = t_TOST(x = data$i2,
                       y = data$i1,
                       paired = TRUE,
                       eqbound_type = eqbound_type,
                       hypothesis = self$options$hypothesis,
                       alpha = alpha,
                       bias_correction = bias_c,
                       low_eqbound = low_eqbound,
                       high_eqbound = high_eqbound)

      tt$setRow(
        rowNo=1,
        values = list(
          `i1` = pair1Name,
          `i2` = pair2Name,
          `b` = "t-test",
          `t` = TOSTres$TOST$t[1],
          `df` = TOSTres$TOST$df[1],
          `p` = TOSTres$TOST$p.value[1]
        )
      )

      tt$setRow(
        rowNo=2,
        values = list(
          `i1` = "",
          `i2` = "",
          `b` = "TOST Lower",
          `t` = TOSTres$TOST$t[2],
          `df` = TOSTres$TOST$df[2],
          `p` = TOSTres$TOST$p.value[2]
        )
      )

      tt$setRow(
        rowNo=3,
        values = list(
          `i1` = "",
          `i2` = "",
          `b` = "TOST Upper",
          `t` = TOSTres$TOST$t[3],
          `df` = TOSTres$TOST$df[3],
          `p` = TOSTres$TOST$p.value[3]
        )
      )

      eqb$setRow(
        rowNo=1,
        values = list(
          `stat` = TOSTres$smd$smd_label,
          `low` = TOSTres$eqb$low_eq[2],
          `high` = TOSTres$eqb$high_eq[2]
        )
      )

      eqb$setRow(
        rowNo=2,
        values = list(
          `stat` = "Raw",
          `low` = TOSTres$eqb$low_eq[1],
          `high` = TOSTres$eqb$high_eq[1]
        )
      )

      effsize$setRow(
        rowNo = 2,
        values = list(
          `stat` = TOSTres$smd$smd_label,
          `est` = TOSTres$smd$d,
          `cil` = TOSTres$smd$dlow,
          `ciu` = TOSTres$smd$dhigh
        )
      )

      effsize$setRow(
        rowNo = 1,
        values = list(
          `stat` = "Raw",
          `est` = TOSTres$effsize$estimate[1],
          `cil` = TOSTres$effsize$lower.ci[1],
          `ciu` = TOSTres$effsize$upper.ci[1]
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

      if (eqbound_type == 'SMD') {

        pr_l_eqb = low_eqbound * TOSTres$smd$d_denom
        pr_h_eqb = high_eqbound * TOSTres$smd$d_denom
      } else if(eqbound_type == "raw") {

        pr_l_eqb = low_eqbound
        pr_h_eqb = high_eqbound
      }

      if(self$options$hypothesis == "EQU"){
        alt_low = "greater"
        alt_high = "less"
        test_hypothesis = "Hypothesis Tested: Equivalence"
        null_hyp = paste0(round(pr_l_eqb,2),
                          " >= (Mean<sub>1</sub> - Mean<sub>2</sub>) or (Mean<sub>1</sub> - Mean<sub>2</sub>) >= ",
                          round(pr_h_eqb,2))
        alt_hyp = paste0(round(pr_l_eqb,2),
                         " < (Mean<sub>1</sub> - Mean<sub>2</sub>) < ",
                         round(pr_h_eqb,2))
      } else if(self$options$hypothesis == "MET"){
        alt_low = "less"
        alt_high = "greater"
        test_hypothesis = "Hypothesis Tested: Minimal Effect"
        null_hyp = paste0(round(pr_l_eqb,2),
                          " <= (Mean<sub>1</sub> - Mean<sub>2</sub>)  <= ",
                          round(pr_h_eqb,2))
        alt_hyp = paste0(round(pr_l_eqb,2),
                         " > (Mean<sub>1</sub> - Mean<sub>2</sub>) or (Mean<sub>1</sub> - Mean<sub>2</sub>)  > ",
                         round(pr_h_eqb,2))
      }

      if(self$options$hypothesis == "EQU"){
        tost_hypt = "equivalence"
      } else{
        tost_hypt = "MET"
      }

      if(grepl(TOSTres$decision$ttest, pattern="non")){
        nhst_text = "&#10060 NHST: don't reject null significance hypothesis that the effect is equal to zero"
      } else{
        nhst_text = "&#9989 NHST: reject null significance hypothesis that the effect is equal to zero"

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

      #print(points)
      plots$setState(TOSTres)

      #indplot <- indplot$get(key = 1)
      colnames(data) = c(pair1Name,pair2Name)

      indplot$setState(data)

      diffplot$setState(list(data2 = data2,
                             TOSTres = TOSTres))


      #}
    },
    .plot=function(image, ggtheme, theme, ...) {

      if (is.null(image$state))
        return(FALSE)

      TOSTres <- image$state

      plotTOSTr = plot(TOSTres) + ggtheme
      print(plotTOSTr)

      return(TRUE)
    },
    .diffplot = function(image, ggtheme, theme, ...) {

      if (is.null(image$state))
        return(FALSE)

      dats = image$state

      data2 = dats$data2

      TOSTres = dats$TOSTres

      eqb = TOSTres$eqb

      data_summary <- function(x) {
        m <- mean(x)
        ymin <- m-sd(x)
        ymax <- m+sd(x)
        return(c(y=m,ymin=ymin,ymax=ymax))
      }


      p2 = ggplot(data2,
                  aes(x = 1,
                      y = diff)) +
        geom_jitter(aes(x = .75),
                    alpha = 0.5,
                    height = 0,
                    width = 0.125) +
        stat_slab(
          fill = "steelblue2",
          side = "top",
          alpha = 0.5,
          color = "black"
        ) +
        stat_summary(aes(x=0.5),
                     fun.data=data_summary,
                     size = 1.5)  +
        labs(y = "Difference",
             x = "",
             caption = "Dot and Line represent Mean and SD respectively") +
        ggtheme+
        coord_flip()

      p3 = p2 +
        geom_hline(yintercept = eqb$low_eq[1],
                   linetype="dashed",
                   alpha = .5,
                   size = 1.5) +
        geom_hline(yintercept = eqb$high_eq[1],
                   linetype="dashed",
                   alpha = .5,
                   size = 1.5) +
        #issues with placement of text
       # geom_text(aes(x=2.5, y = dats$low_eqbound[1]),
      #            hjust = "outward",
      #            vjust=-1.0,
      #            angle = 90,
      #            label='Lower Bound') +
      #  geom_text(aes(x=2.5, y = dats$high_eqbound[1]),
      #            hjust = "inward",
      #            vjust=1.25,
      #            angle = 90,
      #            label='Upper Bound') +
        ggtheme +
        theme(
          legend.position = "top",
        )

      print(p3)
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
      dat$pair = factor(dat$pair,
                        levels = c(pairnames[1], pairnames[2]),
                        ordered = TRUE)


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
        ggtheme +
        theme(
          legend.position = "top"
        )

      print(p)
      return(TRUE)
    })
)
