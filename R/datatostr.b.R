
#' @importFrom stats cor
dataTOSTrClass <- R6::R6Class(
  "dataTOSTrClass",
  inherit = dataTOSTrBase,
  private = list(
    .init = function() {

      ci <- 100 - as.integer(self$options$alpha * 200)

      tt <- self$results$tost
      #eqb <- self$results$eqb
      desc <- self$results$desc
      plots <- self$results$plots

      for (pair in self$options$pairs) {
        tt$setRow(rowKey=pair,  list(i1=pair[[1]], i2=pair[[2]]))
        #eqb$setRow(rowKey=pair, list(i1=pair[[1]], i2=pair[[2]]))
        desc$setRow(rowKey=pair, list(`name[1]`=pair[[1]], `name[2]`=pair[[2]]))
        plots$get(key=pair)$setTitle(paste(pair[[1]], '-', pair[[2]]))
      }

      tt$getColumn('cil')$setSuperTitle(jmvcore::format('{}% Confidence Interval', ci))
      tt$getColumn('ciu')$setSuperTitle(jmvcore::format('{}% Confidence Interval', ci))

    },
    .run = function() {

      tt <- self$results$tost
      #eqb <- self$results$eqb
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
        cor_type = self$options$cor_type
        #r <- suppressWarnings({
        #  stats::cor(i1, i2,
        #             method = cor_type)
        #})
        res <- suppressWarnings({
          stats::cor.test(i1, i2,
                          method = cor_type)
        })

        ci_r = cor_to_ci(cor = res$estimate,
                         n = n,
                         ci = 1-alpha*2,
                         method = cor_type,
                         correction = "fieller")

        if(unname(res$p.value) < alpha){
          sig_res = TRUE
        } else {
          sig_res = FALSE
        }

        #z1<-((log((1+r)/(1-r))/2)-(log((1+low_eqbound_r)/(1-low_eqbound_r))/2))/(sqrt(1/(n-3)))
        #z2<-((log((1+r)/(1-r))/2)-(log((1+high_eqbound_r)/(1-high_eqbound_r))/2))/(sqrt(1/(n-3)))

        if(self$options$hypothesis == "EQU"){
          test_hypothesis = "Hypothesis Tested: Equivalence"
          if(ci_r[1] > low_eqbound_r && ci_r[2] < high_eqbound_r){
            tost_res = TRUE
          } else {
            tost_res = FALSE
          }

        } else if(self$options$hypothesis == "MET"){
          test_hypothesis = "Hypothesis Tested: Minimal Effect"
          if(ci_r[2] < low_eqbound_r || ci_r[1] > high_eqbound_r){
            tost_res = TRUE
          } else {
            tost_res = FALSE
          }
        }


        #ptost<-max(p1,p2)
        #pttest<-2*(1-pt(abs(r)*sqrt(n-2)/sqrt(1-abs(r)^2),n-2))
        #zLL90<-(log((1+r)/(1-r))/2)-qnorm(1-alpha)*sqrt(1/(n-3))
        #zUL90<-(log((1+r)/(1-r))/2)+qnorm(1-alpha)*sqrt(1/(n-3))
        #LL90<-(exp(1)^(2*zLL90)-1)/(exp(1)^(2*zLL90)+1)
        #UL90<-(exp(1)^(2*zUL90)-1)/(exp(1)^(2*zUL90)+1)

        tt$setRow(rowKey=pair, list(
          `r` = unname(res$estimate),
          `p` = unname(res$p.value),
          `cil` = ci_r[1],
          `ciu` = ci_r[2],
          `sig` = as.character(sig_res),
          `tost` = as.character(tost_res)))

        #eqb$setRow(rowKey=pair, list(
        #  `low`=low_eqbound_r, `high`=high_eqbound_r, `cil`=LL90, `ciu`=UL90))

        desc$setRow(rowKey=pair, list(
          `n[1]`=n, `m[1]`=m1, `med[1]`=med1, `sd[1]`=sd1, `se[1]`=se1,
          `n[2]`=n, `m[2]`=m2, `med[2]`=med2, `sd[2]`=sd2, `se[2]`=se2))

        if(cor_type == "pearson"){
          r_lab = "Pearson's r"
        } else if(cor_type == "spearman"){
          r_lab ="Spearman's rho"
        } else if(cor_type == "kendall"){
          r_lab = "Kendall's tau"
        }

        text_res = paste0("Correlation method: ", r_lab,
                          "\n \n",
                          test_hypothesis,
                          "\n \n",
                          "Bounds: ",low_eqbound_r, " - ", high_eqbound_r)
        self$results$text$setContent(text_res)

        plot <- plots$get(key=pair)
        points <- list(r = unname(res$estimate),
                       n = n,
                       cor_type = cor_type,
                       low_eqbound_r = low_eqbound_r,
                       high_eqbound_r = high_eqbound_r,
                       r_lab = r_lab,
                       ciLL = ci_r[1],
                       ciUL = ci_r[2]
                       )
        plot$setState(points)
      }
    },
    .plot=function(image, ...) {

      if (is.null(image$state))
        return(FALSE)

      points = image$state

      res = plot_cor(r = points$r,
               n = points$n,
               method = points$cor_type,
               type = "cd") +
        geom_point(data = data.frame(y = 0,
                                     x = points$r),
                   aes(x = x, y = y),
                   size = 3) +
        annotate("segment",
                 x = points$ciLL,
                 xend = points$ciUL,
                 y = 0, yend = 0,
                 size = 1.5,
                 colour = "black")+
        geom_vline(aes(xintercept = points$low_eqbound_r),
                   linetype="dashed") +
        geom_vline(aes(xintercept = points$high_eqbound_r),
                   linetype="dashed") +
        scale_x_continuous(sec.axis = dup_axis(breaks=c(round(points$low_eqbound_r,2),
                                                        round(points$high_eqbound_r,2)))) +
        labs(y = "",x="")+
        facet_grid(~as.character(points$r_lab)) +
        theme_tidybayes() +
        theme(
          legend.position = "top",
          strip.text = element_text(face = "bold", size = 11),
          legend.text = element_text(face = "bold", size = 11),
          legend.title = element_text(face = "bold", size = 11),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(face = "bold", size = 11),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent",colour = NA),
          legend.background = element_rect(fill = "transparent",colour = NA)
        )

      print(res)

      return(TRUE)
    })
)
