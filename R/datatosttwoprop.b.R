
datatosttwopropClass <- if (requireNamespace('jmvcore')) R6::R6Class(
    "datatosttwopropClass",
    inherit = datatosttwopropBase,
    private = list(
        .init = function() {

          tt <- self$results$tost
          eqb <- self$results$eqb
          desc <- self$results$desc

          varName <- self$options$var
          levelName <- self$options$level
          if ( ! is.null(varName))
            tt$setNote('level', paste0("Testing proportion difference for \'", varName, " = ", levelName, "\'"))

          ci <- 100 - as.integer(self$options$alpha * 200)
          eqb$getColumn('cil')$setSuperTitle(jmvcore::format('{}% Confidence interval', ci))
          eqb$getColumn('ciu')$setSuperTitle(jmvcore::format('{}% Confidence interval', ci))

          groupName <- self$options$group
          groups <- NULL
          if ( ! is.null(groupName))
            groups <- base::levels(self$data[[groupName]])
          if (length(groups) != 2)
            groups <- c('Group 1', 'Group 2')

          desc <- self$results$desc
          desc$setRow(rowNo=1, values=list(
              `name[1]`=groups[1],
              `name[2]`=groups[2]))
        },
        .run = function() {

          if (is.null(self$options$group) || is.null(self$options$var))
            return()

          tt <- self$results$tost
          eqb <- self$results$eqb
          desc <- self$results$desc
          plot <- self$results$plot

          group <- self$data[[self$options$group]]
          group <- as.factor(group)
          group <- droplevels(group)

          groupLevels <- base::levels(group)
          if (length(groupLevels) != 2)
            jmvcore::reject("Grouping variable must have exactly 2 levels", code="grouping_var_must_have_2_levels")

          dep <- self$data[[self$options$var]]
          data <- data.frame(dep=dep, group=group)
          data <- na.omit(data)

          level <- self$options$level

          counts <- tapply(data$dep, data$group, function(x) sum(x == level))
          ns <- tapply(data$dep, data$group, function(x) length(x))
          props <- counts / ns

          prop_dif <- as.numeric(props[1] - props[2])


          prop_se <- as.numeric(sqrt((props[1]*(1-props[1]))/ns[1] + (props[2]*(1-props[2]))/ns[2]))

          low_eqbound <- self$options$low_eqbound
          high_eqbound <- self$options$high_eqbound
          alpha <- self$options$alpha


          # calculating z-statistic
          z1 <- (prop_dif - low_eqbound)/prop_se
          z2 <- (prop_dif - high_eqbound)/prop_se
          z  <- prop_dif / prop_se
          ztest <- 1 - pnorm(abs(z))

          # calculating p-value for both one-sided tests
          # Old just equivalence test
          #p1 <- 1 - pnorm(z1)
          #p2 <- pnorm(z2)

          if(self$options$hypothesis == "EQU"){
              p1 <- 1 - pnorm(z1)
              p2 <- pnorm(z2)
              hyp = "equivalence"
              test_hypothesis = "Hypothesis Tested: Equivalence"
          } else if(self$options$hypothesis == "MET"){
              p1 <- pnorm(z1)
              p2 <- 1-pnorm(z2)
              hyp = "minimal.effect"
              test_hypothesis = "Hypothesis Tested: Minimal Effect"
          }

          nhst_test = twoprop_test(
              p1 = props[1],
              p2 = props[2],
              n1 = ns[1],
              n2 = ns[2],
              null = 0,
              alpha = alpha,
              alternative = "t",
              effect_size = "difference"
          )
          tost_test = twoprop_test(
              p1 = props[1],
              p2 = props[2],
              n1 = ns[1],
              n2 = ns[2],
              null = c(low_eqbound, high_eqbound),
              alpha = alpha,
              alternative = hyp,
              effect_size = "difference"
          )

          # calculating CIs
          CI_lb <- prop_dif - (qnorm(1-alpha) * prop_se)
          CI_ub <- prop_dif + (qnorm(1-alpha) * prop_se)
          CI_lb95 <- prop_dif - (qnorm(1-(alpha/2)) * prop_se)
          CI_ub95 <- prop_dif + (qnorm(1-(alpha/2)) * prop_se)

          decision = ifelse(tost_test$p.value < alpha, " At the desired error rate, it can be stated that the ",
                            " At the desired error rate, it cannot be stated that the ")

          if (tost_test$alternative == "equivalence") {
              alt_state = paste0("true ", names(tost_test$null.value)[1],
                                 " is ", "between", " ", tost_test$null.value[1],
                                 " and ", tost_test$null.value[2])
          }
          if (tost_test$alternative == "minimal.effect") {
              alt_state = paste0("true ", names(tost_test$null.value)[1],
                                 " is ", "less than ", tost_test$null.value[1],
                                 " or ", "greater than ", tost_test$null.value[2])
          }

          text_res = paste0(test_hypothesis,
                            "<br> <br>",
                            decision, alt_state,". ", "<br>")
          self$results$text$setContent(text_res)
          #TOSTER::describe_htest()
          tt$setRow(rowNo=1, list(
            `z[0]`=nhst_test$statistic,  `p[0]`=nhst_test$p.value,
            `z[1]`=tost_test$statistic, `p[1]`=tost_test$p.value))

          eqb$setRow(rowNo=1, list(
            `low`=low_eqbound, `high`=high_eqbound, `cil`=CI_lb, `ciu`=CI_ub))

          desc$setRow(rowNo=1, list(
            `count[1]`=counts[1], `n[1]`=ns[1], `prop[1]`=props[1],
            `count[2]`=counts[2], `n[2]`=ns[2], `prop[2]`=props[2]))

          points <- data.frame(
            m=prop_dif,
            cil=CI_lb,
            ciu=CI_ub,
            low=low_eqbound,
            high=high_eqbound,
            stringsAsFactors=FALSE)
          plot$setState(points)

        },
        .plot = function(image, ggtheme, theme, ...) {

          if (is.null(image$state))
            return(FALSE)

          tostplot(image, ggtheme, theme)

          return(TRUE)
        })
)
