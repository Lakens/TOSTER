
# This file is automatically generated, you probably don't want to edit this

datatosttwopropOptions <- if (requireNamespace("jmvcore", quietly=TRUE)) R6::R6Class(
    "datatosttwopropOptions",
    inherit = jmvcore::Options,
    public = list(
        initialize = function(
            var = NULL,
            level = NULL,
            group = NULL,
            hypothesis = "EQU",
            low_eqbound = -0.1,
            high_eqbound = 0.1,
            alpha = 0.05,
            desc = FALSE,
            plot = FALSE, ...) {

            super$initialize(
                package="TOSTER",
                name="datatosttwoprop",
                requiresData=TRUE,
                ...)

            private$..var <- jmvcore::OptionVariable$new(
                "var",
                var,
                suggested=list(
                    "nominal"),
                permitted=list(
                    "factor"))
            private$..level <- jmvcore::OptionLevel$new(
                "level",
                level,
                variable="(var)")
            private$..group <- jmvcore::OptionVariable$new(
                "group",
                group,
                suggested=list(
                    "nominal"))
            private$..hypothesis <- jmvcore::OptionList$new(
                "hypothesis",
                hypothesis,
                options=list(
                    "EQU",
                    "MET"),
                default="EQU")
            private$..low_eqbound <- jmvcore::OptionNumber$new(
                "low_eqbound",
                low_eqbound,
                default=-0.1)
            private$..high_eqbound <- jmvcore::OptionNumber$new(
                "high_eqbound",
                high_eqbound,
                default=0.1)
            private$..alpha <- jmvcore::OptionNumber$new(
                "alpha",
                alpha,
                min=0,
                max=1,
                default=0.05)
            private$..desc <- jmvcore::OptionBool$new(
                "desc",
                desc,
                default=FALSE)
            private$..plot <- jmvcore::OptionBool$new(
                "plot",
                plot,
                default=FALSE)

            self$.addOption(private$..var)
            self$.addOption(private$..level)
            self$.addOption(private$..group)
            self$.addOption(private$..hypothesis)
            self$.addOption(private$..low_eqbound)
            self$.addOption(private$..high_eqbound)
            self$.addOption(private$..alpha)
            self$.addOption(private$..desc)
            self$.addOption(private$..plot)
        }),
    active = list(
        var = function() private$..var$value,
        level = function() private$..level$value,
        group = function() private$..group$value,
        hypothesis = function() private$..hypothesis$value,
        low_eqbound = function() private$..low_eqbound$value,
        high_eqbound = function() private$..high_eqbound$value,
        alpha = function() private$..alpha$value,
        desc = function() private$..desc$value,
        plot = function() private$..plot$value),
    private = list(
        ..var = NA,
        ..level = NA,
        ..group = NA,
        ..hypothesis = NA,
        ..low_eqbound = NA,
        ..high_eqbound = NA,
        ..alpha = NA,
        ..desc = NA,
        ..plot = NA)
)

datatosttwopropResults <- if (requireNamespace("jmvcore", quietly=TRUE)) R6::R6Class(
    "datatosttwopropResults",
    inherit = jmvcore::Group,
    active = list(
        text = function() private$.items[["text"]],
        tost = function() private$.items[["tost"]],
        eqb = function() private$.items[["eqb"]],
        desc = function() private$.items[["desc"]],
        plot = function() private$.items[["plot"]]),
    private = list(),
    public=list(
        initialize=function(options) {
            super$initialize(
                options=options,
                name="",
                title="TOST Two Proportions")
            self$add(jmvcore::Html$new(
                options=options,
                name="text",
                refs=list(
                    "TOSTER")))
            self$add(jmvcore::Table$new(
                options=options,
                name="tost",
                title="Hypothesis Tests Results",
                rows=1,
                clearWith=list(
                    "level",
                    "alpha",
                    "low_eqbound",
                    "high_eqbound",
                    "var",
                    "group"),
                columns=list(
                    list(
                        `name`="b[0]", 
                        `title`="", 
                        `type`="text", 
                        `content`="NHST"),
                    list(
                        `name`="z[0]", 
                        `title`="Z", 
                        `type`="number"),
                    list(
                        `name`="p[0]", 
                        `title`="p", 
                        `type`="number", 
                        `format`="zto,pvalue"),
                    list(
                        `name`="b[1]", 
                        `title`="", 
                        `type`="text", 
                        `content`="TOST"),
                    list(
                        `name`="z[1]", 
                        `title`="Z", 
                        `type`="number"),
                    list(
                        `name`="p[1]", 
                        `title`="p", 
                        `type`="number", 
                        `format`="zto,pvalue"))))
            self$add(jmvcore::Table$new(
                options=options,
                name="eqb",
                title="Equivalence Bounds",
                rows=1,
                clearWith=list(
                    "level",
                    "alpha",
                    "low_eqbound",
                    "high_eqbound",
                    "var",
                    "group"),
                columns=list(
                    list(
                        `name`="low", 
                        `title`="Low", 
                        `type`="number"),
                    list(
                        `name`="high", 
                        `title`="High", 
                        `type`="number"),
                    list(
                        `name`="estimate", 
                        `title`="Estimate"),
                    list(
                        `name`="cil", 
                        `title`="Lower", 
                        `superTitle`="Confidence interval"),
                    list(
                        `name`="ciu", 
                        `title`="Upper", 
                        `superTitle`="Confidence interval"))))
            self$add(jmvcore::Table$new(
                options=options,
                name="desc",
                title="Descriptives",
                visible="(desc)",
                rows=1,
                clearWith=list(
                    "level",
                    "var",
                    "group"),
                columns=list(
                    list(
                        `name`="name[1]", 
                        `title`="", 
                        `type`="text"),
                    list(
                        `name`="count[1]", 
                        `title`="Count", 
                        `type`="number"),
                    list(
                        `name`="n[1]", 
                        `title`="Total", 
                        `type`="number"),
                    list(
                        `name`="prop[1]", 
                        `title`="Proportion", 
                        `type`="number", 
                        `format`="zto"),
                    list(
                        `name`="name[2]", 
                        `title`="", 
                        `type`="text"),
                    list(
                        `name`="count[2]", 
                        `title`="Count", 
                        `type`="number"),
                    list(
                        `name`="n[2]", 
                        `title`="Total", 
                        `type`="number"),
                    list(
                        `name`="prop[2]", 
                        `title`="Proportion", 
                        `type`="number", 
                        `format`="zto"))))
            self$add(jmvcore::Image$new(
                options=options,
                name="plot",
                title="Plot",
                width=180,
                renderFun=".plot",
                visible="(plot)",
                clearWith=list(
                    "level",
                    "alpha",
                    "low_eqbound",
                    "high_eqbound",
                    "var",
                    "group")))}))

datatosttwopropBase <- if (requireNamespace("jmvcore", quietly=TRUE)) R6::R6Class(
    "datatosttwopropBase",
    inherit = jmvcore::Analysis,
    public = list(
        initialize = function(options, data=NULL, datasetId="", analysisId="", revision=0) {
            super$initialize(
                package = "TOSTER",
                name = "datatosttwoprop",
                version = c(1,0,0),
                options = options,
                results = datatosttwopropResults$new(options=options),
                data = data,
                datasetId = datasetId,
                analysisId = analysisId,
                revision = revision,
                pause = NULL,
                completeWhenFilled = FALSE,
                requiresMissings = FALSE,
                weightsSupport = 'auto')
        }))

#' TOST Two Proportions
#'
#' TOST Two Proportions for jamovi. This function is not meant to be utilized 
#' in R.
#' @param data .
#' @param var .
#' @param level .
#' @param group .
#' @param hypothesis \code{'EQU'} for equivalence (default), or \code{'MET'}
#'   for minimal effects test, the alternative hypothesis.
#' @param low_eqbound a number (default: -0.1) the lower equivalence bounds
#' @param high_eqbound a number (default: 0.1) the upper equivalence bounds
#' @param alpha alpha level (default = 0.05)
#' @param desc \code{TRUE} or \code{FALSE} (default), provide descriptive
#'   statistics
#' @param plot \code{TRUE} or \code{FALSE} (default), provide plot
#' @return A results object containing:
#' \tabular{llllll}{
#'   \code{results$text} \tab \tab \tab \tab \tab a html \cr
#'   \code{results$tost} \tab \tab \tab \tab \tab a table \cr
#'   \code{results$eqb} \tab \tab \tab \tab \tab a table \cr
#'   \code{results$desc} \tab \tab \tab \tab \tab a table \cr
#'   \code{results$plot} \tab \tab \tab \tab \tab an image \cr
#' }
#'
#' Tables can be converted to data frames with \code{asDF} or \code{\link{as.data.frame}}. For example:
#'
#' \code{results$tost$asDF}
#'
#' \code{as.data.frame(results$tost)}
#'
#' @export
datatosttwoprop <- function(
    data,
    var,
    level,
    group,
    hypothesis = "EQU",
    low_eqbound = -0.1,
    high_eqbound = 0.1,
    alpha = 0.05,
    desc = FALSE,
    plot = FALSE) {

    if ( ! requireNamespace("jmvcore", quietly=TRUE))
        stop("datatosttwoprop requires jmvcore to be installed (restart may be required)")

    if ( ! missing(var)) var <- jmvcore::resolveQuo(jmvcore::enquo(var))
    if ( ! missing(group)) group <- jmvcore::resolveQuo(jmvcore::enquo(group))
    if (missing(data))
        data <- jmvcore::marshalData(
            parent.frame(),
            `if`( ! missing(var), var, NULL),
            `if`( ! missing(group), group, NULL))

    for (v in var) if (v %in% names(data)) data[[v]] <- as.factor(data[[v]])

    options <- datatosttwopropOptions$new(
        var = var,
        level = level,
        group = group,
        hypothesis = hypothesis,
        low_eqbound = low_eqbound,
        high_eqbound = high_eqbound,
        alpha = alpha,
        desc = desc,
        plot = plot)

    analysis <- datatosttwopropClass$new(
        options = options,
        data = data)

    analysis$run()

    analysis$results
}

