#' TOST function for t-tests
#' @param x a (non-empty) numeric vector of data values.
#' @param y an optional (non-empty) numeric vector of data values.
#' @param formula a formula of the form lhs ~ rhs where lhs is a numeric variable giving the data values and rhs either 1 for a one-sample or paired test or a factor with two levels giving the corresponding groups. If lhs is of class "Pair" and rhs is 1, a paired test is done.
#' @param data an optional matrix or data frame (or similar: see model.frame) containing the variables in the formula formula. By default the variables are taken from environment(formula).
#' @param paired a logical indicating whether you want a paired t-test.
#' @param low_eqbound lower equivalence bounds
#' @param high_eqbound upper equivalence bounds
#' @param eqbound_type Type of equivalance bound. Can be set to "SMD" for standardized mean difference (i.e., Cohen's d) or  "raw" for the mean difference. Default is "raw".
#' @param alpha alpha level (default = 0.05)
#' @param ...  further arguments to be passed to or from methods.
#' @return Returns TOSTt result object
#' @examples
#' ## TO BE ADDED
#' @importFrom stats pnorm pt qnorm qt
#' @importFrom graphics abline plot points segments title
#' @export


TOSTt = function(x, y = NULL,
                 formula,
                 data,
                 hypothesis = "EQU",
                 paired = FALSE,
                 var.equal = FALSE,
                 low_eqbound,
                 high_eqbound,
                 eqbound_type = "raw",
                 alpha = 0.05, ...){

}

