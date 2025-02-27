## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----warning=FALSE, message=FALSE---------------------------------------------
library(TOSTER)
# Get Data
data("InsectSprays")

# Look at the data structure
head(InsectSprays)

## ----warning=FALSE, message=FALSE---------------------------------------------
# Build ANOVA
aovtest = aov(count ~ spray,
              data = InsectSprays)

# Display overall results
knitr::kable(broom::tidy(aovtest),
            caption = "Traditional ANOVA Test")

## -----------------------------------------------------------------------------
equ_ftest(Fstat = 34.70228,
          df1 = 5,
          df2 = 66,
          eqbound = 0.35)

## -----------------------------------------------------------------------------
equ_anova(aovtest,
          eqbound = 0.35)

## -----------------------------------------------------------------------------
equ_anova(aovtest,
          eqbound = 0.35,
          MET = TRUE)

## ----fig.width=7, fig.height=6------------------------------------------------
plot_pes(Fstat = 34.70228,
         df1 = 5,
         df2 = 66)

## -----------------------------------------------------------------------------
# Simulate data with a small effect
set.seed(123)
groups <- factor(rep(1:3, each = 30))
y <- rnorm(90) + rep(c(0, 0.3, 0.3), each = 30)
small_aov <- aov(y ~ groups)

# Traditional ANOVA
knitr::kable(broom::tidy(small_aov),
            caption = "Traditional ANOVA Test (Small Effect)")

# Equivalence test
equ_anova(small_aov, eqbound = 0.15)

# Visualize
plot_pes(Fstat = 2.36, df1 = 2, df2 = 87)

## -----------------------------------------------------------------------------
power_eq_f(df1 = 2,         # Numerator df (groups - 1) 
           df2 = 60,      # Set to NULL to solve for sample size
           eqbound = 0.15)  # Equivalence bound)     

## -----------------------------------------------------------------------------
power_eq_f(df1 = 2,         # Numerator df (groups - 1)
           df2 = 60,        # Error df (N - groups)
           eqbound = 0.15)  # Equivalence bound

## -----------------------------------------------------------------------------
power_eq_f(df1 = 2,         # Numerator df (groups - 1)
           df2 = 60,        # Error df (N - groups)
           power = 0.8)     # Desired power

