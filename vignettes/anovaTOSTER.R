## ----warning=FALSE, message=FALSE---------------------------------------------
library(TOSTER)
# Get Data
data("InsectSprays")
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
powerTOST_f(df1 = 2, 
            df2 = 60,
            eqbound = .15)

