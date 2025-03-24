## -----------------------------------------------------------------------------
library(TOSTER)
# Base R correlation test
cor.test(mtcars$mpg, mtcars$qsec)

# TOSTER's z-transformed correlation test
z_cor_test(mtcars$mpg, mtcars$qsec)

## -----------------------------------------------------------------------------
# Spearman correlation
z_cor_test(mtcars$mpg,
           mtcars$qsec,
           method = "spear") # Short form accepted; "spearman" also works

# Kendall correlation
z_cor_test(mtcars$mpg,
           mtcars$qsec,
           method = "kendall")

## -----------------------------------------------------------------------------
# Equivalence test with null boundary of 0.4
z_cor_test(mtcars$mpg,
           mtcars$qsec,
           alternative = "e", # e for equivalence
           null = .4)

## -----------------------------------------------------------------------------
# Testing a correlation of 0.121 from a sample of 105 paired observations
corsum_test(r = .121,
            n = 105,
            alternative = "e",
            null = .4)

## -----------------------------------------------------------------------------
set.seed(993) # Setting seed for reproducibility
boot_cor_test(mtcars$mpg,
           mtcars$qsec,
           alternative = "e",
           null = .4)

# Bootstrapped Spearman correlation
boot_cor_test(mtcars$mpg,
           mtcars$qsec,
           method = "spear",
           alternative = "e",
           null = .4)

# Bootstrapped Kendall correlation
boot_cor_test(mtcars$mpg,
           mtcars$qsec,
           method = "ken", # Short form accepted
           alternative = "e",
           null = .4)

## -----------------------------------------------------------------------------
# Winsorized correlation with 10% trimming
boot_cor_test(mtcars$mpg,
           mtcars$qsec,
           method = "win",
           alternative = "e",
           null = .4,
           tr = .1) # Set trim amount (default is 0.2)

# Percentage bend correlation
boot_cor_test(mtcars$mpg,
           mtcars$qsec,
           method = "bend",
           alternative = "e",
           null = .4,
           beta = .15) # Beta parameter controlling resistance to outliers

## -----------------------------------------------------------------------------
# Comparing correlation r1=0.8 from n=40 with r2=0.2 from n=100
compare_cor(r1 = .8,
            df1 = 38,  # df = n-2
            r2 = .2,
            df2 = 98)  # df = n-2

## -----------------------------------------------------------------------------
# Testing equivalence using Fisher's method
compare_cor(r1 = .8,
            df1 = 38,
            r2 = .2,
            df2 = 98,
            null = .2,
            method = "f", # Fisher (can also use "fisher")
            alternative = "e") # Equivalence test

## -----------------------------------------------------------------------------
set.seed(8922) # Setting seed for reproducibility
# Generating example data
x1 = rnorm(40)
y1 = rnorm(40)

x2 = rnorm(100)
y2 = rnorm(100)

# Bootstrap comparison with winsorized correlation
boot_compare_cor(
  x1 = x1,
  x2 = x2,
  y1 = y1,
  y2 = y2,
  null = .2,
  alternative = "e", # Equivalence test
  method = "win" # Winsorized correlation
)

## ----eval=FALSE---------------------------------------------------------------
#  # Customizing the bootstrap procedure
#  boot_cor_test(
#    x = mtcars$mpg,
#    y = mtcars$qsec,
#    method = "pearson",
#    R = 2000,  # Increasing number of bootstrap samples
#    alpha = 0.01,  # Using 99% confidence interval
#    alternative = "t"  # Two-sided test
#  )

## ----eval=FALSE---------------------------------------------------------------
#  # Example with missing data
#  x_with_na <- c(mtcars$mpg, NA, NA)
#  y_with_na <- c(mtcars$qsec, 10, NA)
#  
#  # Default behavior handles NAs with pairwise deletion
#  z_cor_test(x_with_na, y_with_na)

