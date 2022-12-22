## -----------------------------------------------------------------------------
library(TOSTER)
cor.test(mtcars$mpg,
           mtcars$qsec)

z_cor_test(mtcars$mpg,
           mtcars$qsec)

## -----------------------------------------------------------------------------
z_cor_test(mtcars$mpg,
           mtcars$qsec,
           method = "spear") # Don't need to spell full name
z_cor_test(mtcars$mpg,
           mtcars$qsec,
           method = "kendall")

## -----------------------------------------------------------------------------
z_cor_test(mtcars$mpg,
           mtcars$qsec,
           alternative = "e", # e for equivalence
           null = .4)

## -----------------------------------------------------------------------------
corsum_test(r = .121,
            n = 105,
            alternative = "e",
            null = .4)

## -----------------------------------------------------------------------------
set.seed(993)
boot_cor_test(mtcars$mpg,
           mtcars$qsec,
           alternative = "e",
           null = .4)

boot_cor_test(mtcars$mpg,
           mtcars$qsec,
           method = "spear",
           alternative = "e",
           null = .4)

boot_cor_test(mtcars$mpg,
           mtcars$qsec,
           method = "ken",
           alternative = "e",
           null = .4)

## -----------------------------------------------------------------------------
boot_cor_test(mtcars$mpg,
           mtcars$qsec,
           method = "win",
           alternative = "e",
           null = .4,
           tr = .1) # set trim

boot_cor_test(mtcars$mpg,
           mtcars$qsec,
           method = "bend",
           alternative = "e",
           null = .4,
           beta = .15) # bend argument

## -----------------------------------------------------------------------------
compare_cor(r1 = .8,
            df1 = 38,
            r2 = .2,
            df2 = 98)

## -----------------------------------------------------------------------------
compare_cor(r1 = .8,
            df1 = 38,
            r2 = .2,
            df2 = 98)

## -----------------------------------------------------------------------------
compare_cor(r1 = .8,
            df1 = 38,
            r2 = .2,
            df2 = 98,
            null = .2,
            method = "f",
            alternative = "e")

## -----------------------------------------------------------------------------
set.seed(8922)
x1 = rnorm(40)
y1 = rnorm(40)

x2 = rnorm(100)
y2 = rnorm(100)

boot_compare_cor(
  x1 = x1,
  x2 = x2,
  y1 = y1,
  y2 = y2,
  null = .2,
  alternative = "e",
  method = "win"
)

