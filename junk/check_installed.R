library(TOSTER)
data(sleep)
x <- sleep$extra[sleep$group == 2]
y <- sleep$extra[sleep$group == 1]
res <- suppressMessages(ses_calc(x=x, y=y, paired=TRUE, ses="cstat", se_method="agresti"))
cat("Installed TOSTER agresti paired cstat:", as.numeric(res$estimate), "\n")
cat("Label:", names(res$estimate), "\n")
