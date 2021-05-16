## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(TOSTER)
library(ggplot2)
library(ggdist)


## ----fig.width=6, fig.height=6------------------------------------------------
plot_smd(d = .43,
         df = 58,
         lambda = 1.66,
         smd_label = "Cohen's d"
         )

