## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(TOSTER)
library(ggplot2)
library(ggdist)


## -----------------------------------------------------------------------------
smd_calc(formula = extra ~ group,
         data = sleep,
         paired = TRUE,
         smd_ci = "nct",
         bias_correction = F)

## ----fig.width=6, fig.height=6------------------------------------------------
plot_smd(d = .43,
         df = 58,
         lambda = 1.66,
         smd_label = "Cohen's d"
         )

## -----------------------------------------------------------------------------
compare_smd(smd1 = 0.95,
            n1 = 25,
            smd2 = 0.23,
            n2 = 50,
            paired = TRUE)

## -----------------------------------------------------------------------------
set.seed(4522)
diff_study1 = rnorm(25,.95)
diff_study2 = rnorm(50)
boot_test = boot_compare_smd(x1 = diff_study1,
                             x2 = diff_study2,
                             paired = TRUE)

boot_test

# Table of bootstrapped CIs
knitr::kable(boot_test$df_ci, digits = 4)

## -----------------------------------------------------------------------------
library(ggplot2)

list_res = boot_test$boot_res
df1 = data.frame(study = c(rep("original",length(list_res$smd1)),
                           rep("replication",length(list_res$smd2))),
                 smd = c(list_res$smd1,list_res$smd2))
ggplot(df1,
            aes(fill = study, color =smd, x = smd))+ 
 geom_histogram(aes(y=..density..), alpha=0.5, 
                position="identity")+
 geom_density(alpha=.2) +
  labs(y = "", x = "SMD (bootstrapped estimates)") +
  theme_classic()

df2 = data.frame(diff = list_res$d_diff)
ggplot(df2,
            aes(x = diff))+ 
 geom_histogram(aes(y=..density..), alpha=0.5, 
                position="identity")+
 geom_density(alpha=.2) +
  labs(y = "", x = "Difference in SMDs (bootstrapped estimates)") +
  theme_classic()

