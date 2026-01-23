## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(TOSTER)
library(ggplot2)
library(ggdist)


## -----------------------------------------------------------------------------
# For paired tests, use separate vectors
smd_calc(x = sleep$extra[sleep$group == 1],
         y = sleep$extra[sleep$group == 2],
         paired = TRUE,
         smd_ci = "nct",
         bias_correction = F)

# Setting bootstrap replications low to
## reduce compiling time of vignette
boot_smd_calc(x = sleep$extra[sleep$group == 1],
              y = sleep$extra[sleep$group == 2],
         R = 199,
         paired = TRUE,
         boot_ci = "stud",
         bias_correction = F)

## ----fig.width=6, fig.height=6------------------------------------------------
plot_smd(d = .43,
         df = 58,
         sigma = .33,
         smd_label = "Cohen's d",
         smd_ci = "z"
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

