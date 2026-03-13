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

## -----------------------------------------------------------------------------
boot_smd_calc(x = sleep$extra[sleep$group == 1],
              y = sleep$extra[sleep$group == 2],
         R = 199,
         paired = TRUE,
         boot_ci = "stud",
         bias_correction = TRUE,
         null.value = c(-0.5,0.5), 
         alternative = "equ")

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

## -----------------------------------------------------------------------------
set.seed(8484)
group1 <- rnorm(40, mean = 100, sd = 15)
group2 <- rnorm(40, mean = 110, sd = 15)

# Standard Cohen's d
smd_calc(x = group1, y = group2, bias_correction = FALSE)

# Robust version with 20% trimming (Algina et al., 2005)
smd_calc(x = group1, y = group2, bias_correction = FALSE, tr = 0.2)

# Robust Hedges' g with 10% trimming
smd_calc(x = group1, y = group2, bias_correction = TRUE, tr = 0.1)

## -----------------------------------------------------------------------------
boot_smd_calc(x = group1, y = group2, tr = 0.2, R = 999, boot_ci = "perc")

## -----------------------------------------------------------------------------
set.seed(7171)
n <- 50

# Clean data from known populations
x_clean <- rnorm(n, mean = 0, sd = 1)
y_clean <- rnorm(n, mean = 0.5, sd = 1)

# Contaminated version: replace 10% with outliers
n_contam <- ceiling(0.1 * n)
x_contam <- c(x_clean[1:(n - n_contam)], rnorm(n_contam, mean = 0, sd = 10))
y_contam <- c(y_clean[1:(n - n_contam)], rnorm(n_contam, mean = 0.5, sd = 10))

# Standard Cohen's d: attenuated by inflated SD
smd_calc(x = x_contam, y = y_contam, bias_correction = FALSE, smd_ci = "z")

# Trimmed version: closer to the true separation
smd_calc(x = x_contam, y = y_contam, bias_correction = FALSE, tr = 0.2, smd_ci = "z")

# For reference: Cohen's d on the clean data
smd_calc(x = x_clean, y = y_clean, bias_correction = FALSE, smd_ci = "z")

