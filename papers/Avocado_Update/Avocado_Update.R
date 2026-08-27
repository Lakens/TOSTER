## ----setup, echo=FALSE--------------------------------------------------------
knitr::opts_chunk$set(
  warning = FALSE,
  message = FALSE,
  echo = TRUE,
  fig.pos = "H",
  fig.width = 6,
  fig.height = 4
)
knitr::knit_hooks$set(purl = knitr::hook_purl)
library(TOSTER)
library(ggplot2)
library(ggdist)
library(patchwork)

## ----hypplot, fig.width=6, fig.height=2.75, echo=FALSE, message = FALSE, warning = FALSE, fig.show='hold', fig.cap = "Types of hypotheses tests supported by two one-sided tests (TOST) procedures."----

p1 = ggplot() +
  geom_vline(aes(xintercept = -.5),
             linetype = "dashed") +
  geom_vline(aes(xintercept = .5),
             linetype = "dashed") +
  geom_text(aes(
    y = 1,
    x = -0.5,
    vjust = -.9,
    hjust = "middle"
  ),
  angle = 90,
  label = 'Lower Bound') +
  geom_text(aes(
    y = 1,
    x = 0.5,
    vjust = 1.5,
    hjust = "middle"
  ),
  angle = 90,
  label = 'Upper Bound') +
  geom_text(aes(
    y = 1,
    x = 0,
    vjust = 1.5,
    hjust = "middle"
  ),
  #alignment = "center",
  label = "H0"
  ) +
  geom_text(aes(
    y = 1,
    x = 1.5,
    vjust = 1.5,
    hjust = "middle"
  ),
  #alignment = "center",
  label = "H1"
  ) +
  geom_text(aes(
    y = 1,
    x = -1.5,
    vjust = 1.5,
    hjust = "middle"
  ),
  #alignment = "center",
  label = "H1"
  ) +
theme_tidybayes() +
  scale_y_continuous(limits = c(0,1.75)) +
  scale_x_continuous(limits = c(-2,2)) +
  labs(x = "", y = "",
       title="Minimal Effect Test",
       caption = "H1 = Alternative Hypothesis \n H0 = Null Hypothesis") +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

p2 = ggplot() +
  geom_vline(aes(xintercept = -.5),
             linetype = "dashed") +
  geom_vline(aes(xintercept = .5),
             linetype = "dashed") +
  geom_text(aes(
    y = 1,
    x = -0.5,
    vjust = -.9,
    hjust = "middle"
  ),
  angle = 90,
  label = 'Lower Bound') +
  geom_text(aes(
    y = 1,
    x = 0.5,
    vjust = 1.5,
    hjust = "middle"
  ),
  angle = 90,
  label = 'Upper Bound') +
  geom_text(aes(
    y = 1,
    x = 0,
    vjust = 1.5,
    hjust = "middle"
  ),
  #alignment = "center",
  label = "H1"
  ) +
  geom_text(aes(
    y = 1,
    x = 1.5,
    vjust = 1.5,
    hjust = "middle"
  ),
  #alignment = "center",
  label = "H0"
  ) +
  geom_text(aes(
    y = 1,
    x = -1.5,
    vjust = 1.5,
    hjust = "middle"
  ),
  #alignment = "center",
  label = "H0"
  ) +
theme_tidybayes() +
  scale_y_continuous(limits = c(0,1.75)) +
  scale_x_continuous(limits = c(-2,2)) +
  labs(x = "",
       y = "",
       title="Equivalence Test",
       caption = "H1 = Alternative Hypothesis \n H0 = Null Hypothesis") +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

p1 + p2 + plot_annotation(tag_levels = 'A')

## -----------------------------------------------------------------------------
data('sleep')
library(jmv)
data('bugs')

## -----------------------------------------------------------------------------
res1 = t_TOST(
  # difference in "group" of outcome "extra"
  formula = extra ~ group,
  # data
  data = sleep,
  # equivalence bound
  eqb = .5,
  # sets type of SMD confidence interval
  smd_ci = "z"
)

## -----------------------------------------------------------------------------
print(res1)

## ----results = "asis"---------------------------------------------------------
cat(describe(res1))

## ----genericplot,fig.width=6, fig.height=5,fig.cap="Results of the independent-samples equivalence test visualized with the default plot method. The top panel displays the raw mean difference with its 90% confidence interval, and the bottom panel displays the standardized effect size (Hedges's g~av~) with its 90% confidence interval. In both panels, the dashed vertical lines indicate the equivalence bounds (raw: -0.5 to 0.5; standardized: approximately -0.3 to 0.3), and the point estimate is marked by a filled circle. A 90% confidence interval is used because it corresponds to the two one-sided tests (TOST) procedure at an alpha level of 0.05."----
plot(res1)

## ----cdplot,fig.width=6, fig.height=5,fig.cap="Consonance density plot for the independent-samples equivalence test. The top panel displays the standardized effect size (Hedges's g~av~) and the bottom panel displays the raw mean difference. In each panel, the consonance density curve represents the distribution of confidence intervals across all levels, with color-coded regions corresponding to the 68%, 90%, 95%, and 99.9% confidence intervals. The point estimate is marked by a filled circle along the x-axis. The dashed vertical lines indicate the equivalence bounds (standardized: -0.26 to 0.26; raw: -0.5 to 0.5). This visualization conveys the full range of effect sizes compatible with the data at varying confidence levels, rather than relying on a single confidence interval. Wider, less certain intervals (e.g., 99.9%) extend further from the point estimate, while narrower intervals (e.g., 68%) cluster near it."----
plot(res1, type = "cd")

## -----------------------------------------------------------------------------
res4 = t_TOST(
  # single sample, vector from which to estimate mean
  x = bugs$LDHF,
  # set nil significance test to 7.5
  mu = 7.5,
  # set TOST bounds to 5.5 and 8.5
  eqb = c(5.5, 8.5)
)
res4

## -----------------------------------------------------------------------------
res_tsum = tsum_TOST(
  # sample mean
  m1 = mean(bugs$LDHF, na.rm = TRUE),
  # sample standard deviation
  sd1 = sd(bugs$LDHF, na.rm = TRUE),
  # sample size
  n1 = length(na.omit(bugs$LDHF)),
  # equivalence rather than minimal effects test
  hypothesis = "EQU",
  # sets type of SMD confidence interval
  smd_ci = "t",
  # equivalence bounds
  eqb = c(5.5, 8.5)
)

res_tsum

## -----------------------------------------------------------------------------
res_sh1 = simple_htest(
  extra ~ group,
  data = sleep,
  paired = TRUE,
  # equivalence hypothesis
  alternative = "equivalence",
  # symmetric bounds of -0.5 and 0.5
  mu = 0.5
)
res_sh1

## -----------------------------------------------------------------------------
res_sh2 = simple_htest(
  extra ~ group,
  data = sleep,
  paired = TRUE,
  # non-inferiority rather than equivalence
  alternative = "greater",
  # margin of non-inferiority
  mu = -0.5
)
res_sh2

## ----results = "asis"---------------------------------------------------------
cat(describe_htest(res_sh1))

## ----shplot1, fig.width=6, fig.height=3.5, fig.cap = "Estimate plot for an equivalence test via simple_htest."----
plot_htest_est(res_sh1)

## -----------------------------------------------------------------------------
smd_calc(
  formula = extra ~ group,
  data = sleep,
  paired = TRUE,
  # apply Hedges' correction
  bias_correction = TRUE,
  # apply the repeated measures correction
  rm_correction = TRUE,
  # return an htest object rather than a data frame
  output = "htest"
)

## -----------------------------------------------------------------------------
smd_calc(
  formula = extra ~ group,
  data = sleep,
  paired = TRUE,
  # apply Hedges' correction
  bias_correction = TRUE,
  # apply the repeated measures correction
  rm_correction = TRUE,
  # equivalence hypothesis
  alternative = "equivalence",
  # equivalence bounds on the SMD
  null.value = c(-0.5, 0.5),
  # return an htest object rather than a data frame
  output = "htest"
)

## -----------------------------------------------------------------------------
boot_smd_calc(
  formula = extra ~ group,
  data = sleep,
  paired = TRUE,
  # apply Hedges' correction
  bias_correction = TRUE,
  # apply the repeated measures correction
  rm_correction = TRUE,
  # equivalence hypothesis
  alternative = "equivalence",
  # equivalence bounds on the SMD
  null.value = c(-0.5, 0.5)
)

## -----------------------------------------------------------------------------
z_cor_test(
  x = bugs$LDLF,
  y = bugs$LDHF,
  # type of correlation coefficient
  method = "pearson",
  # two-sided nil hypothesis test
  alternative = "t"
)

## -----------------------------------------------------------------------------
res_boot_cor = boot_cor_test(
  x = bugs$LDLF,
  y = bugs$LDHF,
  # type of correlation coefficient
  method = "pearson",
  # equivalence hypothesis
  alternative = "equivalence",
  # equivalence bounds of -0.4 and 0.4
  null = 0.4
)

print(res_boot_cor)

## ----results = "asis"---------------------------------------------------------
cat(describe_htest(res_boot_cor))

## -----------------------------------------------------------------------------
# Significance test
compare_smd(
  # SMD and sample size of the original study
  smd1 = 0.95,
  n1 = 25,
  # SMD and sample size of the replication
  smd2 = 0.23,
  n2 = 50,
  # both estimates come from paired designs
  paired = TRUE,
  # nil hypothesis test
  null = 0,
  # significance hypothesis
  alternative = "two.sided"
)

# Equivalence test
compare_smd(
  # SMD and sample size of the original study
  smd1 = 0.95,
  n1 = 25,
  # SMD and sample size of the replication
  smd2 = 0.23,
  n2 = 50,
  # both estimates come from paired designs
  paired = TRUE,
  # equivalence bounds on the difference in SMDs
  null = .25,
  # equivalence hypothesis
  alternative = "equivalence"
)

## -----------------------------------------------------------------------------
compare_cor(
  # correlation and degrees of freedom of the first study
  r1 = 0.45,
  df1 = 48,
  # correlation and degrees of freedom of the second study
  r2 = 0.25,
  df2 = 78,
  # equivalence hypothesis
  alternative = "equivalence",
  # equivalence bounds on the difference in correlations
  null = 0.25
)

## -----------------------------------------------------------------------------
test1 = wilcox_TOST(
  formula = extra ~ group,
  data = sleep,
  paired = TRUE,
  # equivalence bounds on the pseudo-median
  eqb = .5
)
print(test1)

## -----------------------------------------------------------------------------
# symmetry plot for sleep paired differences
d_sleep <- sleep$extra[sleep$group == 2] - sleep$extra[sleep$group == 1]
m <- median(d_sleep)
sorted <- sort(d_sleep)
upper <- sorted[sorted > m] - m
lower <- m - sorted[sorted < m]
k <- min(length(upper), length(lower))
sp_df <- data.frame(lower = sort(lower)[1:k], upper = sort(upper)[1:k])
mx <- max(sp_df$lower, sp_df$upper)

ggplot(sp_df, aes(x = lower, y = upper)) +
  geom_point(size = 2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  coord_equal(xlim = c(0, mx), ylim = c(0, mx)) +
  labs(x = "Distance below median",
       y = "Distance above median") +
  theme_minimal()

## -----------------------------------------------------------------------------
bm_test = brunner_munzel(
  formula = extra ~ group,
  data = sleep,
  # equivalence hypothesis
  alternative = "equivalence",
  # equivalence bounds on the probability scale
  mu = c(0.3, 0.7)
)
print(bm_test)

## -----------------------------------------------------------------------------
set.seed(4522)
boot_t_test(
  formula = extra ~ group,
  data = sleep,
  paired = TRUE,
  # equivalence hypothesis
  alternative = "equivalence",
  # equivalence bounds
  mu = c(-0.5, 0.5),
  # number of bootstrap resamples
  R = 999
)

## -----------------------------------------------------------------------------
set.seed(891111)
test1 = boot_t_TOST(
  formula = extra ~ group,
  data = sleep,
  paired = TRUE,
  # equivalence bounds
  eqb = .5,
  # number of bootstrap resamples
  R = 999
)
print(test1)

## -----------------------------------------------------------------------------
set.seed(8812)
perm_t_test(
  formula = extra ~ group,
  data = sleep,
  paired = TRUE,
  # equivalence hypothesis
  alternative = "equivalence",
  # equivalence bounds
  mu = c(-0.5, 0.5),
  # number of random permutations
  R = 999
)

## ----error=FALSE--------------------------------------------------------------
log_TOST(mpg ~ am, data = mtcars)

## ----error=FALSE--------------------------------------------------------------
boot_log_TOST(mpg ~ am, data = mtcars, R=999)

## ----warning=FALSE, message=FALSE---------------------------------------------
data("InsectSprays")
aovtest = aov(count ~ spray, data = InsectSprays)
anova(aovtest)


## -----------------------------------------------------------------------------
equ_ftest(
  # F statistic and degrees of freedom from the ANOVA table
  Fstat = 34.70228,
  df1 = 5,
  df2 = 66,
  # equivalence bound on partial eta-squared
  eqb = 0.35
)

## -----------------------------------------------------------------------------
# Example using a purely within-subjects design
# (Maxwell & Delaney, 2004, Chapter 12, Table 12.5, p. 578):
library(afex)
data(md_12.1)
aovtest2 = aov_ez(
  # participant identifier
  "id",
  # dependent variable
  "rt",
  # data
  md_12.1,
  # within-subjects factors
  within = c("angle", "noise"),
  anova_table = list(correction = "none", es = "none")
)
equ_anova(aovtest2, eqb = 0.35)

## -----------------------------------------------------------------------------
power_t_TOST(
  # true mean difference
  delta = 0,
  # assumed standard deviation
  sd = 1,
  # equivalence bounds
  eqb = 0.5,
  alpha = 0.05,
  # leaving n unspecified solves for the n giving 80% power
  power = 0.8,
  type = "two.sample"
)

## -----------------------------------------------------------------------------
power_z_cor(
  # true correlation
  rho = 0,
  # leaving n unspecified solves for the n giving 80% power
  power = 0.8,
  # equivalence bounds of -0.3 and 0.3
  null = 0.3,
  alpha = 0.05,
  # equivalence hypothesis
  alternative = "equivalence"
)

