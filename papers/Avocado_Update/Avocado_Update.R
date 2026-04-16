## ----setup, echo=FALSE--------------------------------------------------------
knitr::opts_chunk$set(
  warning = FALSE,
  message = FALSE,
  echo = TRUE,
  fig.pos = "H"
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
res1 = t_TOST(formula = extra ~ group, # difference in "group" of outcome "extra"
              data = sleep, # data
              eqb = .5, # equivalence bound
              smd_ci = "z") # sets type of SMD confidence interval calculation

## -----------------------------------------------------------------------------
print(res1)

## ----genericplot,fig.width=6, fig.height=5,fig.cap="Results of the independent-samples equivalence test visualized with the default plot method. The top panel displays the raw mean difference with its 90% confidence interval, and the bottom panel displays the standardized effect size (Hedges's g~av~) with its 90% confidence interval. In both panels, the dashed vertical lines indicate the equivalence bounds (raw: -0.5 to 0.5; standardized: approximately -0.3 to 0.3), and the point estimate is marked by a filled circle. A 90% confidence interval is used because it corresponds to the two one-sided tests (TOST) procedure at an alpha level of 0.05."----
plot(res1)

## ----cdplot,fig.width=6, fig.height=5,fig.cap="Consonance density plot for the independent-samples equivalence test. The top panel displays the standardized effect size (Hedges's g~av~) and the bottom panel displays the raw mean difference. In each panel, the consonance density curve represents the distribution of confidence intervals across all levels, with color-coded regions corresponding to the 68%, 90%, 95%, and 99.9% confidence intervals. The point estimate is marked by a filled circle along the x-axis. The dashed vertical lines indicate the equivalence bounds (standardized: -0.26 to 0.26; raw: -0.5 to 0.5). This visualization conveys the full range of effect sizes compatible with the data at varying confidence levels, rather than relying on a single confidence interval. Wider, less certain intervals (e.g., 99.9%) extend further from the point estimate, while narrower intervals (e.g., 68%) cluster near it."----
plot(res1, type = "cd")

## -----------------------------------------------------------------------------
res4 = t_TOST(x = bugs$LDHF, # simple sample, vector from which to estimate mean
              mu = 7.5, # set nil significance test to 7.5
              eqb = c(5.5,8.5) # set TOST bounds to 5.5 and 8.5
              )
res4

## -----------------------------------------------------------------------------
res_tsum = tsum_TOST(
  m1 = mean(bugs$LDHF, na.rm=TRUE), 
  sd1 = sd(bugs$LDHF, na.rm=TRUE),
  n1 = length(na.omit(bugs$LDHF)),
  hypothesis = "EQU", smd_ci = "t", 
  eqb = c(5.5, 8.5)
)

res_tsum

## -----------------------------------------------------------------------------
res_sh1 <- simple_htest(extra ~ group, 
                        data = sleep,
                        alternative = "equivalence", 
                        mu = 0.5)
res_sh1

## -----------------------------------------------------------------------------
describe_htest(res_sh1)

## ----shplot1, fig.width=6, fig.height=3.5, fig.cap = "Estimate plot for an equivalence test via simple_htest."----
plot_htest_est(res_sh1)

## -----------------------------------------------------------------------------
smd_calc(formula = extra ~ group,
         data = sleep,
         paired = TRUE,
         bias_correction = TRUE,
         rm_correction = TRUE,
         output = "htest")

## -----------------------------------------------------------------------------
smd_calc(formula = extra ~ group,
         data = sleep,
         paired = TRUE,
         bias_correction = TRUE,
         rm_correction = TRUE,
         alternative = "equivalence",
         null.value = c(-0.5,0.5),
         output = "htest")

## -----------------------------------------------------------------------------
boot_smd_calc(formula = extra ~ group,
              data = sleep,
              paired = TRUE,
              bias_correction = TRUE,
              rm_correction = TRUE,
              alternative = "equivalence",
              null.value = c(-0.5,0.5))

