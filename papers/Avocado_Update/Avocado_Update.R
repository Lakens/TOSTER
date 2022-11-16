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

## ----hypplot, fig.width=6, fig.height=2.75, echo=FALSE, message = FALSE, warning = FALSE, fig.show='hold', fig.cap = "Type of Hypothesis"----

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
head(sleep,2)

## -----------------------------------------------------------------------------
# Formula Interface
res1 = t_TOST(formula = extra ~ group, data = sleep, 
              eqb = .5, smd_ci = "t")
# x & y Interface
res1a = t_TOST(x = subset(sleep,group==1)$extra,
               y = subset(sleep,group==2)$extra, eqb =.5)

## -----------------------------------------------------------------------------
print(res1)

## ----cdplot,fig.width=6, fig.height=5,fig.cap="Example of consonance density plot."----
plot(res1, type = "cd")

## ----shadeplot,fig.width=6, fig.height=5, fig.cap = "Demonstrating the shading in plot method."----
plot(res1, type = "cd",
     ci_shades = c(.9,.95))

## ----conplot,fig.width=6, fig.height=5, fig.cap = "Example of consonance plot."----
plot(res1, type = "c",
     ci_lines =  c(.9,.95))

## -----------------------------------------------------------------------------
res2 = t_TOST(formula = extra ~ group,
              data = sleep,
              paired = TRUE,
              eqb = .5)
res2

## -----------------------------------------------------------------------------
res3 = t_TOST(x = bugs$LDHF,
              y = bugs$LDLF,
              paired = TRUE,
              eqb = 1)
res3

## -----------------------------------------------------------------------------
res3a = t_TOST(x = bugs$LDHF,
               y = bugs$LDLF,
               paired = TRUE,
               hypothesis = "MET",
               eqb = 1)
res3a

## -----------------------------------------------------------------------------
res4 = t_TOST(x = bugs$LDHF,
              hypothesis = "EQU",
              mu = 7.5,
              eqb = c(5.5,8.5))
res4

## -----------------------------------------------------------------------------
res_tsum = tsum_TOST(
  m1 = mean(bugs$LDHF, na.rm=TRUE), sd1 = sd(bugs$LDHF, na.rm=TRUE),
  n1 = length(na.omit(bugs$LDHF)),
  hypothesis = "EQU", smd_ci = "t", eqb = c(5.5, 8.5)
)

res_tsum

## -----------------------------------------------------------------------------
test1 = wilcox_TOST(formula = extra ~ group,
                      data = sleep,
                      paired = FALSE,
                      eqb = .5)
print(test1)

## -----------------------------------------------------------------------------
set.seed(891111)
test1 = boot_t_TOST(formula = extra ~ group,
                    data = sleep,
                    paired = TRUE,
                    eqb = .5,
                    R = 999)


print(test1)

## -----------------------------------------------------------------------------
x = 7; y = 10.5
log(y) - log(x)
log(y/x)
exp(log(y) - log(x))
y/x

## ---- error=FALSE-------------------------------------------------------------
log_TOST(mpg ~ am, data = mtcars)

## ---- error=FALSE-------------------------------------------------------------
boot_log_TOST(mpg ~ am, data = mtcars, R=999)

## ----warning=FALSE, message=FALSE---------------------------------------------
data("InsectSprays")
aovtest = aov(count ~ spray, data = InsectSprays)
anova(aovtest)


## -----------------------------------------------------------------------------
equ_ftest(Fstat = 34.70228,  df1 = 5, df2 = 66,  eqb = 0.35)

## -----------------------------------------------------------------------------
# Example using a purely within-subjects design 
# (Maxwell & Delaney, 2004, Chapter 12, Table 12.5, p. 578):
library(afex)
data(md_12.1)
aovtest2 = aov_ez("id", "rt", md_12.1, within = c("angle", "noise"), 
       anova_table=list(correction = "none", es = "none"))
equ_anova(aovtest2,
          eqb = 0.35)

## -----------------------------------------------------------------------------
compare_smd(smd1 = 0.95,
            n1 = 25,
            smd2 = 0.23,
            n2 = 50,
            paired = TRUE)

## -----------------------------------------------------------------------------
compare_smd(smd1 = 0.95, n1 = 25, smd2 = 0.23,n2 = 50,
            paired = TRUE, TOST = TRUE, null = .25)

## -----------------------------------------------------------------------------
set.seed(4522)
boot_test = boot_compare_smd(x1 = rnorm(25,.95), x2 = rnorm(50), 
                             paired = TRUE, alpha = .1)
boot_test

