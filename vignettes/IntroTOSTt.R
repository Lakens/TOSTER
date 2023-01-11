## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(TOSTER)
library(ggplot2)
library(ggdist)


## ----tostplots,echo=FALSE, message = FALSE, warning = FALSE, fig.show='hold'----

ggplot() +
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

ggplot() +
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

## -----------------------------------------------------------------------------
data('sleep')
data('iris')

## -----------------------------------------------------------------------------
head(sleep)

## -----------------------------------------------------------------------------
res1 = t_TOST(formula = extra ~ group,
              data = sleep,
              eqb = .5,
              smd_ci = "goulet")

res1a = t_TOST(x = subset(sleep,group==1)$extra,
               y = subset(sleep,group==2)$extra,
               eqb = .5)

## -----------------------------------------------------------------------------
# Simple htest

res1b = simple_htest(formula = extra ~ group,
                     data = sleep,
                     mu = .5, # set equivalence bound
                     alternative = "e")



## -----------------------------------------------------------------------------

# t_TOST
print(res1)

# htest

print(res1b)

## ----fig.width=6, fig.height=6------------------------------------------------
plot(res1, type = "cd")

## ----fig.width=6, fig.height=6, eval=FALSE------------------------------------
#  # Set to shade only the 90% and 95% CI areas
#  plot(res1, type = "cd",
#       ci_shades = c(.9,.95))

## ----fig.width=6, fig.height=6------------------------------------------------
plot(res1, type = "c",
     ci_lines =  c(.9,.95))

## -----------------------------------------------------------------------------
res2 = t_TOST(formula = extra ~ group,
              data = sleep,
              paired = TRUE,
              eqb = .5)
res2

res2b = simple_htest(
  formula = extra ~ group,
  data = sleep,
  paired = TRUE,
  mu = .5,
  alternative = "e")
res2b

## -----------------------------------------------------------------------------
res3 = t_TOST(x = iris$Sepal.Length,
              y = iris$Sepal.Width,
              paired = TRUE,
              eqb = 1)
res3

res3a = simple_htest(
  x = iris$Sepal.Length,
  y = iris$Sepal.Width,
  paired = TRUE,
  mu = 1,
  alternative = "e"
)
res3a

## -----------------------------------------------------------------------------
res_met = t_TOST(x = iris$Sepal.Length,
              y = iris$Sepal.Width,
               paired = TRUE,
               hypothesis = "MET",
               eqb = 1,
              smd_ci = "goulet")
res_met

res_metb = simple_htest(x = iris$Sepal.Length,
                       y = iris$Sepal.Width,
                       paired = TRUE,
                       mu = 1,
                       alternative = "minimal.effect")
res_metb

## -----------------------------------------------------------------------------
res4 = t_TOST(x = iris$Sepal.Length,
              hypothesis = "EQU",
              eqb = c(5.5,8.5),
              smd_ci = "goulet")
res4

## -----------------------------------------------------------------------------
res_tsum = tsum_TOST(
  m1 = mean(iris$Sepal.Length, na.rm=TRUE),
  sd1 = sd(iris$Sepal.Length, na.rm=TRUE),
  n1 = length(na.omit(iris$Sepal.Length)),
  hypothesis = "EQU",
  eqb = c(5.5,8.5)
)

res_tsum

## ----fig.width=6, fig.height=6------------------------------------------------
plot(res_tsum)

## -----------------------------------------------------------------------------
power_t_TOST(n = NULL,
  delta = 1,
  sd = 2.5,
  eqb = 2.5,
  alpha = .025,
  power = .95,
  type = "two.sample")

