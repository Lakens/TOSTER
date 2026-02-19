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
              eqb = .5,  # equivalence bounds of ±0.5 hours
              smd_ci = "t")  # t-distribution for SMD confidence intervals

# Alternative syntax with separate vectors
res1a = t_TOST(x = subset(sleep, group==1)$extra,
               y = subset(sleep, group==2)$extra,
               eqb = .5)

## -----------------------------------------------------------------------------
# Simple htest approach
res1b = simple_htest(formula = extra ~ group,
                     data = sleep,
                     mu = .5,  # equivalence bound
                     alternative = "e")  # "e" for equivalence

## -----------------------------------------------------------------------------
# Comprehensive t_TOST output
print(res1)

# Concise htest output
print(res1b)

## ----fig.width=6, fig.height=6------------------------------------------------
plot(res1, type = "simple", layout = "combined")

## ----fig.width=6, fig.height=6, eval=TRUE-------------------------------------
# Shade the 90% and 95% CI areas
plot(res1, type = "cd",
     ci_shades = c(.9, .95))

## ----fig.width=6, fig.height=6------------------------------------------------
plot(res1, type = "c",
     ci_lines = c(.9, .95))

## ----fig.width=6, fig.height=6------------------------------------------------
plot(res1, type = "tnull")

## ----eval = FALSE-------------------------------------------------------------
# describe(res1)
# describe_htest(res1b)

## -----------------------------------------------------------------------------
# For paired tests, use separate vectors rather than formula notation
res2 = t_TOST(x = sleep$extra[sleep$group == 1],
              y = sleep$extra[sleep$group == 2],
              paired = TRUE,  # specify paired analysis
              eqb = .5)
res2

res2b = simple_htest(
  x = sleep$extra[sleep$group == 1],
  y = sleep$extra[sleep$group == 2],
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
              hypothesis = "MET",  # Change to minimal effect test
              eqb = 1,
              smd_ci = "t")
res_met

res_metb = simple_htest(x = iris$Sepal.Length,
                       y = iris$Sepal.Width,
                       paired = TRUE,
                       mu = 1,
                       alternative = "minimal.effect")
res_metb

## ----eval = FALSE-------------------------------------------------------------
# describe(res_met)
# describe_htest(res_metb)

## ----error=TRUE---------------------------------------------------------------
try({
set.seed(221)
dat1 = rnorm(30)
dat2 = rnorm(30)

test = t_TOST(
  x = dat1,
  y = dat2,
  eqbound_type = "SMD",
  eqb = .2,
  hypothesis = "MET"
)

test
})

## -----------------------------------------------------------------------------
res4 = t_TOST(x = iris$Sepal.Length,
              hypothesis = "EQU",
              eqb = c(5.5, 8.5),  # lower and upper bounds
              smd_ci = "t")
res4

## -----------------------------------------------------------------------------
res_tsum = tsum_TOST(
  m1 = mean(iris$Sepal.Length, na.rm=TRUE),  # sample mean
  sd1 = sd(iris$Sepal.Length, na.rm=TRUE),   # sample standard deviation
  n1 = length(na.omit(iris$Sepal.Length)),  # sample size
  hypothesis = "EQU",
  eqb = c(5.5, 8.5)
)

res_tsum

## ----fig.width=6, fig.height=6------------------------------------------------
plot(res_tsum)

## -----------------------------------------------------------------------------
describe(res_tsum)

## -----------------------------------------------------------------------------
power_t_TOST(n = NULL,
  delta = 1,          # assumed true difference
  sd = 2.5,           # assumed standard deviation
  eqb = 2.5,          # equivalence bounds
  alpha = .025,       # significance level
  power = .95,        # desired power
  type = "two.sample") # test type

