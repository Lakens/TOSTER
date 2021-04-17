## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(TOSTER)
library(ggplot2)
library(ggdist)


## ----echo=FALSE, warning = FALSE, fig.show='hold'-----------------------------

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
  label = expression('H'[0])) +
  geom_text(aes(
    y = 1,
    x = 1.5,
    vjust = 1.5,
    hjust = "middle"
  ),
  #alignment = "center",
  label = expression('H'["A"])) +
  geom_text(aes(
    y = 1,
    x = -1.5,
    vjust = 1.5,
    hjust = "middle"
  ),
  #alignment = "center",
  label = expression('H'["A"])) +
theme_tidybayes() + 
  scale_y_continuous(limits = c(0,1.75)) +
  scale_x_continuous(limits = c(-2,2)) +
  labs(x = "", y = "",title="Minimal Effect Test") +
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
  label = expression('H'["A"])) +
  geom_text(aes(
    y = 1,
    x = 1.5,
    vjust = 1.5,
    hjust = "middle"
  ),
  #alignment = "center",
  label = expression('H'["0"])) +
  geom_text(aes(
    y = 1,
    x = -1.5,
    vjust = 1.5,
    hjust = "middle"
  ),
  #alignment = "center",
  label = expression('H'["0"])) +
theme_tidybayes() + 
  scale_y_continuous(limits = c(0,1.75)) +
  scale_x_continuous(limits = c(-2,2)) +
  labs(x = "", y = "",title="Equivalence Test") +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank() 
  )

## -----------------------------------------------------------------------------
data('sleep')
library(jmv)
data('bugs')

## -----------------------------------------------------------------------------
head(sleep)

## -----------------------------------------------------------------------------
res1 = t.TOST(formula = extra ~ group,
              data = sleep,
              low_eqbound = -.5,
              high_eqbound = .5)

res1a = t.TOST(x = subset(sleep,group==1)$extra,
               y = subset(sleep,group==2)$extra,
               low_eqbound = -.5,
               high_eqbound = .5)

## -----------------------------------------------------------------------------
print(res1)

## ----fig.width=6, fig.height=6------------------------------------------------
plot(res1)

## -----------------------------------------------------------------------------
res2 = t.TOST(formula = extra ~ group,
              data = sleep,
              paired = TRUE,
              low_eqbound = -.5,
              high_eqbound = .5)
res2

## -----------------------------------------------------------------------------
res3 = t.TOST(x = bugs$LDHF,
              y = bugs$LDLF,
              paired = TRUE,
              low_eqbound = -1,
              high_eqbound = 1)
res3

## -----------------------------------------------------------------------------
res3a = t.TOST(x = bugs$LDHF,
               y = bugs$LDLF,
               paired = TRUE,
               hypothesis = "MET",
               low_eqbound = -1,
               high_eqbound = 1)
res3a

## ----fig.width=6, fig.height=6------------------------------------------------
plot(res3a)

## -----------------------------------------------------------------------------
res4 = t.TOST(x = bugs$LDHF,
              hypothesis = "EQU",
              low_eqbound = 5.5,
              high_eqbound = 8.5)
res4

## ----fig.width=6, fig.height=6------------------------------------------------
plot(res4)

