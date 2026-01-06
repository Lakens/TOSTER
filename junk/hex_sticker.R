
# install.packages("hexSticker") # if needed
library(ggplot2)
library(hexSticker)

# Parameters for equivalence bounds (adjust to your SE and SESOI)


set.seed(56843)
res1 = t_TOST(rnorm(300),
              rnorm(300),
              eqb = .5,  # equivalence bounds of ±0.5 hours
              smd_ci = "t")  # t-distribution for SMD confidence intervals

plotted_null =  plot(res1, type = "tnull") +
  labs(caption = "")
# Create the hex sticker
sticker(
  plotted_null,
  package    = "TOSTER",
  p_size     = 30,            # size of 'TOSTER' label
  p_color    = "#2D2A26",     # label color
  s_x        = 1.0,           # motif x-position
  s_y        = 0.75,          # motif y-position
  s_width    = 1.4,           # motif width scaling
  s_height   = 1.1,           # motif height scaling
  h_fill     = "#FAF3E0",     # hex fill (cream)
  h_color    = "#7A4F24",     # hex border (toast brown)
  # url      = "https://cran.r-project.org/package=TOSTER",
  # u_color  = "#2D2A26",
  # u_size   = 4,
  # To embed a toast icon, uncomment and supply a file:
  # subplot  = "toast.png",    # transparent PNG of a toast slice
  filename   = "mess/TOSTER_hex.png",
  dpi        = 600
)


# TOSTER Hex Sticker - Option A with Toast Background
# ====================================================

library(ggplot2)
library(hexSticker)
library(png)
library(grid)
library(cowplot)

# Color palette (toast theme)
colors <- list(
  cream = "#FAF3E0",
  toast_brown = "#7A4F24",
  dark_text = "#2D2A26",
  toast_medium = "#B8844C",
  green_reject = "#6B8E23"
)


# Step 1: Create simplified tnull-style distributions -------

# Parameters for the null distributions
se <- 0.15
lower_bound <- -0.5
upper_bound <- 0.5
alpha <- 0.05

# Create density data for two t-distributions centered at each bound
x_seq <- seq(-1.2, 1.2, length.out = 500)

df_lower <- data.frame(
  x = x_seq,
  y = dt((x_seq - lower_bound) / se, df = 50) / se,
  group = "lower"
)

df_upper <- data.frame(
  x = x_seq,
  y = dt((x_seq - upper_bound) / se, df = 50) / se,
  group = "upper"
)

# Critical values for rejection regions
crit_lower <- lower_bound + qt(alpha, df = 50) * se
crit_upper <- upper_bound + qt(1 - alpha, df = 50) * se

# Rejection region data
reject_lower <- df_lower[df_lower$x <= crit_lower, ]
reject_upper <- df_upper[df_upper$x >= crit_upper, ]

# Build the minimal plot
p_distributions <- ggplot() +
  # Lower null distribution
  geom_area(data = df_lower,
            aes(x = x, y = y),
            fill = colors$toast_medium,
            alpha = 0.5) +
  geom_line(data = df_lower,
            aes(x = x, y = y),
            color = colors$toast_brown,
            linewidth = 0.8) +
  # Upper null distribution
  geom_area(data = df_upper,
            aes(x = x, y = y),
            fill = colors$toast_medium,
            alpha = 0.5) +
  geom_line(data = df_upper,
            aes(x = x, y = y),
            color = colors$toast_brown,
            linewidth = 0.8) +
  # Rejection regions (green)
  geom_area(data = reject_lower,
            aes(x = x, y = y),
            fill = colors$green_reject,
            alpha = 0.7) +
  geom_area(data = reject_upper,
            aes(x = x, y = y),
            fill = colors$green_reject,
            alpha = 0.7) +
  # Equivalence bounds (dashed vertical lines)
  geom_vline(xintercept = lower_bound,
             linetype = "dashed",
             color = colors$dark_text,
             linewidth = 0.7) +
  geom_vline(xintercept = upper_bound,
             linetype = "dashed",
             color = colors$dark_text,
             linewidth = 0.7) +
  # Point estimate + CI at bottom
  annotate("segment",
           x = -0.15, xend = 0.15,
           y = -0.15, yend = -0.15,
           linewidth = 1.2,
           color = colors$dark_text) +
  annotate("point",
           x = 0, y = -0.15,
           size = 2,
           color = colors$dark_text) +
  # Completely blank theme
  theme_void() +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none"
  ) +
  coord_cartesian(xlim = c(-1.1, 1.1), ylim = c(-0.3, 3), expand = FALSE)


# Step 2: Load toast image and create combined plot ---------------


# Read your toast PNG
toast_img <- readPNG("junk/toast.png")
toast_grob <- rasterGrob(toast_img, interpolate = TRUE)

# Create combined plot with toast background
p_combined <- ggplot() +
  annotation_custom(toast_grob,
                    xmin = -Inf, xmax = Inf,
                    ymin = -Inf, ymax = Inf) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA)
  ) +
  coord_fixed(ratio = 1)

# Overlay distributions on toast using cowplot
p_final <- ggdraw() +
  draw_plot(p_combined) +
  draw_plot(p_distributions,
            x = 0.05, y = 0.1,      # position adjustments
            width = 0.9, height = 0.8)


# Step 3: Create hex sticker ------------------


sticker(
  p_final,
  package = "🎉TOSTER",
  p_size = 28,
  p_color = colors$dark_text,
  p_y = 1.4,                    # package name position
  s_x = 1.0,
  s_y = 0.75,
  s_width = 1.5,
  s_height = 1.2,
  h_fill = colors$cream,
  h_color = colors$toast_brown,
  h_size = 1,
  filename = "TOSTER_hex.png",
  dpi = 600
)

# Also save a version for web/preview
sticker(
  p_final,
  package = "TOSTER",
  p_size = 28,
  p_color = colors$dark_text,
  p_y = 1.4,
  s_x = 1.0,
  s_y = 0.75,
  s_width = 1.5,
  s_height = 1.2,
  h_fill = colors$cream,
  h_color = colors$toast_brown,
  h_size = 1.5,
  filename = "TOSTER_hex_preview.png",
  dpi = 150
)

message("Hex sticker saved!")
