# TOSTER Hex Sticker - Toast Shadow Version (Final)
# ==================================================
library(ggplot2)
library(hexSticker)

# Color palette (toast theme)
colors <- list(
  cream = "#FAF3E0",
  toast_brown = "#7A4F24",
  dark_text = "#2D2A26",
  toast_medium = "#B8844C",
  toast_light = "#D4A574",
  toast_inner = "#F5DEB3",
  green_reject = "#6B8E23"
)

# Parameters for null distributions
se <- 0.18
lower_bound <- -0.5
upper_bound <- 0.5
alpha <- 0.05

# Only generate x values WITHIN equivalence bounds
x_seq <- seq(lower_bound, upper_bound, length.out = 500)

df_lower <- data.frame(
  x = x_seq,
  y = dt((x_seq - lower_bound) / se, df = 50) / se
)

df_upper <- data.frame(
  x = x_seq,
  y = dt((x_seq - upper_bound) / se, df = 50) / se
)

# Critical values for rejection regions (within bounds)
crit_lower <- lower_bound + qt(1 - alpha, df = 50) * se
crit_upper <- upper_bound + qt(alpha, df = 50) * se

# Rejection regions: tails within equivalence bounds
reject_lower <- df_lower[df_lower$x >= crit_lower, ]
reject_upper <- df_upper[df_upper$x <= crit_upper, ]

# Scale distributions
scale_factor <- 0.28
df_lower$y <- df_lower$y * scale_factor
df_upper$y <- df_upper$y * scale_factor
reject_lower$y <- reject_lower$y * scale_factor
reject_upper$y <- reject_upper$y * scale_factor

peak_height <- max(df_lower$y)

# Toast shape
n_top <- 50
top_x <- seq(-1.05, 1.05, length.out = n_top)
top_y <- 1.35 + 0.18 * sqrt(pmax(0, 1 - (top_x / 1.05)^2))

right_x <- c(1.05, 1.08, 1.1, 1.08, 1.05, 1.0, 0.95, 0.92)
right_y <- c(1.35, 1.1, 0.7, 0.3, 0.0, -0.2, -0.35, -0.4)

bottom_x <- seq(0.92, -0.92, length.out = 25)
bottom_y <- rep(-0.42, 25) + 0.02 * cos(pi * bottom_x / 0.92)

left_x <- c(-0.92, -0.95, -1.0, -1.05, -1.08, -1.1, -1.08, -1.05)
left_y <- c(-0.4, -0.35, -0.2, 0.0, 0.3, 0.7, 1.1, 1.35)

toast_x <- c(top_x, right_x, bottom_x, left_x)
toast_y <- c(top_y, right_y, bottom_y, left_y)

toast_outer <- data.frame(x = toast_x, y = toast_y)

center_x <- 0
center_y <- 0.45
toast_inner <- data.frame(
  x = center_x + (toast_x - center_x) * 0.88,
  y = center_y + (toast_y - center_y) * 0.88
)
toast_center <- data.frame(
  x = center_x + (toast_x - center_x) * 0.75,
  y = center_y + (toast_y - center_y) * 0.75
)

(p_toast <- ggplot() +
    # Outer crust
    geom_polygon(data = toast_outer,
                 aes(x = x, y = y),
                 fill = colors$toast_brown) +
    # Inner bread
    geom_polygon(data = toast_inner,
                 aes(x = x, y = y),
                 fill = colors$toast_light) +
    # Lighter center
    geom_polygon(data = toast_center,
                 aes(x = x, y = y),
                 fill = colors$toast_inner,
                 alpha = 0.8) +
    # Lower null distribution (truncated at bounds)
    geom_area(data = df_lower,
              aes(x = x, y = y),
              fill = colors$toast_medium,
              alpha = 0.75) +
    geom_line(data = df_lower,
              aes(x = x, y = y),
              color = colors$toast_brown,
              linewidth = 1.1) +
    # Upper null distribution (truncated at bounds)
    geom_area(data = df_upper,
              aes(x = x, y = y),
              fill = colors$toast_medium,
              alpha = 0.75) +
    geom_line(data = df_upper,
              aes(x = x, y = y),
              color = colors$toast_brown,
              linewidth = 1.1) +
    # Rejection regions
    geom_area(data = reject_lower,
              aes(x = x, y = y),
              fill = colors$green_reject,
              outline.type = "full",
              alpha = 0.85) +
    geom_area(data = reject_upper,
              aes(x = x, y = y),
              fill = colors$green_reject,
              alpha = 0.85) +
    # Equivalence bounds
    geom_segment(aes(x = lower_bound, xend = lower_bound,
                     y = 0, yend = peak_height),
                 linetype = "dashed",
                 color = colors$dark_text,
                 linewidth = 0.5) +
    geom_segment(aes(x = upper_bound, xend = upper_bound,
                     y = 0, yend = peak_height),
                 linetype = "dashed",
                 color = colors$dark_text,
                 linewidth = 0.5) +
    # Point estimate + CI
    annotate("segment",
             x = -0.12, xend = 0.12,
             y = -0.08, yend = -0.08,
             linewidth = 0.75,
             color = colors$dark_text) +
    annotate("point",
             x = 0, y = -0.08,
             size = 1.5,
             color = colors$dark_text) +
    theme_void() +
    theme(
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background = element_rect(fill = "transparent", color = NA)
    ) +
    coord_fixed(ratio = 1, xlim = c(-1.3, 1.3), ylim = c(-0.55, 1.65)))

library(showtext)
library(sysfonts)

# Add Google Fonts (many options available)
font_add_google("Roboto", "roboto")
font_add_google("Montserrat", "montserrat")
font_add_google("Open Sans", "opensans")
font_add_google("Lato", "lato")
font_add_google("Poppins", "poppins")

# Or add a system font like Arial
font_add("arial", "/usr/share/fonts/truetype/msttcorefonts/Arial.ttf")
# On Mac: font_add("arial", "/Library/Fonts/Arial.ttf")
# On Windows: font_add("arial", "C:/Windows/Fonts/arial.ttf")

# Enable showtext
showtext_auto()

# Create hex sticker
s = sticker(
  p_toast,
  package = "TOSTER",
  p_family = "montserrat",
  p_fontface = "bold",
  p_size = 26,
  p_color = colors$dark_text,
  p_x = 1,
  p_y = 0.42,
  s_x = 1.02,
  s_y = 1.02,
  s_width = 1.45,
  s_height = 1.45,
  h_fill = colors$cream,
  h_color = colors$toast_brown,
  h_size = 3,
  filename = "junk/TOSTER_hex_v2_toast_shadow.png",
  dpi = 600
)

message("Sticker saved!")
