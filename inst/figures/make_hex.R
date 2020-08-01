library(ggplot2)
library(dplyr)
library(hexSticker)
library(showtext)

set.seed(1)
samps <- rnorm(1e4, mean = 0, sd = 1/3)
post_samps <- tibble(samples = samps)

post_samps %>%
  ggplot(aes(samples)) +
  geom_histogram(bins = 30,
                 col = "black",
                 fill = "#FFC72C",
                 size = 1) +
  # CrI
  geom_segment(x = -0.8,
               xend = 0.8,
               y = 0,
               yend = 0,
               col = "black",
               size = 4) +
  # Point
  geom_point(x = 0,
             y = 0,
             size = 3.5,
             col = "white") +
  # ROPE
  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed",
             col = "white",
             size = 1) +
  # Plot settings
  labs(x = "",
       y = "") +
  scale_x_continuous(labels = NULL) +
  scale_y_continuous(labels = NULL) +
  theme_void() +
  theme_transparent() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 36, color = "white")
  )  -> rope_plot

# rope_plot
ggsave(
  "inst/figures/hex_plot.png",
  rope_plot,
  bg = "transparent",
  width = 15,
  height = 7.5,
  units = "cm",
  dpi = 300)

# fonts
# font_add_google("Fira Code", "fira")
# showtext_auto()

# rope_plot
sticker(
  "inst/figures/hex_plot.png",
  # package settings
  package = "2*bayes+lin=com",
  p_family = "fira",
  p_color = "white",
  p_size = 12,
  p_y = 1.5,
  # plot settings
  s_x = 1.02,
  s_y = 0.9,
  s_width =  0.9,
  s_height = 0.9,
  # hex settings
  h_fill = "#046A38",
  h_color = "#FFC72C",
  h_size = 2.25,
  filename = "inst/figures/hex_sticker.png",
  dpi = 320
)

# a test
