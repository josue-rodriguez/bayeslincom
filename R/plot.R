plot.ggmtest <- function(x, bins = 30, hist_col = "black", hist_fill = "gray", bar_col = "tomato") {

  plot_data <- data.frame(samples = x$post)

  ggplot(plot_data) +
    geom_histogram(aes(samples),
                   col = hist_col,
                   fill = hist_fill,
                   bins = bins) +
    geom_segment(x = x$CrI[[1]],
                 xend = x$CrI[[2]],
                 y = 0,
                 yend = 0,
                 col = bar_col,
                 size = 4) +
    geom_vline(xintercept = x$rope,
               linetype = "dashed")
}

# tst
# plot(tst) +
#   theme_bw()

