plot.hypothesis <- function(x, bins = 30, hist_col = "black", hist_fill = "gray", bar_col = "tomato") {

  if ("post_samples" %in% names(x)) {
    plot_data <- data.frame(samples = x$post_samples$pcors)
  }
  else {
    plot_data <- data.frame(samples = x$samples)
  }

  p <-
    ggplot(plot_data) +
    geom_histogram(aes(samples),
                   col = hist_col,
                   fill = hist_fill,
                   bins = bins) +
    geom_segment(x = x$ci[[1]],
                 xend = x$ci[[2]],
                 y = 0,
                 yend = 0,
                 col = bar_col,
                 size = 4)

    if (!is.null(x$rope)) {
      p <- p + geom_vline(xintercept = x$rope, linetype = "dashed")
    }
  return(p)
}
