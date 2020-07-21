#' Perform a linear combination of posterior samples
#'
#' @param hypothesis A number
#' @param obj A number
#' @param obj An object of class \code{BGGM}, \code{bbcor}, or \code{data.frame}
#' @return An object of class \code{hypothesis}
#' @importFrom  ggplot2 ggplot aes geom_histogram geom_point geom_segment geom_vline
#' @examples
#' add(1, 1)
#' @export plot.bayeslincom
#' @export
plot.bayeslincom <- function(x, bins = 30,
                             point_col = "black", hist_col = "black",
                             hist_fill = "gray", bar_col = "steelblue") {


  if (is(x, "BGGM")) {
    plot_data <- data.frame(samples = x$post_samples$pcors)
  } else {
    plot_data <- data.frame(samples = x$samples)
  }

  p <-
    ggplot(plot_data) +
    geom_histogram(aes(samples),
                   col = hist_col,
                   fill = hist_fill,
                   bins = bins) +
    geom_segment(x = x$cri[[1]],
                 xend = x$cri[[2]],
                 y = 0,
                 yend = 0,
                 col = bar_col,
                 size = 4) +
    geom_point(x = x$mean_samples,
               y = 0,
               col = point_col,
               size = 4)

    if (!is.null(x$rope)) {
      p <- p +
        geom_vline(xintercept = x$rope,
                   linetype = "dashed")
    }
  return(p)
}
