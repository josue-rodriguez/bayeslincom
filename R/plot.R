#' Perform a linear combination of posterior samples
#'
#' @param x An object of class \code{bayeslincom}
#' @param bins Number of bins
#' @param point_col Color for point indicating mean of posterior
#' @param hist_col Color for histogram edges
#' @param hist_fill Color for histogram bars
#' @param bar_col Color of bar for credible interval
#' @return A plot
#' @examples
#'Y <- BGGM::ptsd
#'colnames(Y) <- letters[1:20]
#'est <- BGGM::estimate(Y)
#'bggm_comb <- lin_comb("a--c + a--d > b--c + b--d",
#'                     obj = est,
#'                     cri_level = 0.90,
#'                     rope = c(-0.1, 0.1))
#'plot(bggm_comb)
#' @importFrom  ggplot2 ggplot aes geom_histogram geom_point geom_segment geom_vline
#' @export plot.bayeslincom
#' @export
plot.bayeslincom <- function(x, bins = 30,
                             point_col = "black", hist_col = "black",
                             hist_fill = "gray", bar_col = "steelblue") {
  if (is(x$samples, "list")) {
    plot_data <- data.frame(samples = x$samples$pcors)
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
