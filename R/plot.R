#' Perform a linear combination of posterior samples
#'
#' @param x An object of class \code{bayeslincom}
#' @param bins Number of bins
#' @param point_col Color for point indicating mean of posterior
#' @param hist_col Color for histogram edges
#' @param hist_fill Color for histogram bars
#' @param bar_col Color of bar for credible interval
#' @param display_comb_strings If \code{TRUE}, displays full strings for
#'                             combinations in \code{ggplot} facets  when there
#'                             is more than one combination in \code{x}
#' @param ... Currently ignored
#' @return An object of class \code{ggplot}
#' @examples
#'Y <- BGGM::ptsd
#'colnames(Y) <- letters[1:20]
#'est <- BGGM::estimate(Y)
#'bggm_comb <- lin_comb("a--c + a--d > b--c + b--d",
#'                     obj = est,
#'                     cri_level = 0.90,
#'                     rope = c(-0.1, 0.1))
#'plot(bggm_comb)
#' @importFrom  ggplot2 ggplot aes_string geom_histogram geom_point geom_segment
#'                      geom_vline facet_wrap
#' @importFrom stats reshape
#' @export plot.bayeslincom
#' @export
plot.bayeslincom <- function(x,
                             point_col = "black",
                             hist_col = "black",
                             hist_fill = "gray",
                             bar_col = "steelblue",
                             bins = 30,
                             display_comb_strings = TRUE,
                             ...) {
  # ---- Prep plot data ----
  res <- x$results

  # Extract & clean means + CrI dat a
  raw_cri <- extract_list_items(res, "cri", as_df = TRUE)
  raw_means <- extract_list_items(res, "mean_samples", as_df = TRUE)
  raw_means$comb <- row.names(raw_means)
  names(raw_means) <- c("mean", "comb")

  cri_data <- reshape(
    raw_cri,
    ids = c("lb", "ub"),
    direction = "long",
    timevar = "comb",
    v.names = "bounds",
    times = names(raw_cri),
    varying = list(names(raw_cri))
  )

  segment_data <-
    reshape(
      cri_data,
      direction = "wide",
      idvar = "comb",
      v.names = "bounds",
      timevar = "id",
      ids = row.names(cri_data)
    )

  # merge data for CrI bar and mean point
  plot_data <- merge(segment_data, raw_means, by = "comb")

  # Extract & clean samples
  raw_samples <- extract_list_items(res, "samples", as_df = TRUE)
  sample_data <- reshape(
    raw_samples,
    ids = row.names(raw_samples),
    direction = "long",
    timevar = "comb",
    v.names = "samples",
    times = names(raw_samples),
    varying = list(names(raw_samples))
  )

  if (display_comb_strings) {
    comb_strings <- extract_list_items(res, "lin_comb")
    plot_data$comb <- as.factor(plot_data$comb)
    sample_data$comb <- as.factor(sample_data$comb)

    levels(plot_data$comb) <- comb_strings
    levels(sample_data$comb) <- comb_strings
  }

  # ---- Begin plot ----
  p <-
    ggplot() +
    geom_histogram(
      data = sample_data,
      aes_string(x = "samples"),
      col = hist_col,
      fill = hist_fill,
      bins = bins
    ) +
    geom_segment(
      data = plot_data,
      aes_string(x = "bounds.lb",
                 xend = "bounds.ub",
                 y = 0,
                 yend = 0,
                 group = "comb"),
      size = 4,
      col = bar_col
    ) +
    geom_point(
      data = plot_data,
      aes_string(x = "mean",
                 y = 0,
                 group = "comb"),
      col = point_col,
      size = 4
    )

  if (!is.null(x$rope)) {
    p <- p +
      geom_vline(xintercept = x$rope,
                 linetype = "dashed")
  }

  if (length(res) > 1) {
      p <- p + facet_wrap(~ comb, scales = "free")
  }
  return(p)
}
