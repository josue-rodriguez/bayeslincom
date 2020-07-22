#' Print formatted summary of a \code{bayeslincom} object
#'
#' @param x An object of class \code{bayeslincom}
#' @param ... Other arguments to be passed to \code{print}
#' @return An object of class \code{bayeslincom}
#' @examples
#' add(1, 1)
#' @export print.bayeslincom
#' @export
print.bayeslincom <- function(x, ...) {
    cri <- round(x$cri, 2)

    print_df <- data.frame(
      mean = round(x$mean_samples, 2),
      sd = round(x$sd_samples, 2),
      cri_lb = cri[[1]],
      cri_ub = cri[[2]]
    )

    cat("bayeslincom: linear combinations of posterior samples\n")
    cat("Call:\n")
    print(x$call)
    cat("------ \n")


    cat("combination:", x$lin_comb, "\n")

    if (!is.null(x$rope)) {
      cat("rope: [", x$rope[[1]], ",", x$rope[[2]], "] \n")

      print_df$overlap <- x$rope_overlap
    }
    cat("------ \n")

    cat("summary of posterior:\n")
    print(print_df, row.names = FALSE, right = T)

    if (!is.null(x$support)) {
    cat("------ \n")
    cat(x$support)
    }
}
