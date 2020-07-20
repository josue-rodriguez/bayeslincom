print.hypothesis <- function(x, ...) {
    ci <- round(x$ci, 2)

    print_df <- data.frame(
      mean = round(x$mean_samples, 2),
      sd = round(x$sd_samples, 2),
      ci_lb = ci[[1]],
      ci_ub = ci[[2]]
    )

    cat("lhInt: Testing Linear Combinations with Intervals\n")
    cat("Call:\n")
    print(x$call)
    cat("------ \n")


    cat("hypothesis:", x$hypothesis, "\n")

    if (!is.null(x$rope)) {
      cat("rope: [", x$rope[[1]], ",", x$rope[[2]], "] \n")

      print_df$overlap <- x$rope_overlap
    }
    cat("------ \n")


    print(print_df, row.names = FALSE, right = T)

    if (!is.null(x$support)) {
    cat("------ \n")
    cat(x$support)
    }

}
