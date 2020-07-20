print.hypothesis <- function(x, ...) {
    cat("lhtInt: Testing Linear Hypotheses with Intervals\n")
    cat("Call:\n")
    print(x$call)
    cat("------ \n")


    cat("Hypothesis:", x$hypothesis, "\n")

    if (!is.null(x$rope)) {
      cat("ROPE: [", x$rope[[1]], ",", x$rope[[2]], "] \n")
      cat(x$support, "\n")
    }
    cat("------ \n")

    ci <- round(x$ci, 2)

    cat(paste0(x$ci_level*100, "%"), "CrI for difference: [", ci[[1]], ",", ci[[2]], "] \n")

    cat("Mean of Difference:", round(x$mean_samples, 2), "\n")
    cat("SD of Difference:", round(x$sd_samples, 2), "\n")

    cat("------ \n")
}
