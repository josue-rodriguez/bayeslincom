print.hypothesis <- function(x, ...) {
  if (is(x, "hypothesis")) {
    cat("NAME: Testing Hypotheses in GGMs with Credible Intervals \n\n")
    cat("Call:", deparse(x$call), "\n")

    cat("------ \n")

    cat("Hypothesis:", x$hypothesis, "\n")

    cat("------ \n")

    cri <- round(x$CrI, 2)

    cat(paste0(x$cred*100, "%"), "CrI of the difference: [", cri[[1]], ",", cri[[2]], "] \n")
    cat("ROPE: [", x$rope[[1]], ",", x$rope[[2]], "] \n")


    cat("Posterior Mean:", round(x$post_mean, 2), "\n")
    cat("Posterior SD:", round(x$post_sd, 2), "\n")

    cat("------ \n")

    cat(x$support)
  }
}
