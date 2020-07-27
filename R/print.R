#' Print formatted summary of a \code{bayeslincom} object
#'
#' @param x An object of class \code{bayeslincom}
#' @param ... Other arguments to be passed to \code{print}
#' @return A formatted summary of posterior samples
#' @export print.bayeslincom
#' @export
print.bayeslincom <- function(x, ...) {
    cri <- round(x$cri, 2)

    print_df <- data.frame(
      post_mean = round(x$mean_samples, 2),
      post_sd = round(x$sd_samples, 2),
      cred_lb = cri[[1]],
      cred_ub = cri[[2]]
    )

    cat("bayeslincom: linear combinations of posterior samples\n")
    cat("------ \n")
    cat("Call:\n")
    print(x$call)

    cat("------ \n")

    cat("combination:", x$lin_comb, "\n")

    if (!is.null(x$rope)) {
      cat("rope: [", x$rope[[1]], ",", x$rope[[2]], "] \n")
      print_df$Pr.in <- x$rope_overlap

      # note for ROPE
      note <- "Pr.in: Posterior probability in ROPE"
    } else {
      print_df$pr_less <- round(1 - x$prob_greater, 2)
      print_df$pr_greater<- round(x$prob_greater, 2)

      note <- paste0("Pr.less: Posterior probability less than zero\n",
                     "Pr.greater: Posterior probability less than zero")
    }

    cat("------ \n")

    cat("Posterior Summary:\n\n")
    print(print_df, row.names = FALSE, right = T)
    cat("------ \n")
    cat(paste0("note:\n", note))
    }
