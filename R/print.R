#' Print formatted summary of a \code{bayeslincom} object
#'
#' @param x An object of class \code{bayeslincom}
#' @param ... Other arguments to be passed to \code{print}
#' @return A formatted summary of posterior samples
#' @export print.bayeslincom
#' @export
print.bayeslincom <- function(x, ...) {
  res <- x$results

  cri_raw <- extract_list_items(res, "cri")
  cri <- round(cri_raw, 2)

  Post.mean_raw <- extract_list_items(res, "mean_samples")
  Post.mean <- round(Post.mean_raw, 2)

  Post.sd_raw <- extract_list_items(res, "sd_samples")
  Post.sd <- round(Post.mean_raw, 2)

  print_df <- data.frame(
    Post.mean = Post.mean,
    Post.sd = Post.sd,
    Cred.lb = cri[1, ],
    Cred.ub = cri[2, ]
  )
  row.names(print_df) <- names(x$results)

  # ---- Begin pasting output ----
  cat("bayeslincom: Linear Combinations of Posterior Samples\n")
  cat("------ \n")
  cat("Call:\n")
  print(x$call)

  cat("------ \n")

  cat("Combinations:\n")
  comb_list <- extract_list_items(res, "lin_comb")

  for (comb in seq_along(comb_list)) {
    cat(paste0(" C", comb, ":"), comb_list[[comb]], "\n")
  }
  cat("------ \n")

  cat("Posterior Summary:\n\n")

  if (!is.null(x$rope)) {
    cat("ROPE: [", x$rope[[1]], ",", x$rope[[2]], "] \n\n")
    print_df$Pr.in <- extract_list_items(res, "rope_overlap")

    # note for ROPE
    note <- "Pr.in: Posterior probability in ROPE"
  } else {
    prob_greater <- extract_list_items(res, "prob_greater")
    print_df$Pr.less <- round(1 - prob_greater, 2)
    print_df$Pr.greater<- round(prob_greater, 2)

    note <- paste0("Pr.less: Posterior probability less than zero\n",
                   "Pr.greater: Posterior probability greater than zero")
  }

  print(print_df, right = T)
  cat("------ \n")
  cat(paste0("Note:\n", note))
  }
