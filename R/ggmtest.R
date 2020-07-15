ggmtest <- function(hypothesis,
                    obj,
                    cred = 0.90,
                    rope = c(-0.1, 0.1),
                    type = c("default", "ei", "bei")) # not used right now
  {

  stopifnot(length(hypothesis) == 1L & is.character(hypothesis))

  all_vars <- dimnames(obj$Y)[[2]]
  all_cors <- sapply(all_vars,
                     function(x) paste(all_vars, x, sep = "--"))


  h <- gsub("[ \t\r\n]", "", hypothesis)
  sign <- get_matches("=|<|>", h)

  stopifnot(length(sign) == 1L & sign %in% c(">", "<", "="))

  # left and right hand sides
  lr <- get_matches("[^=<>]+", h)

  h <- paste0("(", lr[1], ")")

  h <- paste0(h,
              ifelse(lr[2] != "0",
                     yes = paste0("-(", lr[2], ")"),
                     no = ""))
  vars <- find_vars(h)
  vars_cors <- get_matches("[[:alnum:]]+--[[:alnum:]]+", h)
  # vars_cors <- regmatches(h,
  #                         gregexpr("[[:alnum:]]+--[[:alnum:]]+", h))
  vars_cors <- unlist(vars_cors)

  # add backticks for evaluation of hypotheses
  h_eval <- h
  for (vc in unique(vars_cors)) {
    h_eval <- gsub(pattern = vc,
              replacement = paste0("`", vc, "`"),
              x = h_eval)
  }

  # all parameters in hypothesis are valid
  miss_pars <- setdiff(vars_cors, all_cors)
  if (length(miss_pars)) {
      miss_pars <- paste(miss_pars, collapse = ",")
      stop(paste("Some variables not found in graph: \n", miss_pars))
    }


  # store posterior samples in a matrix
  p <- obj$p

  upper_tri <- upper.tri(matrix(0, p, p))

  iter <- obj$iter

  post <- matrix(
    data = obj$post_samp$fisher_z[,,51:(iter + 50)][upper_tri],
    nrow = iter,
    ncol = p*(p-1)*0.5,
    byrow = TRUE
    )

  # sample prior and store in matrix
  # prior_samps <- BGGM:::sample_prior(
  #   Y = matrix(0, p, p),
  #   iter = iter,
  #   delta = BGGM:::delta_solve(obj$call$prior_sd),
  #   epsilon = 0.01,
  #   prior_only = 1,
  #   explore = 0,
  #   progress = 1
  # )

  # prior <- matrix(
  #   data = prior_samps$fisher_z[,,1:iter][upper_tri],
  #   nrow = iter,
  #   ncol = p*(p-1)*0.5,
  #   byrow = TRUE
  # )

  # name all columns
  dimnames(post)[[2]] <- all_cors[upper_tri]

  # evaluate transformations in post/prior
  post_transform_z <- eval(str2expression(h_eval), as.data.frame(post))
  post_transform <- BGGM:::fisher_z_to_r(post_transform_z)
  # prior_transform <- eval(str2expression(h_eval), as.data.frame(prior))

  # bounds for credible interval
  a <- (1 - cred) / 2
  cri_bound <- c(a, 1 - a)

  # compute credible invervals
  cri_z <- quantile(post_transform_z, cri_bound)
  cri <- quantile(post_transform, cri_bound)
  # q_prior <- quantile(prior_transform, cri_bound)

  # summary statistics
  post_mean <- mean(post_transform)
  post_sd <- sd(post_transform)

  # decision rule
  excludes_rope <- excludes_rope(cri, rope, sign)

  if (sign != "=") {
    support <- ifelse(excludes_rope,
                      paste0("Hypothesis is supported"),
                      paste0("Hypothesis is not supported"))
  }
  else {
    support <- ifelse(excludes_rope,
                      paste0("Hypothesis is not supported"),
                      paste0("Hypothesis is supported"))
  }

  out <- list(hypothesis = hypothesis,
              rope = rope,
              post = post_transform,
              post_z = post_transform_z,
              CrI = cri,
              CrI_z = cri_z,
              support = as.character(support),
              cred = cred,
              call = match.call())
  class(out) <- "ggmtest"
  return(out)
}


print.ggmtest <- function(x, ...) {
  if (is(x, "ggmtest")) {
    cat("ggmtest: Testing Hypotheses in GGMs with Credible Intervals \n\n")
    cat("Call:", deparse(x$call), "\n")
    cat("------ \n")
    cat("Hypothesis:", x$hypothesis, "\n")

    cat("------ \n")


    cri <- round(x$CrI, 2)
    cat(paste0(x$cred*100, "%"), "CrI of the difference: [", cri[[1]], ",", cri[[2]], "] \n")


    cat("ROPE: [", x$rope[[1]], ",", x$rope[[2]], "] \n")

    cat("------ \n")

    cat(x$support)
  }
}

# tst <- ggmtest("2*a--b > a--b",
#                obj = est,
#                cred = 0.95,
#                rope = c(-0.1, 0.1))
#
# tst

