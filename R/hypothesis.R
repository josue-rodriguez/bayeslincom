

# ----- TO DO --------
# + Output in a dataframe
# + Add summary stats for lhs & rhs
# + Allow option to simply take data.frame
# + Expand to lm, BB, etc.
# + Add comments


hypothesis <- function(hypothesis,
                       obj,
                       interval ,
                       rope) {
  UseMethod("hypothesis", obj)
}


############################
# ---- BGGM method ----
############################
hypothesis.BGGM <- function(hypothesis,
                            obj,
                            interval = 0.90,
                            rope = NULL) {

  stopifnot(length(hypothesis) == 1L & is.character(hypothesis))

  # Extract variable names
  all_vars <- dimnames(obj$Y)[[2]]

  # Create p by p matrix of name combinations
  all_cors <- sapply(all_vars,
                     function(x) paste(all_vars, x, sep = "--"))

  # remove whitespace
  h <- gsub("[\\s\t\r\n]", "", hypothesis)

  # extract sign
  sign <- get_matches("=|<|>", h)

  stopifnot(length(sign) == 1L & sign %in% c("=", "<", ">"))

  # left and right hand sides
  lr <- get_matches("[^=<>]+", h)

  # write wrap lhs and rhs with parentheses
  h <- paste0("(", lr[1], ")")

  h <- paste0(h,
              ifelse(lr[2] != "0",
                     yes = paste0("-(", lr[2], ")"),
                     no = ""))

  # extract variables in hyp
  # vars <- find_vars(h)

  # extract all correlations in hyp
  vars_cors_list <- get_matches("[[:alnum:]]+--[[:alnum:]]+", h)
  vars_cors <- unlist(vars_cors_list)

  # add backticks for evaluation of hypotheses
  h_eval <- h
  for (vc in unique(vars_cors)) {
    h_eval <- gsub(pattern = vc,
              replacement = paste0("`", vc, "`"),
              x = h_eval)
  }

  # check all parameters in hypothesis are valid
  miss_pars <- setdiff(vars_cors, all_cors)
  if (length(miss_pars)) {
      miss_pars <- paste(miss_pars, collapse = ",")
      stop(paste("Some variables not found in graph: \n", miss_pars))
    }


  # columns
  p <- obj$p

  # upper triangle indices
  upper_tri <- upper.tri(diag(p))

  # number of iterations in original object
  iter <- obj$iter

  # place posterior samples in matrix
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
  a <- (1 - interval) / 2
  int_bounds <- c(a, 1 - a)

  # compute credible invervals
  int_z <- quantile(post_transform_z, int_bounds)
  int <- quantile(post_transform, int_bounds)

  # summary statistics
  post_mean <- mean(post_transform)
  post_sd <- sd(post_transform)


  if (!is.null(rope)) {
    # decision rule
    excludes_rope <- excludes_rope(interval, rope, sign)

    if (sign != "=") {
      support <- ifelse(excludes_rope,
                        paste0("Test is supported"),
                        paste0("Test is not supported"))
    }
    else {
      support <- ifelse(excludes_rope,
                        paste0("Test is not supported"),
                        paste0("Test is supported"))
    }

  }
  else {
    support <- NULL
  }


  out <- list(hypothesis = hypothesis,
              rope = rope,
              post = post_transform,
              post_z = post_transform_z,
              interval_range = int,
              interval_range_z = int_z,
              post_mean = post_mean,
              post_sd = post_sd,
              support = as.character(support),
              interval = interval,
              call = match.call())
  class(out) <- "hypothesis"

  return(out)
}



############################
# ---- GGMnonreg method ----
############################

hypothesis.GGM_bootstrap <- function(hypothesis,
                                     obj,
                                     interval = 0.90,
                                     rope = NULL)
{

  stopifnot(length(hypothesis) == 1L & is.character(hypothesis))

  # Extract variable names
  all_vars <- dimnames(obj$dat)[[2]]

  # Create p by p matrix of name combinations
  all_cors <- sapply(all_vars,
                     function(x) paste(all_vars, x, sep = "--"))

  # remove whitespace
  h <- gsub("[\\s\t\r\n]", "", hypothesis)

  # extract sign
  sign <- get_matches("=|<|>", h)

  stopifnot(length(sign) == 1L & sign %in% c("=", "<", ">"))

  # left and right hand sides
  lr <- get_matches("[^=<>]+", h)

  # write wrap lhs and rhs with parentheses
  h <- paste0("(", lr[1], ")")

  h <- paste0(h,
              ifelse(lr[2] != "0",
                     yes = paste0("-(", lr[2], ")"),
                     no = ""))

  # extract all correlations in hyp
  vars_cors_list <- get_matches("[[:alnum:]]+--[[:alnum:]]+", h)
  vars_cors <- unlist(vars_cors_list)

  # add backticks for evaluation of hypotheses
  h_eval <- h
  for (vc in unique(vars_cors)) {
    h_eval <- gsub(pattern = vc,
                   replacement = paste0("`", vc, "`"),
                   x = h_eval)
  }

  # check all parameters in hypothesis are valid
  miss_pars <- setdiff(vars_cors, all_cors)
  if (length(miss_pars)) {
    miss_pars <- paste(miss_pars, collapse = ",")
    stop(paste("Some variables not found in graph: \n", miss_pars))
  }

  # bootstrap samples
  boot_samps_list <- lapply(x$boot_results, function(x) x$upper_tri)
  boot_samps <- do.call(rbind, boot_samps_list)

  # get upper.tri indices
  p <- obj$p
  upper_tri <- upper.tri(diag(p))
  # name all columns
  dimnames(boot_samps)[[2]] <- all_cors[upper_tri]

  # evaluate combinations
  boot_comb <- eval(str2expression(h_eval), as.data.frame(boot_samps))

  # bounds for interval
  a <- (1 - interval) / 2
  int_bounds <- c(a, 1 - a)

  # compute  interval
  int <- quantile(boot_comb, int_bounds)

  # summary statistics
  mean_samples <- mean(boot_comb)
  sd_samples <- sd(boot_comb)

  if (!is.null(rope)) {
    # decision rule
    excludes_rope <- excludes_rope(int, rope, sign)

    if (sign != "=") {
      support <- ifelse(excludes_rope,
                        paste0("Test is supported"),
                        paste0("Test is not supported"))
    }
    else {
      support <- ifelse(excludes_rope,
                        paste0("Test is not supported"),
                        paste0("Test is supported"))
    }

  }
  else {
    support <- NULL
  }

  out <- list(hypothesis = hypothesis,
              rope = rope,
              boot_comb = boot_comb,
              interval_range = int,
              mean_samples = mean_samples,
              sd_samples = sd_samples,
              support = as.character(support),
              interval = interval,
              call = match.call())


  class(out) <- "hypothesis"
  return(out)
}

############################
# ---- data.frame method ----
############################
hypothesis.data.frame <- function(hypothesis,
                                  obj,
                                  cred = 0.90,
                                  rope = NULL)
{

  stopifnot(length(hypothesis) == 1L & is.character(hypothesis))

  # Extract variable names
  all_vars <- dimnames(obj)[[2]]

  # Create p by p matrix of name combinations
  all_cors <- sapply(all_vars,
                     function(x) paste(all_vars, x, sep = "--"))

  # remove whitespace
  h <- gsub("[\\s\t\r\n]", "", hypothesis)

  # extract sign
  sign <- get_matches("=|<|>", h)

  stopifnot(length(sign) == 1L & sign %in% c("=", "<", ">"))

  # left and right hand sides
  lr <- get_matches("[^=<>]+", h)

  # write wrap lhs and rhs with parentheses
  h <- paste0("(", lr[1], ")")

  h <- paste0(h,
              ifelse(lr[2] != "0",
                     yes = paste0("-(", lr[2], ")"),
                     no = ""))

  # extract all correlations in hyp
  vars_cors_list <- get_matches("[[:alnum:]]+--[[:alnum:]]+", h)
  vars_cors <- unlist(vars_cors_list)

  # add backticks for evaluation of hypotheses
  h_eval <- h
  for (vc in unique(vars_cors)) {
    h_eval <- gsub(pattern = vc,
                   replacement = paste0("`", vc, "`"),
                   x = h_eval)
  }

  # check all parameters in hypothesis are valid
  miss_pars <- setdiff(vars_cors, all_cors)
  if (length(miss_pars)) {
    miss_pars <- paste(miss_pars, collapse = ",")
    stop(paste("Some variables not found in graph: \n", miss_pars))
  }

  # bootstrap samples
  boot_samps_list <- lapply(x$boot_results, function(x) x$upper_tri)
  boot_samps <- do.call(rbind, boot_samps_list)

  # get upper.tri indices
  p <- obj$p
  upper_tri <- upper.tri(diag(p))
  # name all columns
  dimnames(boot_samps)[[2]] <- all_cors[upper_tri]

  # evaluate combinations
  boot_comb <- eval(str2expression(h_eval), as.data.frame(boot_samps))

  # bounds for interval
  a <- (1 - cred) / 2
  cri_bound <- c(a, 1 - a)

  # compute  invervals
  cri <- quantile(boot_comb, cri_bound)

  # summary statistics
  mean_samples <- mean(boot_comb)
  sd_samples <- sd(boot_comb)

  if (!is.null(rope)) {
    # decision rule
    excludes_rope <- excludes_rope(cri, rope, sign)

    if (sign != "=") {
      support <- ifelse(excludes_rope,
                        paste0("Test is supported"),
                        paste0("Test is not supported"))
    }
    else {
      support <- ifelse(excludes_rope,
                        paste0("Test is not supported"),
                        paste0("Test is supported"))
    }

  }
  else {
    support <- NULL
  }

  out <- list(hypothesis = hypothesis,
              rope = rope,
              boot_comb = boot_comb,
              CrI = cri,
              mean_samples = mean_samples,
              sd_samples = sd_samples,
              support = as.character(support),
              cred = cred,
              call = match.call())


  class(out) <- "hypothesis"
  return(out)
}


#########################
# ---- BBcor method ----
#########################

hypothesis.bbcor <- function(hypothesis,
                             obj,
                             interval = 0.90,
                             rope = NULL)
{

  stopifnot(length(hypothesis) == 1L & is.character(hypothesis))

  # Extract variable names
  all_vars <- dimnames(obj$Y)[[2]]

  # Create p by p matrix of name combinations
  all_cors <- sapply(all_vars,
                     function(x) paste(all_vars, x, sep = "--"))

  # remove whitespace
  h <- gsub("[\\s\t\r\n]", "", hypothesis)

  # extract sign
  sign <- get_matches("=|<|>", h)

  stopifnot(length(sign) == 1L & sign %in% c("=", "<", ">"))

  # left and right hand sides
  lr <- get_matches("[^=<>]+", h)

  # write wrap lhs and rhs with parentheses
  h <- paste0("(", lr[1], ")")

  h <- paste0(h,
              ifelse(lr[2] != "0",
                     yes = paste0("-(", lr[2], ")"),
                     no = ""))

  # extract all correlations in hyp
  vars_cors_list <- get_matches("[[:alnum:]]+--[[:alnum:]]+", h)
  vars_cors <- unlist(vars_cors_list)

  # add backticks for evaluation of hypotheses
  h_eval <- h
  for (vc in unique(vars_cors)) {
    h_eval <- gsub(pattern = vc,
                   replacement = paste0("`", vc, "`"),
                   x = h_eval)
  }

  # check all parameters in hypothesis are valid
  miss_pars <- setdiff(vars_cors, all_cors)
  if (length(miss_pars)) {
    miss_pars <- paste(miss_pars, collapse = ",")
    stop(paste("Some variables not found in graph: \n", miss_pars))
  }

  p <- ncol(obj$Y)
  iter <- obj$iter
  # get upper.tri indices
  upper_tri <- upper.tri(diag(p))

  # bootstrap samples
  boot_samps_list <- sapply(1:iter,
                              function(s) obj$samps[,,s][upper_tri])
  boot_samps <- t(boot_samps_list)

  # name all columns
  dimnames(boot_samps)[[2]] <- all_cors[upper_tri]

  # evaluate combinations
  boot_eval <- eval(str2expression(h_eval), as.data.frame(boot_samps))

  # bounds for interval
  a <- (1 - interval) / 2
  int_bound <- c(a, 1 - a)

  # compute  invervals
  int <- quantile(boot_eval, int_bound)

  # summary statistics
  mean_samples <- mean(boot_eval)
  sd_samples <- sd(boot_eval)

  if (!is.null(rope)) {
    # decision rule
    excludes_rope <- excludes_rope(int, rope, sign)

    if (sign != "=") {
      support <- ifelse(excludes_rope,
                        paste0("Test is supported"),
                        paste0("Test is not supported"))
    }
    else {
      support <- ifelse(excludes_rope,
                        paste0("Test is not supported"),
                        paste0("Test is supported"))
    }

  }
  else {
    support <- NULL
  }

  out <- list(hypothesis = hypothesis,
              rope = rope,
              samples_= boot_eval,
              interval_rannge = int,
              mean_samples = mean_samples,
              sd_samples = sd_samples,
              support = as.character(support),
              interval = interval,
              call = match.call())


  class(out) <- "hypothesis"
  return(out)
}

#########################
# ---- default method ----
#########################
hypothesis.default <- function(hypothesis,
                               obj,
                               interval = 0.90,
                               rope = NULL)
{

  stopifnot(length(hypothesis) == 1L & is.character(hypothesis))

  # Extract variable names
  all_vars <- dimnames(obj)[[2]]


  # remove whitespace
  h <- gsub("[\\s\t\r\n]", "", hypothesis)

  # extract sign
  sign <- get_matches("=|<|>", h)

  stopifnot(length(sign) == 1L & sign %in% c("=", "<", ">"))

  # left and right hand sides
  lr <- get_matches("[^=<>]+", h)

  # write wrap lhs and rhs with parentheses
  h <- paste0("(", lr[1], ")")

  h <- paste0(h,
              ifelse(lr[2] != "0",
                     yes = paste0("-(", lr[2], ")"),
                     no = ""))

  # extract all correlations in hyp
  vars <- find_vars(h)
  # vars_cors_list <- get_matches("[[:alnum:]]+--[[:alnum:]]+", h)
  # vars_cors <- unlist(vars_cors_list)

  # add backticks for evaluation of hypotheses
  h_eval <- h
  for (var in unique(vars)) {
    h_eval <- gsub(pattern = var,
                   replacement = paste0("`", var, "`"),
                   x = h_eval)
  }

  # check all parameters in hypothesis are valid
  miss_pars <- setdiff(vars_cors, all_cors)
  if (length(miss_pars)) {
    miss_pars <- paste(miss_pars, collapse = ",")
    stop(paste("Some variables not found in graph: \n", miss_pars))
  }

  # bootstrap samples
  boot_samps_list <- lapply(x$boot_results, function(x) x$upper_tri)
  boot_samps <- do.call(rbind, boot_samps_list)

  # get upper.tri indices
  p <- obj$p
  upper_tri <- upper.tri(matrix(0, p, p))
  # name all columns
  dimnames(boot_samps)[[2]] <- all_cors[upper_tri]

  # evaluate combinations
  boot_comb <- eval(str2expression(h_eval), as.data.frame(boot_samps))

  # bounds for interval
  a <- (1 - interval) / 2
  int_bounds <- c(a, 1 - a)

  # compute  interval
  int <- quantile(boot_comb, int_bounds)

  # summary statistics
  mean_samples <- mean(boot_comb)
  sd_samples <- sd(boot_comb)

  if (!is.null(rope)) {
    # decision rule
    excludes_rope <- excludes_rope(int, rope, sign)

    if (sign != "=") {
      support <- ifelse(excludes_rope,
                        paste0("Test is supported"),
                        paste0("Test is not supported"))
    }
    else {
      support <- ifelse(excludes_rope,
                        paste0("Test is not supported"),
                        paste0("Test is supported"))
    }
  }
  else {
    support <- NULL
  }

  out <- list(hypothesis = hypothesis,
              rope = rope,
              boot_comb = boot_comb,
              interval_range = int,
              mean_samples = mean_samples,
              sd_samples = sd_samples,
              support = as.character(support),
              interval = interval,
              call = match.call())


  class(out) <- "hypothesis"
  return(out)
}
