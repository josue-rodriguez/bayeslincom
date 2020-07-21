#' Perform a linear combination of posterior samples
#'
#' @param hypothesis A number
#' @param obj A number
#' @param obj An object of class \code{BGGM}, \code{bbcor}, or \code{data.frame}
#' @return An object of class \code{hypothesis}
#' @examples
#' add(1, 1)
#' @export
hypothesis <- function(obj, ...) {
  UseMethod("hypothesis", obj)
}


############################
# ---- BGGM method ----
############################

#' Hypothesis method for BGGM objects
#'
#' @export hypothesis.BGGM
#' @export
hypothesis.BGGM <- function(hypothesis,
                            obj,
                            cri_level = 0.90,
                            rope = NULL) {

  stopifnot(length(hypothesis) == 1L & is.character(hypothesis))

  # Extract variable names
  all_vars <- dimnames(obj$Y)[[2]]

  # Create p by p matrix of name combinations
  all_cors <- sapply(all_vars,
                     function(x) paste(all_vars, x, sep = "--"))

  # remove whitespace, space in bracekets is intentional
  h <- gsub("[ \t\r\n]", "", hypothesis)

  # extract sign
  sign <- get_matches("=|<|>", h)

  stopifnot(length(sign) == 1L & sign %in% c("=", "<", ">"))

  # left and right hand sides of hypothesis
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


  # columns
  p <- obj$p

  # upper triangle indices
  upper_tri <- upper.tri(diag(p))

  # number of iterations in original object
  iter <- obj$iter

  # place posterior samples in matrix
  post <- matrix(
    data = obj$post_samp$pcors[,,51:(iter + 50)][upper_tri],
    nrow = iter,
    ncol = p*(p-1)*0.5,
    byrow = TRUE
    )

  # name all columns
  dimnames(post)[[2]] <- all_vars[upper_tri]

  # evaluate combinations in post/prior
  post_eval <- eval(str2expression(h_eval), as.data.frame(post))
  post_eval_z <- BGGM:::fisher_r_to_z(post_eval)

  prior_eval <- eval(str2expression(h_eval), as.data.frame(prior))
  prior_eval_z <- BGGM:::fisher_r_to_z(prior_eval)

  # bounds for credible interval
  a <- (1 - cri_level) / 2
  cri_bounds <- c(a, 1 - a)

  # compute credible invervals
  cri_z <- quantile(post_eval_z, cri_bounds)
  cri <- quantile(post_eval, cri_bounds)

  # summary statistics
  post_mean <- mean(post_eval)
  post_sd <- sd(post_eval)

  if (!is.null(rope)) {
    # decision rule
    excludes_rope <- excludes_rope(cri, rope, sign)
    rope_overlap <- sum(rope[[1]] < post_eval & post_eval < rope[[2]]) / length(post_eval)

    if (sign != "=") {
      support <- ifelse(excludes_rope,
                        paste0("Test is supported"),
                        paste0("Test is not supported"))
      support <- as.character(support)
    } else {
      support <- ifelse(excludes_rope,
                        paste0("Test is not supported"),
                        paste0("Test is supported"))
      support <- as.character(support)
    }

  } else {
    support <- NULL
  }


  out <- list(hypothesis = hypothesis,
              rope = rope,
              rope_overlap = rope_overlap,
              samples = list(pcors = post_eval,
                             fisher_z = post_eval_z),
              cri = cri,
              cri_z = cri_z,
              mean_samples = post_mean,
              sd_samples = post_sd,
              support = support,
              cri_level = cri_level,
              call = match.call())
  class(out) <- c("bayeslincom", "BGGM")

  return(out)
}

#########################
# ---- BBcor method ----
#########################

#' Hypothesis method for bbcor objects
#'
#' @param hypothesis A number
#' @export hypothesis.bbcor
#' @export
hypothesis.bbcor <- function(hypothesis,
                             obj,
                             cri_level = 0.90,
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
  a <- (1 - cri_level) / 2
  cri_bounds <- c(a, 1 - a)

  # compute  invervals
  cri <- quantile(boot_eval, cri_bounds)

  # summary statistics
  mean_samples <- mean(boot_eval)
  sd_samples <- sd(boot_eval)

  if (!is.null(rope)) {
    # decision rule
    excludes_rope <- excludes_rope(cri, rope, sign)

    if (sign != "=") {
      support <- ifelse(excludes_rope,
                        paste0("test is supported"),
                        paste0("test is not supported"))
      support <- as.character(support)
    } else {
      support <- ifelse(excludes_rope,
                        paste0("test is not supported"),
                        paste0("test is supported"))
      support <- as.character(support)
    }
    rope_overlap <- sum(rope[[1]] < boot_eval & boot_eval < rope[[2]]) / length(boot_eval)
  } else {
    support <- NULL
  }

  out <- list(hypothesis = hypothesis,
              rope = rope,
              samples = boot_eval,
              cri  = cri,
              mean_samples = mean_samples,
              sd_samples = sd_samples,
              support = support,
              cri_level = cri_level,
              call = match.call())


  class(out) <- c("bayeslincom", "bbcor")
  return(out)
}

############################
# ---- data.frame method ----
############################

#' #' Hypothesis method for data.frame objects
#'
#' @param hypothesis A number
#' @export hypothesis.data.frame
#' @export
hypothesis.data.frame <- function(hypothesis,
                                  obj,
                                  cri_level = 0.90,
                                  rope = NULL) {

  stopifnot(length(hypothesis) == 1L & is.character(hypothesis))

  # Extract variable names
  all_vars <- dimnames(obj)[[2]]

  # Create p by p matrix of name combinations
  # all_cors <- sapply(all_vars,
  #                    function(x) paste(all_vars, x, sep = "--"))

  # remove whitespace from hypothesis
  h <- gsub("[ \t\r\n]", "", hypothesis)

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

  # extract all variables in hyp
  vars <- find_vars(h)

  # add backticks for evaluation of hypotheses
  h_eval <- h
  for (v in unique(vars)) {
    h_eval <- gsub(pattern = v,
                   replacement = paste0("`", v, "`"),
                   x = h_eval)
  }

  # check all parameters in hypothesis are valid
  miss_pars <- setdiff(vars, all_vars)
  if (length(miss_pars)) {
    miss_pars <- paste(miss_pars, collapse = ",")
    stop(paste("Some variables not found in graph: \n", miss_pars))
  }

  # evaluate combinations
  samples_eval <- eval(str2expression(h_eval), obj)

  # bounds for interval
  a <- (1 - cri_level) / 2
  cri_bounds <- c(a, 1 - a)

  # compute  invervals
  cri <- quantile(samples_eval, cri_bounds)

  # summary statistics
  mean_samples <- mean(samples_eval)
  sd_samples <- sd(samples_eval)

  if (!is.null(rope)) {
    # decision rule
    excludes_rope <- excludes_rope(cri, rope, sign)

    if (sign != "=") {
      support <- ifelse(excludes_rope,
                        paste0("Test is supported"),
                        paste0("Test is not supported"))
      support <- as.character(support)
    } else {
      support <- ifelse(excludes_rope,
                        paste0("Test is not supported"),
                        paste0("Test is supported"))
      support <- as.character(support)
    }
    rope_overlap <- sum(rope[[1]] < samples_eval & samples_eval < rope[[2]]) / length(samples_eval)
  } else {
    support <- NULL
  }

  out <- list(hypothesis = hypothesis,
              rope = rope,
              samples = samples_eval,
              cri = cri,
              mean_samples = mean_samples,
              sd_samples = sd_samples,
              support = support,
              cri_level = cri_level,
              call = match.call())


  class(out) <- c("bayeslincom")
  return(out)
}
