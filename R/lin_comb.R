#' Perform a linear combination of posterior samples
#'
#' @param hypothesis A number
#' @param obj A number
#' @param obj An object of class \code{BGGM}, \code{bbcor}, or \code{data.frame}
#' @return An object of class \code{hypothesis}
#' @examples
#' add(1, 1)
#' @export
lin_comb <- function(lin_comb,
                     obj,
                     cri_level = 0.90,
                     rope = NULL) {
  if (length(lin_comb) != 1L) stop("argument 'lin_comb' must have length of 1")
  if (!is.character(lin_comb)) stop("argument 'lin_comb' must be of type 'character'")

  not_supported <- !any(class(obj) %in% c("data.frame", "BGGM", "bbcor"))
  if (not_supported) stop("Currently only objects of type 'data.frame', 'BGGM', and 'bbcor' are supported")

  if (is(obj, "data.frame")) {
    out <- lin_comb.data.frame(lin_comb, obj, cri_level, rope)
  }

  if (is(obj, "BGGM")) {
    out <- lin_comb.BGGM(lin_comb, obj, cri_level, rope)
  }

  if (is(obj, "bbcor")) {
    out <- lin_comb.bbcor(lin_comb, obj, cri_level, rope)
  }

  return(out)
}


############################
# ---- BGGM method ----
############################

lin_comb.BGGM <- function(lin_comb,
                          obj,
                          cri_level = 0.90,
                          rope = NULL) {


  # Extract variable names
  all_vars <- extract_var_names(obj)

  comb <- clean_comb(lin_comb)

  # extract all correlations in hyp
  comb_vars_list <- get_matches("[[:alnum:]]+--[[:alnum:]]+", comb)
  comb_vars <- unlist(comb_vars_list)

  # add backticks for evaluation of hypotheses
  comb_eval <- comb
  for (cv in unique(comb_vars)) {
    comb_eval <- gsub(pattern = cv,
                      replacement = paste0("`", cv, "`"),
                      x = comb_eval)
  }

  # check all parameters in hypothesis are valid
  miss_pars <- setdiff(comb_vars, all_vars)
  if (length(miss_pars)) {
      miss_pars <- paste(miss_pars, collapse = ",")
      stop(paste("Some variables not found in 'obj' \n", miss_pars))
    }

  post_samps <- get_corr_samples(obj, all_vars)

  # evaluate combinations in post/prior
  post_eval <- eval(str2expression(comb_eval), post_samps)
  post_eval_z <- BGGM:::fisher_r_to_z(post_eval)

  # bounds for credible interval
  a <- (1 - cri_level) / 2
  cri_bounds <- c(a, 1 - a)

  cri_z <- quantile(post_eval_z, cri_bounds)
  cri <- quantile(post_eval, cri_bounds)

  post_mean <- mean(post_eval)
  post_sd <- sd(post_eval)

  rope_info <- rope_helper(rope, lin_comb, cri, post_eval)

  out <- list(lin_comb = lin_comb,
              rope = rope,
              rope_overlap = rope_info$rope_overlap,
              samples = list(pcors = post_eval,
                             fisher_z = post_eval_z),
              cri = cri,
              cri_z = cri_z,
              mean_samples = post_mean,
              sd_samples = post_sd,
              support = rope_info$support,
              cri_level = cri_level,
              call = match.call())

  class(out) <- c("bayeslincom")
  return(out)
}

#########################
# ---- BBcor method ----
#########################

#' Hypothesis method for bbcor objects
#'
#' @param hypothesis A number
lin_comb.bbcor <- function(lin_comb,
                           obj,
                           cri_level = 0.90,
                           rope = NULL) {
  all_vars <- extract_var_names(obj)

  comb <- clean_comb(lin_comb)

  # extract all correlations in hyp
  comb_vars_list <- get_matches("[[:alnum:]]+--[[:alnum:]]+", comb)
  comb_vars <- unlist(comb_vars_list)

  # add backticks for evaluation of hypotheses
  comb_eval <- comb
  for (cv in unique(comb_vars)) {
    comb_eval <- gsub(pattern = cv,
                      replacement = paste0("`", cv, "`"),
                      x = comb_eval)
  }

  # check all parameters in hypothesis are valid
  miss_pars <- setdiff(comb_vars, all_vars)
  if (length(miss_pars)) {
    miss_pars <- paste(miss_pars, collapse = ",")
    stop(paste("Some variables not found in 'obj' \n", miss_pars))
  }

  post_samps <- get_corr_samples(obj, all_vars)

  # evaluate combinations in post/prior
  post_eval <- eval(str2expression(comb_eval), post_samps)

  # bounds for credible interval
  a <- (1 - cri_level) / 2
  cri_bounds <- c(a, 1 - a)

  cri <- quantile(post_eval, cri_bounds)

  post_mean <- mean(post_eval)
  post_sd <- sd(post_eval)

  rope_info <- rope_helper(rope, lin_comb, cri, post_eval)

  out <- list(lin_comb = lin_comb,
              rope = rope,
              rope_overlap = rope_info$rope_overlap,
              samples = post_eval,
              cri = cri,
              mean_samples = post_mean,
              sd_samples = post_sd,
              support = rope_info$support,
              cri_level = cri_level,
              call = match.call())

  class(out) <- c("bayeslincom")
  return(out)
}

############################
# ---- data.frame method ----
############################

#' #' Hypothesis method for data.frame objects
#'
#' @param hypothesis A number
lin_comb.data.frame <- function(lin_comb,
                                obj,
                                cri_level = 0.90,
                                rope = NULL) {

  all_vars <- extract_var_names(obj)

  comb <- clean_comb(lin_comb)

  # extract all correlations in hyp
  comb_vars <- find_vars(comb)

  # add backticks for evaluation of hypotheses
  comb_eval <- comb
  for (cv in unique(comb_vars)) {
    comb_eval <- gsub(pattern = cv,
                      replacement = paste0("`", cv, "`"),
                      x = comb_eval)
  }

  # check all parameters in hypothesis are valid
  miss_pars <- setdiff(comb_vars, all_vars)
  if (length(miss_pars)) {
    miss_pars <- paste(miss_pars, collapse = ",")
    stop(paste("Some variables not found in 'obj' \n", miss_pars))
  }

  # evaluate combinations in post/prior
  post_eval <- eval(str2expression(comb_eval), as.data.frame(obj))

  # bounds for credible interval
  a <- (1 - cri_level) / 2
  cri_bounds <- c(a, 1 - a)

  cri <- quantile(post_eval, cri_bounds)

  post_mean <- mean(post_eval)
  post_sd <- sd(post_eval)

  rope_info <- rope_helper(rope, lin_comb, cri, post_eval)

  out <- list(lin_comb = lin_comb,
              rope = rope,
              rope_overlap = rope_info$rope_overlap,
              samples = post_eval,
              cri = cri,
              mean_samples = post_mean,
              sd_samples = post_sd,
              support = rope_info$support,
              cri_level = cri_level,
              call = match.call())

  class(out) <- c("bayeslincom")
  return(out)
}
