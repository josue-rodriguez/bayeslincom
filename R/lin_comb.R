#' Perform a linear combination of posterior samples
#'
#' @param lin_comb A character of length one experessing a linear combination
#' @param obj An object of class \code{BGGM}, \code{bbcor}, or \code{data.frame}
#' @param cri A numeric for the 100*(1-a)% credible interval
#' @return An object of class \code{bayeslincom}
#' @examples
#' add(1, 1)
#' @export
lin_comb <- function(lin_comb,
                     obj,
                     cri_level = 0.90,
                     rope = NULL) {
  # blc chcks
  if (length(lin_comb) != 1L) stop("argument 'lin_comb' must have length of 1")
  if (!is.character(lin_comb)) stop("argument 'lin_comb' must be of type 'character'")

  # Extract variable names
  all_vars <- extract_var_names(obj)

  # check if samples are (partial) correlations
  is_corr <- is(obj, "BGGM") || is(obj, "bbcor")

  # Create p by p matrix of name combinations for cors
  if (is_corr) {
    all_vars <- sapply(all_vars,
                       function(x) paste(all_vars, x, sep = "--"))
  }

  # remove whitespace
  # - space in brackets is intentional
  comb <- gsub("[ \t\r\n]", "", lin_comb)

  # extract sign
  sign <- get_matches("=|<|>", comb)

  if (length(sign) != 1L & sign %in% c("=", "<", ">")) {
    stop("LHS and RHS of 'lin_comb' must be separated by '=', '<', or '>'")
  }

  # left and right hand sides of hypothesis
  lr <- get_matches("[^=<>]+", comb)

  # write wrap lhs and rhs with parentheses
  comb <- paste0("(", lr[1], ")")

  comb <- paste0(comb,
              ifelse(lr[2] != "0",
                     yes = paste0("-(", lr[2], ")"),
                     no = ""))

  # extract all variables in combination
  if (is_corr) {
    comb_vars_list <- get_matches("[[:alnum:]]+--[[:alnum:]]+", comb)
    comb_vars <- unlist(comb_vars_list)
  } else {
    comb_vars <- find_vars(comb)
  }

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
    stop(paste("Some variables not found in graph: \n", miss_pars))
  }

  if (is_corr) {
    post_samps <- get_corr_samples(obj, all_vars)
  } else {
    post_samps <- obj
  }

  # evaluate combinations in post/prior
  post_eval <- eval(str2expression(comb_eval), as.data.frame(post_samps))

  # compute credible interval
  a <- (1 - cri_level) / 2
  cri_bounds <- c(a, 1 - a)
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
    rope_overlap <- NULL
  }

  out <- list(lin_comb = lincomb,
              cri = cri,
              cri_level = cri_level,
              rope = rope,
              rope_overlap = rope_overlap,
              samples = post_samps,
              mean_samples = post_mean,
              sd_samples = post_sd,
              support = support,
              call = match.call())

  class(out) <- c("bayeslincom")

  return(out)
}
