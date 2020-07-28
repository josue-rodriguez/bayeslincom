#' @importFrom stats quantile sd

############################
# ---- BGGM method ---
############################

# Hypothesis method for bbcor objects
lin_comb.BGGM <- function(lin_comb,
                          obj,
                          cri_level = 0.90,
                          rope = NULL) {

  if(!requireNamespace("BGGM", quietly = TRUE)) {
    stop("Please install the 'BGGM' package.")
  }

  # out <- vector("list", length(lin_comb) + 3)
  out <- lc_list <-  list()
  for (lc_ind in seq_along(lin_comb)) {
    # Extract variable names
    all_vars <- extract_var_names(obj)
    comb <- clean_comb(lin_comb[lc_ind])

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

    # evaluate combinations in posterior samples
    comb_expr <- str2expression(comb_eval)
    post_eval <- eval(comb_expr, post_samps)

    # bounds for credible interval
    a <- (1 - cri_level) / 2
    cri_bounds <- c(a, 1 - a)

    cri <- quantile(post_eval, cri_bounds)

    post_mean <- mean(post_eval)
    post_sd <- sd(post_eval)
    prob_greater <- mean(post_eval > 0)

    rope_info <- rope_helper(rope, lin_comb[lc_ind], cri, post_eval)

    out_lc <- list(lin_comb = lin_comb[lc_ind],
                   rope_overlap = rope_info$rope_overlap,
                   samples = post_eval,
                   cri = cri,
                   mean_samples = post_mean,
                   sd_samples = post_sd,
                   prob_greater = prob_greater,
                   support = rope_info$support)
    lc_list[[lc_ind]] <- out_lc

  }
  names(lc_list) <- paste0("C", 1:length(lin_comb))

  out$results <- lc_list
  out$rope <- rope
  out$cri_level <- cri_level
  out$call <- match.call()

  class(out) <- "bayeslincom"
  return(out)
}

#########################
# ---- BBcor method ----
#########################

# Hypothesis method for bbcor objects
lin_comb.bbcor <- function(lin_comb,
                           obj,
                           cri_level = 0.90,
                           rope = NULL) {

  if(!requireNamespace("BBcor", quietly = TRUE)) {
    stop("Please install the BBcor package.")
  }

  all_vars <- extract_var_names(obj)

  # initialize empty lists for storing individual results and returned object
  out <- lc_list <- list()

  for (lc_ind in seq_along(lin_comb)) {
    comb <- clean_comb(lin_comb[lc_ind])

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
    comb_expr <- str2expression(comb_eval)
    post_eval <- eval(comb_expr, post_samps)

    # bounds for credible interval
    a <- (1 - cri_level) / 2
    cri_bounds <- c(a, 1 - a)

    cri <- quantile(post_eval, cri_bounds)

    post_mean <- mean(post_eval)
    post_sd <- sd(post_eval)
    prob_greater <- mean(post_eval > 0)

    rope_info <- rope_helper(rope, lin_comb[lc_ind], cri, post_eval)

    out_lc <- list(lin_comb = lin_comb[lc_ind],
                   rope_overlap = rope_info$rope_overlap,
                   samples = post_eval,
                   cri = cri,
                   mean_samples = post_mean,
                   sd_samples = post_sd,
                   prob_greater = prob_greater,
                   support = rope_info$support)
    lc_list[[lc_ind]] <- out_lc
  }
  names(lc_list) <- paste0("C", 1:length(lin_comb))

  out$results <- lc_list
  out$rope <- rope
  out$cri_level <- cri_level
  out$call <- match.call()

  class(out) <- "bayeslincom"
  return(out)
}

############################
# ---- data.frame method ----
############################

# Hypothesis method for data.frame objects
lin_comb.data.frame <- function(lin_comb,
                                obj,
                                cri_level = 0.90,
                                rope = NULL) {
  all_vars <- extract_var_names(obj)

  out <- lc_list <- list()
  for (lc_ind in seq_along(lin_comb)) {
    comb <- clean_comb(lin_comb[lc_ind])

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
    comb_expr <- str2expression(comb_eval)
    post_eval <- eval(comb_expr, as.data.frame(obj))

    # bounds for credible interval
    a <- (1 - cri_level) / 2
    cri_bounds <- c(a, 1 - a)

    cri <- quantile(post_eval, cri_bounds)

    post_mean <- mean(post_eval)
    post_sd <- sd(post_eval)
    prob_greater <- mean(post_eval > 0)

    rope_info <- rope_helper(rope, lin_comb[lc_ind], cri, post_eval)

    out_lc <- list(lin_comb = lin_comb[lc_ind],
                   rope_overlap = rope_info$rope_overlap,
                   samples = post_eval,
                   cri = cri,
                   mean_samples = post_mean,
                   sd_samples = post_sd,
                   prob_greater = prob_greater,
                   support = rope_info$support)
    lc_list[[lc_ind]] <- out_lc
  }
  names(lc_list) <- paste0("C", 1:length(lin_comb))

  out$results <- lc_list
  out$rope <- rope
  out$cri_level <- cri_level
  out$call <- match.call()

  class(out) <- "bayeslincom"
  return(out)
}

# check_lin_comb <- function(lin_comb){
#
#   if (length(lin_comb) != 1L) {
#     stop("argument 'lin_comb' must have length of 1")
#   }
#   if (!is.character(lin_comb)) {
#     stop("argument 'lin_comb' must be of type 'character'")
#   }
# }

##############
# ROPE helpers
##############

# ---- Determine ROPE info ----
rope_helper <- function(rope, lin_comb, cri, post_eval) {
  if (!is.null(rope)) {
    # decision rule
    sign <- get_sign(lin_comb)
    excludes_rope <- excludes_rope(cri, rope, sign)
    rope_overlap <-
      mean(rope[[1]] < post_eval & post_eval < rope[[2]])

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

  out <- list(support = support, rope_overlap = rope_overlap)
  return(out)
}

# ---- Checks if interval excludes ROPE ----
excludes_rope <- function(quantiles, rope, sign) {
  if (sign == ">") {
    excludes_rope <- quantiles[1] > rope[2]
  }
  else if (sign == "<") {
    excludes_rope <- quantiles[2] < rope[1]
  }
  else {
    excludes_rope <- !(rope[1] < quantiles[1] & quantiles[2] < rope[2])
  }
  return(excludes_rope)
}

####################
# Hypothesis parsers
####################

# ---- Returns sign in combination string ----
get_sign <- function(x) get_matches("=|<|>", x)

# format combination string
clean_comb <- function(lin_comb) {
  # remove whitespace
  # - space is intentional
  comb <- gsub("[ \t\r\n]", "", lin_comb)

  # extract sign
  sign <- get_sign(comb)

  if (length(sign) != 1L & sign %in% c("=", "<", ">")) {
    stop("LHS and RHS of 'lin_comb' must be separated by '=', '<', or '>'")
  }

  # left and right hand sides of hypothesis
  lr <- get_matches("[^=<>]+", comb)

  # wrap lhs and rhs with parentheses
  comb <- paste0("(", lr[1], ")")

  comb <- paste0(comb,
                 ifelse(lr[2] != "0",
                        yes = paste0("-(", lr[2], ")"),
                        no = ""))
  return(comb)
}

# ---- Extract samples from objects of type 'BGGM' or 'bbcor' ----
get_corr_samples <- function(obj, all_vars) {
  p <- ncol(obj$Y)
  upper_tri <- upper.tri(diag(p))
  iter <- obj$iter

  # place posterior samples in matrix
  if (is(obj, "BGGM")) {
    post_samps <- matrix(
      data = obj$post_samp$pcors[,,51:(iter + 50)][upper_tri],
      nrow = iter,
      ncol = p*(p-1)*0.5,
      byrow = TRUE
    )

    # name all columns
    dimnames(post_samps)[[2]] <- all_vars[upper_tri]

  } else {
    post_samps_list <- sapply(1:iter,
                              function(s) obj$samps[,,s][upper_tri])
    post_samps <- t(post_samps_list)

    dimnames(post_samps)[[2]] <- all_vars[upper_tri]
  }
  post_samps <- as.data.frame(post_samps)
  return(post_samps)
}

# ---- Extract variable names ----
extract_var_names <- function(obj, is_corr) {
  # check if samples are (partial) correlations
  is_corr <- is(obj, "BGGM") || is(obj, "bbcor")

  if (is_corr) {
    vars <- dimnames(obj$Y)[[2]]
    var_names <- sapply(vars,
                       function(x) paste(vars, x, sep = "--"))
  } else if (is(obj, "data.frame")) {
    var_names <-  dimnames(obj)[[2]]
  } else {
    stop("Currently only objects of type 'data.frame', 'BGGM', and 'bbcor' are supported")
  }

  return(var_names)
}

# ---- Code from brms:::find_vars ----
find_vars <- function (x) {
  regex_all <- paste0("([^([:digit:]|[:punct:])]",
                      "|\\.", ")", "[[:alnum:]_\\:",
                      "\\.", "]*","(\\[[^],]+(,[^],]+)*\\])?")
  pos_all <- gregexpr(regex_all, x)[[1]]
  regex_fun <- paste0("([^([:digit:]|[:punct:])]",
                      "|\\.", ")", "[[:alnum:]_",
                      "\\.", "]*\\(")
  pos_fun <- gregexpr(regex_fun, x)[[1]]
  pos_decnum <- gregexpr("\\.[[:digit:]]+", x)[[1]]
  keep <- !pos_all %in% c(pos_fun, pos_decnum)
  pos_var <- pos_all[keep]
  attr(pos_var, "match.length") <- attributes(pos_all)$match.length[keep]
  if (length(pos_var)) {
    out <- unique(unlist(regmatches(x, list(pos_var))))
  }
  else {
    out <- character(0)
  }
  return(out)
}

# ---- Code from brms:::find_vars ----
get_matches <- function(pattern, text) {
  match_data <- gregexpr(pattern, text)
  x <- regmatches(text, match_data)
  x <- unlist(x)
  return(x)
}

# ---- Extract and reshape data ----
extract_list_items <- function(x, item, as_df = FALSE) {
  out <- sapply(x, "[[", item)
  if (as_df) out <- as.data.frame(out)
  return(out)
}

globalVariables("samples")
globalVariables("comb")
globalVariables("bounds.lb")
globalVariables("bounds.ub")
