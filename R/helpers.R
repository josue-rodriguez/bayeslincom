
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

  if (length(sign) != 1L || !sign %in% c("=", "<", ">")) {
    stop("LHS and RHS of 'lin_comb' must be separated by a single '=', '<', or '>'")
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
