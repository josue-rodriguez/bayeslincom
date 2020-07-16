
# Code taken from brms:::find_vars
# find_vars <- function (x, dot = TRUE, brackets = TRUE) {
#   regex_all <- paste0("([^([:digit:]|[:punct:])]", if (dot)
#     "|\\.", ")", "[[:alnum:]_\\:", if (dot)
#       "\\.", "]*", if (brackets)
#         "(\\[[^],]+(,[^],]+)*\\])?")
#   pos_all <- gregexpr(regex_all, x)[[1]]
#   regex_fun <- paste0("([^([:digit:]|[:punct:])]", if (dot)
#     "|\\.", ")", "[[:alnum:]_", if (dot)
#       "\\.", "]*\\(")
#   pos_fun <- gregexpr(regex_fun, x)[[1]]
#   pos_decnum <- gregexpr("\\.[[:digit:]]+", x)[[1]]
#   keep <- !pos_all %in% c(pos_fun, pos_decnum)
#   pos_var <- pos_all[keep]
#   attr(pos_var, "match.length") <- attributes(pos_all)$match.length[keep]
#   if (length(pos_var)) {
#     out <- unique(unlist(regmatches(x, list(pos_var))))
#   }
#   else {
#     out <- character(0)
#   }
#   out
# }

# returns strings matching pattern
# code taken from brms:::find_vars
get_matches <- function(pattern, text) {
  match_data <- gregexpr(pattern, text)
  x <- regmatches(text, match_data)
  x <- unlist(x)
  return(x)
}

# checks if interval excludes ROPE
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


# creates string representing a node's
# 'Expected Influence' (Robinaugh et al. 2016)
make_ei_hyp <- function(node, data) {

  # number of columns
  p <- ncol(data)

  # names of columns
  node_names <- colnames(data)

  # index of primary node
  node_idx <- which(node_names == node)

  if (node_idx == 1) {
    hypothesis <- paste(node, node_names[-node_idx], sep = "--", collapse = "+")
  }
  else if (node_idx == p) {
    hypothesis <- paste(node_names[-node_idx], node, sep = "--", collapse = "+")
  }
  else {
    before_idx <- 1:(node_idx-1)
    after_idx <- (node_idx + 1):p

    hyp1 <- paste(node_names[before_idx], node, sep = "--", collapse = "+")
    hyp2 <- paste(node, node_names[after_idx], sep = "--", collapse = "+")
    hypothesis <- paste(hyp1, hyp2, sep = "+")
  }
  return(hypothesis)
}
