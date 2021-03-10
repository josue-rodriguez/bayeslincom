#' @importFrom stats quantile sd

############################
# ---- BGGM method ---
############################

lin_comb.BGGM <- function(lin_comb,
                          obj,
                          ci = 0.90,
                          rope = NULL,
                          contrast = NULL) {

  if(!requireNamespace("BGGM", quietly = TRUE)) {
    stop("Please install the 'BGGM' package.")
  }

  # container lists to fill
  out <- lc_list <-  list()

  # Extract variable names
  all_vars <- extract_var_names(obj)

  # -- handle contrast matrix --
  if (!is.null(contrast)) {
    post_samps_all <- get_corr_samples(obj, all_vars)
    post_samps <- post_samps_all[, lin_comb]
    post_samps_mat <- as.matrix(post_samps)
    post_eval <- t( contrast %*% t(post_samps_mat) )

    a <- (1 - ci) / 2
    cri_bounds <- c(a, 1 - a)
    for (i in 1:ncol(post_eval)) {
      cri <- quantile(post_eval[, i], cri_bounds)

      post_mean <- mean(post_eval[, i])
      post_sd <- sd(post_eval[, i])
      prob_greater <- mean(post_eval[, i] > 0)

      rope_info <- rope_helper(rope,
                               lin_comb[i],
                               quantiles = cri,
                               post_eval = post_eval[, i],
                               ci = ci)
      out_lc <- list(lin_comb = paste0("C", i),
                     rope_overlap = rope_info$rope_overlap,
                     samples = post_eval[, i],
                     ci = cri,
                     mean_samples = post_mean,
                     sd_samples = post_sd,
                     prob_greater = prob_greater,
                     support = rope_info$support)
      lc_list[[i]] <- out_lc
    }
  } else { # -- If no contrast matrix --
      for (lc_ind in seq_along(lin_comb)) {

        # Extract all correlations in hyp
        comb <- clean_comb(lin_comb[lc_ind])
        comb_vars_list <- get_matches("[[:alnum:]]+--[[:alnum:]]+", comb)
        comb_vars <- unlist(comb_vars_list)

        # add backticks for evaluation of hypotheses
        comb_eval <- comb
        for (cv in unique(comb_vars)) {
          comb_eval <- gsub(pattern = cv,
                            replacement = paste0("`", cv, "`"),
                            x = comb_eval)
        }

        # check all parameters in combination are valid
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
        a <- (1 - ci) / 2
        cri_bounds <- c(a, 1 - a)
        cri <- quantile(post_eval, cri_bounds)

        post_mean <- mean(post_eval)
        post_sd <- sd(post_eval)
        prob_greater <- mean(post_eval > 0)

        rope_info <- rope_helper(rope,
                                 lin_comb[lc_ind],
                                 quantiles = cri,
                                 post_eval = post_eval,
                                 ci = ci)

        out_lc <- list(lin_comb = lin_comb[lc_ind],
                       rope_overlap = rope_info$rope_overlap,
                       samples = post_eval,
                       ci = cri,
                       mean_samples = post_mean,
                       sd_samples = post_sd,
                       prob_greater = prob_greater,
                       support = rope_info$support)
        lc_list[[lc_ind]] <- out_lc
      }
  }

  if (!is.null(contrast)) {
    names(lc_list) <- paste0("C", 1:ncol(post_eval))
  } else {
    names(lc_list) <- paste0("C", seq_along(lin_comb))
    }

  out$results <- lc_list
  out$rope <- rope
  out$ci <- ci
  out$call <- match.call()

  class(out) <- "bayeslincom"
  return(out)
}

#########################
# ---- BBcor method ----
#########################

lin_comb.bbcor <- function(lin_comb,
                           obj,
                           ci = 0.90,
                           rope = NULL,
                           contrast = NULL) {

  if(!requireNamespace("BBcor", quietly = TRUE)) {
    stop("Please install the BBcor package.")
  }

  all_vars <- extract_var_names(obj)

  # initialize empty lists for storing individual results and returned object
  out <- lc_list <- list()

  if (!is.null(contrast)) {
    post_samps_all <- get_corr_samples(obj, all_vars)
    post_samps <- post_samps_all[, lin_comb]
    post_samps_mat <- as.matrix(post_samps)
    post_eval <- t( contrast %*% t(post_samps_mat) )

    a <- (1 - ci) / 2
    cri_bounds <- c(a, 1 - a)
    for (i in 1:ncol(post_eval)) {
      cri <- quantile(post_eval[, i], cri_bounds)

      post_mean <- mean(post_eval[, i])
      post_sd <- sd(post_eval[, i])
      prob_greater <- mean(post_eval[, i] > 0)

      rope_info <- rope_helper(rope,
                               lin_comb[i],
                               quantiles = cri,
                               post_eval = post_eval[, i],
                               ci = ci)
      out_lc <- list(lin_comb = paste0("C", i),
                     rope_overlap = rope_info$rope_overlap,
                     samples = post_eval[, i],
                     ci = cri,
                     mean_samples = post_mean,
                     sd_samples = post_sd,
                     prob_greater = prob_greater,
                     support = rope_info$support)
      lc_list[[i]] <- out_lc
      }
    } else {
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

          # check all parameters in combination are valid
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
          a <- (1 - ci) / 2
          ci_bounds <- c(a, 1 - a)
          cri <- quantile(post_eval, ci_bounds)

          post_mean <- mean(post_eval)
          post_sd <- sd(post_eval)
          prob_greater <- mean(post_eval > 0)

          rope_info <- rope_helper(rope,
                                   lin_comb[lc_ind],
                                   quantiles = cri,
                                   post_eval = post_eval,
                                   ci = ci)

          out_lc <- list(lin_comb = lin_comb[lc_ind],
                         rope_overlap = rope_info$rope_overlap,
                         samples = post_eval,
                         ci = cri,
                         mean_samples = post_mean,
                         sd_samples = post_sd,
                         prob_greater = prob_greater,
                         support = rope_info$support)
          lc_list[[lc_ind]] <- out_lc
        }
    }

  if (!is.null(contrast)) {
    names(lc_list) <- paste0("C", 1:ncol(post_eval))
  } else {
    names(lc_list) <- paste0("C", seq_along(lin_comb))
  }

  out$results <- lc_list
  out$rope <- rope
  out$ci <- ci
  out$call <- match.call()

  class(out) <- "bayeslincom"
  return(out)
  }

############################
# ---- data.frame method ----
############################

lin_comb.data.frame <- function(lin_comb,
                                obj,
                                ci = 0.90,
                                rope = NULL,
                                contrast = NULL) {
  all_vars <- extract_var_names(obj)

  out <- lc_list <- list()


  if (!is.null(contrast)) {
    post_samps <- obj[, lin_comb]
    post_samps_mat <- as.matrix(post_samps)
    post_eval <- t( contrast %*% t(post_samps_mat) )

    a <- (1 - ci) / 2
    cri_bounds <- c(a, 1 - a)
    for (i in 1:ncol(post_eval)) {
      cri <- quantile(post_eval[, i], cri_bounds)

      post_mean <- mean(post_eval[, i])
      post_sd <- sd(post_eval[, i])
      prob_greater <- mean(post_eval[, i] > 0)

      rope_info <- rope_helper(rope,
                               lin_comb[i],
                               quantiles = cri,
                               post_eval = post_eval[, i],
                               ci = ci)
      out_lc <- list(lin_comb = paste0("C", i),
                     rope_overlap = rope_info$rope_overlap,
                     samples = post_eval[, i],
                     ci = cri,
                     mean_samples = post_mean,
                     sd_samples = post_sd,
                     prob_greater = prob_greater,
                     support = rope_info$support)
      lc_list[[i]] <- out_lc
    }
  } else {
    for (lc_ind in seq_along(lin_comb)) {
      # extract all correlations in hyp
      comb <- clean_comb(lin_comb[lc_ind])
      comb_vars <- find_vars(comb)

      # add backticks for evaluation of hypotheses
      comb_eval <- comb
      for (cv in unique(comb_vars)) {
        comb_eval <- gsub(pattern = cv,
                          replacement = paste0("`", cv, "`"),
                          x = comb_eval)
      }

      # check all parameters in combination are valid
      miss_pars <- setdiff(comb_vars, all_vars)
      if (length(miss_pars)) {
        miss_pars <- paste(miss_pars, collapse = ",")
        stop(paste("Some variables not found in 'obj' \n", miss_pars))
      }

      # evaluate combinations in post/prior
      comb_expr <- str2expression(comb_eval)
      post_eval <- eval(comb_expr, as.data.frame(obj))

      # bounds for credible interval
      a <- (1 - ci) / 2
      cri_bounds <- c(a, 1 - a)
      cri <- quantile(post_eval, cri_bounds)

      post_mean <- mean(post_eval)
      post_sd <- sd(post_eval)
      prob_greater <- mean(post_eval > 0)

      rope_info <- rope_helper(rope,
                               lin_comb[lc_ind],
                               quantiles = cri,
                               post_eval = post_eval,
                               ci = ci)

      out_lc <- list(lin_comb = lin_comb[lc_ind],
                     rope_overlap = rope_info$rope_overlap,
                     samples = post_eval,
                     ci = cri,
                     mean_samples = post_mean,
                     sd_samples = post_sd,
                     prob_greater = prob_greater,
                     support = rope_info$support)
      lc_list[[lc_ind]] <- out_lc
    }
  }

  if (!is.null(contrast)) {
    names(lc_list) <- paste0("C", 1:ncol(post_eval))
  } else {
    names(lc_list) <- paste0("C", seq_along(lin_comb))
  }

  out$results <- lc_list
  out$rope <- rope
  out$ci <- ci
  out$call <- match.call()

  class(out) <- "bayeslincom"
  return(out)
}
