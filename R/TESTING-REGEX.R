# # 01. Set up params
# all_vars <- c("a--c", "a_b--d", "b--c", "b--d")
# post_samps <- as.data.frame(replicate(4, rnorm(100)))
# colnames(post_samps) <- all_vars
#
# st <- "a--c + a_b--d > b--c + b--d"
#
#
# # 02. Begin cleaning combination expression
# comb <- clean_comb(st) # add parentheses around string
# # comb <- st
# comb_vars_list <- get_matches("\\(.+?--.+?\\)", comb) # separates out variables
#
# comb_vars_list <- gsub("[()]", "", comb_vars_list)
#
# comb_vars_list1 <- strsplit(comb_vars_list, "[\\+<>\\/\\*]", perl = TRUE)
# #
#
# # GOLD:: look behind positive
# # https://stackoverflow.com/questions/2973436/regex-lookahead-lookbehind-and-atomic-groups
# comb_vars_list2 <- strsplit(unlist(comb_vars_list1), "(?<!-)-(?!-)", perl = TRUE)
# comb_vars <- unlist(comb_vars_list2)
#
# # remove numbers
# suppressWarnings( numeric_idx <- -which(!is.na(as.numeric(comb_vars))))
# if (length(numeric_idx) != 0) {
#   comb_vars <- comb_vars[-numeric_idx]
# }
#
#
#
#
#
#
# # add backticks for evaluation of hypotheses
# comb_eval <- comb
# for (cv in unique(comb_vars)) {
#   comb_eval <- gsub(pattern = cv,
#                     replacement = paste0("`", cv, "`"),
#                     x = comb_eval)
# }
#
# comb_eval
#
# # check all parameters in combination are valid
# miss_pars <- setdiff(unique(comb_vars), all_vars)
# if (length(miss_pars)) {
#   miss_pars <- paste(miss_pars, collapse = ",")
#   stop(paste("Some variables not found in 'obj' \n", miss_pars))
# }
#
# # evaluate combinations in posterior samples
# comb_expr <- str2expression(comb_eval)
# post_eval <- eval(comb_expr, post_samps)
#
#
#
# a <- unique(comb_vars)[3]
#
# gsub(pattern = a,
#      replacement = paste0("`", a, "`"),
#      x = comb_eval)
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# # #-------------------------------------------------
# # Y <- BGGM::ptsd
# #
# # # names
# # colnames(Y) <- letters[1:20]
# #
# # # estimate model
# # est <- BGGM::estimate(Y)
# #
# # # test
# #
# # st <- "a--c + a--d > b--c + b--d"
# # bggm_comb <- lin_comb(st,
# #                        obj = est,
# #                        ci = 0.90,
# #                        rope = c(-0.1, 0.1))
# #
# # # print
# # bggm_comb
# #
# #
# # #===============
# # library(dplyr)
# #
# # path <- 'C:/Users/josue/Downloads/osf_S1_data.csv'
# #
# # d <- read.csv(path) %>% janitor::clean_names()
# #
# # cor_vars <-
# #   d %>%
# #   select(contains("_china"),
# #          contains("_italy"),
# #          contains("_roma"),
# #          china_covid_conspiracy,
# #          generic_covid_conspiracy)
# #
# #
# # Cvec <- c(0, 0, 0, 1, 0, 0, -1, 0, 0, 0,
# #           0, 0, 0, 0, 1, 0, 0, -1, 0, 0,
# #           0, 0, 0, 0, 0, 1, 0, 0, -1, 0 )
# #
# # Cmat <- matrix(Cvec, nrow = 3, byrow = TRUE)
# #
# # # all vars
# # bb_tau <- BBcor::bbcor(cor_vars, method = "pearson", iter = 100, cores = 3)
# # comp_string <- "negative_feelings_china--social_distance_china > social_distance_china--negative_feelings_italy"
# # lin_comb.bbcor(comp_string, obj = bb_tau)
