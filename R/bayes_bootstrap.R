# library(pbapply)
# # need to install pbapply
# bb_weights <- function(n){
#   wts <- stats::rgamma(n, 1, 1)
#   wts <- wts / sum(wts)
#   return(wts)
# }
#
# bb_lm <- function(formula, data){
#   data <- stats::na.omit(mtcars)
#   n <- nrow(data)
#   wts <- bb_weights(n)
#   fit <- lm(formula, data, weights = wts)
#   return(coefficients(fit))
# }
#
# bb_glm <- function(formula, data,...){
#   data <- stats::na.omit(mtcars)
#   n <- nrow(data)
#   wts <- bb_weights(n)
# }
#
# iter <- 5000
# cl <- parallel::makeCluster(4)
# n <- nrow(mtcars)
# parallel::clusterExport(cl, c("mtcars",
#                               "bb_lm",
#                               "bb_glm",
#                               "bb_weights",
#                               "n"))
# test_lm <- pbapply::pbreplicate(n = iter,
#                                 bb_lm(formula = cbind(mpg, hp) ~ vs,  mtcars),
#                                 cl = cl)
#
# test_glm <- pbapply::pbreplicate(n = iter,
#                                  bb_glm(formula = vs  ~ hp, mtcars, family = "binomial"),
#                                  cl = cl)
#
#
# bb_lm(formula = hp ~ vs,  data = mtcars)
#
