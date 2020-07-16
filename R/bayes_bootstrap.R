# need to install pbapply
bb_weights <- function(n){
  wts <- stats::rgamma(n, 1, 1)
  wts <- wts / sum(wts)
  return(wts)
}

bb_lm <- function(formula, data){
  data <- stats::na.omit(mtcars)
  n <- nrow(data)
  return(coefficients(stats::lm(formula = formula,
                                data = data,
                                weights = bb_weights(n))))
}

bb_glm <- function(formula, data,...){
  data <- stats::na.omit(mtcars)
  n <- nrow(data)
  return(coefficients(stats::glm(formula = formula,
                                 data = data,
                                 weights = bb_weights(n), ...)))
}

iter <- 5000
cl <- parallel::makeCluster(4)
n <- nrow(mtcars)
parallel::clusterExport(cl, c("mtcars",
                              "bb_lm",
                              "bb_glm",
                              "bb_weights",
                              "n"))
test_lm <- pbapply::pbreplicate(n = iter,
                                bb_lm(formula = cbind(mpg, hp)  ~ vs,  mtcars),
                                cl = cl)

test_glm <- pbapply::pbreplicate(n = iter,
                                 bb_glm(formula = vs  ~ hp, mtcars, family = "binomial"),
                                 cl = cl)
