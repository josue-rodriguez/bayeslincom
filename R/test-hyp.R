#====== TEST
library(BGGM)
library(dplyr)
library(ggplot2)

Y <- MASS::mvrnorm(n = 500,
                   mu = rep(0, 16),
                   Sigma = ptsd_cor4)

colnames(Y) <- letters[1:16]

est <- estimate(
  Y,
  prior_sd = 0.5,
  type = "continuous",
  seed = NULL,
  iter = 25000,
  progress = F
)

hyps <- sapply(c("a", "b"), function(x) make_ei_hyp(x, Y))
hyp <- paste(hyps, collapse = ">")
tst <- ggmtest("2*a--b > a--b",
                  obj = est,
                  cred = 0.90,
                  rope = c(-0.1, 0.1))
str(tst)
tst

plot(tst) +
  theme_classic()
