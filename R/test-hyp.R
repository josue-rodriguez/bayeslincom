#====== TEST
library(BGGM)
library(dplyr)
library(ggplot2)

Y <- MASS::mvrnorm(n = 500,
                   mu = rep(0, 16),
                   Sigma = ptsd_cor4)

colnames(Y) <- letters[1:16]

est <- explore(
  Y,
  prior_sd = 0.5,
  type = "continuous",
  seed = NULL,
  iter = 25000,
  progress = F
)

hyps <- sapply(c("a", "b"), function(x) make_ei_hyp(x, Y))
hyp <- paste(hyps, collapse = ">")
tst <- hypothesis("2*a--b > a--b",
                  obj = est,
                  interval = 0.90,
                  rope = NULL)
str(tst)
tst


library(GGMnonreg)
x <- GGMnonreg::GGM_bootstrap(Y)

tst <- hypothesis("2*a--b > 0",
           obj = x,
           interval = 0.90,
           rope = c(-0.1, 0.1))
str(tst)
tset

#===========
library(BBcor)
Y <- mtcars[,1:5]

# sample posterior
bb_sample <- bbcor(Y, method = "spearman")
class(bb_sample)


# correlation matrix
dimnames(bb_sample$samps)[[2]] <- letters[1:5]
str(bb_sample$samps)

#===========
X <- matrix(rnorm(10 * 5000), nrow = 5000, ncol = 10)
df <- as.data.frame(X)
head(df)
dimnames(df)[[2]]

tst_df <- hypothesis.default("2*V1 > V2",
                              obj = df,
                              interval = 0.90,
                              rope = NULL)
str(tst_df)
