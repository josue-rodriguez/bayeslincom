#====== TEST
library(BGGM)
library(dplyr)
library(ggplot2)
library(GGMnonreg)
library(BBcor)

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
tst_bggm <- hypothesis("a--b > 0",
                  obj = est,
                  ci_level = 0.90,
                  # rope = NULL,
                  rope = c(-0.1, 0.1)
                  )
str(tst_bggm)
tst_bggm

#================
x <- GGMnonreg::GGM_bootstrap(Y)

tst_nr <- hypothesis("2*a--b > 0",
           obj = x,
           ci_level = 0.90,
           rope = NULL
           # rope = c(-0.1, 0.1)
           )
str(tst_nr)
tst_nr

#===========


# sample posterior
bb_sample <- bbcor(Y, method = "spearman")

tst_bb <- hypothesis("2*a--b > 0",
                      obj = bb_sample,
                      ci_level = 0.90,
                     # rope = NULL,
                     rope = c(-0.1, 0.1)
                     )
str(tst_bb)
tst_bb
#===========

df_mat <- matrix(
  est$post_samp$pcors[,,51:25050][upper.tri(diag(16))],
  nrow = 25000,
  ncol = 16*15/2,
  byrow = TRUE)
df <- as.data.frame(df_mat)

names(df) <-
  sapply(colnames(Y), function(v) paste0(colnames(Y), "", v))[upper.tri(diag(16))]

tst_df <- hypothesis("ab > 0",
                      obj = df,
                      ci_level = 0.90,
                      rope = NULL
                      # rope = c(-0.1, 0.1)
                     )

str(tst_df)
tst_df



tst_list <- list(tst_bggm, tst_nr, tst_bb, tst_df)
