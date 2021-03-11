context("lin_comb_methods")

test_that("lin_comb methods produce expected errors", {
          Y <- data.frame(x = rnorm(100), y = rnorm(100), z = rnorm(100))

          est <- BGGM::estimate(Y, iter = 200, progress = FALSE)
          # bb <- BBcor::bbcor(Y, iter = 200, cores = 1)

          # valid variables error
          expect_error(bayeslincom:::lin_comb.BGGM("x--a > 2", est),
                       paste("Some variables not found in 'obj' \n", "x--a"))
          # expect_error(bayeslincom:::lin_comb.bbcor("x--a > 2", bb),
                       # paste("Some variables not found in 'obj' \n", "x--a"))
          expect_error(bayeslincom:::lin_comb.data.frame("x > a", Y),
                       paste("Some variables not found in 'obj' \n", "a"))

          # sign error
          sign_error <- "LHS and RHS of 'lin_comb' must be separated by a single '=', '<', or '>'"
          expect_error(bayeslincom:::lin_comb.BGGM("x--y > 2 > 3", est),
                       sign_error)
          # expect_error(bayeslincom:::lin_comb.bbcor("x--a > 2 = x--z", bb),
          #              sign_error)
          expect_error(bayeslincom:::lin_comb.data.frame("x >= y", Y),
                       sign_error)
          })
