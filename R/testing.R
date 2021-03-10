Y = BGGM::bfi[1:50, 1:10]

bggm = BGGM::estimate(Y)
bb = BBcor::bbcor(Y)
obj = bb
ci = 0.90


lin_comb1 = c("A1--A2 > A1--A3",
              "A1--A2 > A2--A3",
              "A1--A3 > A2--A3")
lin_comb = lin_comb1

lin_comb2 = c("A1--A2", "A1--A3", "A2--A3")
contrast = matrix(c(1, -1, 0,
                    1, 0, -1,
                    0, 1, -1), byrow = TRUE, nrow = 3)

tmp <- lin_comb.bbcor(lin_comb2,
                obj,
                contrast = contrast,
                rope = c(-0.1, 0.1))

tmp2 <- lin_comb.BGGM(lin_comb1,
                     obj,
                     rope = c(-0.1, 0.1))

tmp
tmp2
