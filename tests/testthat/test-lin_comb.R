context("lin_comb")

test_that("lin_comb produces expected errors", {
  Y <- data.frame(x = rnorm(100), y = rnorm(100))

  expect_error(lin_comb("2*x > y", obj = "a string"),
               "Object class not supported. Must be 'BGGM', 'BBcor', or 'data.frame'")
})
