#' Perform a linear combination of posterior samples
#'
#' @param lin_comb A string specifying a linear combination of variables, or a list of variable names if using \code{contrast}.
#'
#' @param obj An object of class \code{BGGM}, \code{bbcor}, or a \code{data.frame} of posterior samples.
#'
#' @param ci The level for which a credible interval should be computed.
#'
#' @param rope Specify a ROPE. Optional.
#'
#' @param contrast A contrast matrix specifying which combinations to test.
#'
#' @return An object of class \code{lin_comb}
#'
#' @examples
#' # data
#' if (require(BGGM)) library(BGGM)
#' Y <- ptsd
#'
#' # names
#' colnames(Y) <- letters[1:20]
#'
#' # estimate model
#' est <- estimate(Y)
#'
#' # test
#' bggm_comb <- lin_comb("a--c + a--d > b--c + b--d",
#'                        obj = est,
#'                        ci = 0.90,
#'                        rope = c(-0.1, 0.1))
#'
#' # print
#' bggm_comb
#'
#' # Using a contrast matrix to test pairwise differences
#' vars <- c("a--c", "a--d", "b--c")
#'
#' contrast_mat <- matrix(c(1, -1, 0,
#'                          1, 0, -1,
#'                          0, 1, -1), nrow = 3, byrow = TRUE)
#'
#' bggm_comb <- lin_comb(vars,
#'                       obj = est,
#'                       ci = 0.90,
#'                       contrast = contrast_mat)
#'# print
#'bggm_comb
#'
#' @export
#' @importFrom methods is
lin_comb <- function(lin_comb,
                     obj,
                     ci = 0.90,
                     rope = NULL,
                     contrast = NULL) {

  if (methods::is(obj, "data.frame")) {
    out <- lin_comb.data.frame(lin_comb, obj, ci, rope, contrast)
  } else if (methods::is(obj, "BGGM")) {
    out <- lin_comb.BGGM(lin_comb, obj, ci, rope, contrast)
  } else if (methods::is(obj, "bbcor")) {
    out <- lin_comb.bbcor(lin_comb, obj, ci, rope, contrast)
  } else {
    stop("Object class not supported. Must be 'BGGM', 'BBcor', or 'data.frame'")
    }
  return(out)
}


