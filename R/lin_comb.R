#' Perform a linear combination of posterior samples
#'
#' @param lin_comb A string specifying a linear combination of variables
#'
#' @param obj An object of class \code{BGGM}, \code{bbcor}, or a \code{data.frame} of posterior samples
#'
#' @param cri_level The level for which a credible interval should be computed
#'
#' @param rope Specify a ROPE. Optional.
#'
#' @return An object of class \code{lin_comb}
#'
#' @examples
#' # data
#' Y <- BGGM::ptsd
#'
#' # names
#' colnames(Y) <- letters[1:20]
#'
#' ###########################
#' ######### BGGM ############
#' ###########################
#'
#' # estimate model
#' est <- BGGM::estimate(Y)
#'
#' # test
#' bggm_comb <- lin_comb("a--c + a--d > b--c + b--d",
#'                     obj = est,
#'                     cri_level = 0.90,
#'                     rope = c(-0.1, 0.1))
#'
#' # print
#' bggm_comb
#' @export
#' @importFrom methods is
lin_comb <- function(lin_comb,
                     obj,
                     cri_level = 0.90,
                     rope = NULL) {

  check <- check_lin_comb(lin_comb)

  if (methods::is(obj, "data.frame")) {
    out <- lin_comb.data.frame(lin_comb, obj, cri_level, rope)
  } else if (methods::is(obj, "BGGM")) {
    out <- lin_comb.BGGM(lin_comb, obj, cri_level, rope)
  } else if (methods::is(obj, "bbcor")) {
    out <- lin_comb.bbcor(lin_comb, obj, cri_level, rope)
  } else {
    stop("object class not supported. must be 'BGGM', 'BBcor', or 'data.fram'")
    }
  return(out)
}


