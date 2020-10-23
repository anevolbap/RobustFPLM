#' Define the default range for spline basis size exploration
#'
#' A possible way to define the possible sizes of the spline basis is after the
#' rates of convergence given in Theorem 3.2 in Boente et al. (2020) (see
#' Section 2.2) Note that for cubic splines the smallest possible number of
#' knots is 4.
#'
#' @param sample_size Sample size.
#' @param spl_order Integer for spline order (spl_order = 4 for cubic splines).
#'
#' @return vector of basis sizes.
#'
#' @references Boente, Graciela & Salibian-Barrera, Matias & Vena,
#'     Pablo. (2020). Robust estimation for semi-functional linear regression
#'     models. Computational Statistics & Data
#'     Analysis. 152. 107041. 10.1016/j.csda.2020.107041.
#'
#' @examples
#' # For cubic splines and 100 observations the sizes range from 4 to 13:
#' define_default_range(sample_size = 100, spl_order = 4)
#' @export
define_default_range <- function(sample_size, spl_order) {
    lower <- floor(max(sample_size^(1 / 5), spl_order))
    upper <- floor(2 * (spl_order + sample_size^(1 / 5)))
    return(seq.int(lower, upper))
}
