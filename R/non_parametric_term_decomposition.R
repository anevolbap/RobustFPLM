#' Projection of the non-parametric term in a finite dimensional B-spline basis
#'
#' Given the size of a B-splines basis computes the corresponding coefficients
#' of projection of the non-parametric term.
#'
#' @param u Numeric vector with values of the explanatory variable that enters
#'     the model non-parametrically.
#' @param w Numeric vector with varying coefficients for the non-parametric term
#'     (w = 1 for non varying coefficients).
#' @param spl_n basis size.
#' @param spl_order spline order (spl_order = 4 for cubic splines).
#'
#' @return Matrix with spline basis evaluated on the same values as the
#'     functional observations.
#' 
decompose_non_parametric_term_splines <- function(u, w, spl_n, spl_order) {

    knots <- seq(min(u), max(u), length = spl_n - spl_order + 2)
    spl_basis <- create.bspline.basis(rangeval = range(u),
                                      norder = spl_order,
                                      breaks = knots)
    spl_basis_eval <- getbasismatrix(u, spl_basis) * w

    return(spl_basis_eval)
}
