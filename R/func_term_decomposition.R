#' Projection of the functional data in a finite dimensional B-spline basis
#'
#' Given the size of a B-splines basis computes the corresponding coefficients
#' of the functional observations.
#'
#' @param x Matrix of the functional covariates, where each row contains the
#'     functions evaluated on a (common) grid.
#' @param t Numeric grid over which the functional covariates were evaluated.
#' @param spl_n basis size.
#' @param spl_order spline order (spl_order = 4 for cubic splines).
#'
#' @return A list with components:
#' \describe{
#'   \item{x_coeff}{Matrix with coefficients for projected functional observations,}
#'   \item{spl_basis_eval}{Matrix with spline basis evaluated on the same values as the functional observations.}
#' }
#' 
decompose_functional_term_splines <- function(x, t, spl_n, spl_order) {

    ## Integration step
    dt <- min(diff(t)) # width of grid

    ## Center the observations
    xcenter <- x
    
    ## Spline basis decomposition
    knots_spl <- seq(min(t), max(t), length = spl_n - spl_order + 2)
    spl_basis <- create.bspline.basis(rangeval = range(t),
                                      norder = spl_order,
                                      breaks = knots_spl)
    spl_basis_eval <- getbasismatrix(t, spl_basis)[, 1:spl_n]
    
    ## Estimated Fourier coefficients (by row)
    x_coeff <- xcenter %*% spl_basis_eval * dt

    return(list(x_coeff = x_coeff,
                spl_basis_eval = spl_basis_eval))
    }
