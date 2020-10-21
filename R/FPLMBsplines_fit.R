#' FPLM fit with fixed B-splines basis sizes.
#'
#' Given the sizes of both B-splines basis computes the corresponding FPLM fit.
#'
#' @param y Numeric vector of scalar responses.
#' @param x Matrix of the functional covariates, where each row contains the
#'     functions evaluated on a (common) grid.
#' @param u Numeric vector with values of the explanatory variable that enters
#'     the model non-parametrically.
#' @param t Numeric grid over which the functional covariates were evaluated.
#' @param k_ft basis size for the functional regression term.
#' @param k_npt basis size for the non-parametric term.
#' @param norder Integer for spline order (norder = 4 for cubic splines).
#' @param loss_fun String specifying the loss function. 'ls' for least squares,
#'     'huang' for Huber, 'lmrob' for MM-estimator.
#'
#' @return A list with components:
#' \describe{
#'   \item{estimates_npt}{estimate for the non parametric term,}
#'   \item{estimates_ft}{estimate for the functional slope,}
#'   \item{value}{the minimum value of the loss,}
#'   \item{scale}{the scale estimate.}
#'   \item{slope_fun}{functional slope estimate}.
#' }
#'
#' @import fda robustbase
#' 
#' @references Boente, Graciela & Salibian-Barrera, Matias & Vena,
#'     Pablo. (2020). Robust estimation for semi-functional linear regression
#'     models. Computational Statistics & Data
#'     Analysis. 152. 107041. 10.1016/j.csda.2020.107041.
#' @export
FPLMBsplines_fit <- function(y, x, u, t, w, k_ft, k_npt, norder, loss_fun) {

    ## Integration step
    dt <- min(diff(t)) # width of grid
    xcenter <- x

    ## Spline basis decomposition
    nodos_spl <- seq(min(t), max(t), length = k_ft - norder + 2)
    base_spl <- create.bspline.basis(rangeval = range(t),
                                     norder = norder,
                                     breaks = nodos_spl)
    beta_spl <- getbasismatrix(t, base_spl)
    beta_spl_red <- beta_spl[, 1:k_ft]
    
    ## Estimated Fourier coefficients (by row)
    x_coeff_ft <- xcenter %*% beta_spl_red * dt

    ## Parameter estimation
    est <- minimize(y, x_coeff_ft, u, w, k_npt, k_ft, loss_fun, norder)
    est$slope_fun <- beta_spl_red %*% est$estimates_ft + est$intercept
    
    return(est)
}
