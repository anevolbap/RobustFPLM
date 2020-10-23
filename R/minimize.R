#' Minimization routine
#'
#' This function is the minimization routine used internally by
#' \code{FPLMBsplines_fit()}, not meant for 'public' use.
#'
#' @param y Numeric vector of scalar responses.
#' @param x_coeff_ft Matrix with coefficients for the functional term
#'     decomposition.
#' @param u Numeric vector with values of the explanatory variable that enters
#'     the model non-parametrically.
#' @param w Numeric vector with varying coefficients for the non-parametric
#'     term (w = 1 for non varying coefficients).
#' @param k_npt Number of splines for the non-parametric term.
#' @param k_ft Number of splines for the functional term.
#' @param loss_fun String specifying the loss function. 'ls' for least squares,
#'     'huang' for Huber, 'lmrob' for MM-estimator.
#' @param norder Integer for spline order (norder = 4 for cubic splines).
#' @param rob.control Optional parameters if \code{loss_fun == "lmbrob"}.
#'
#' @return A list with components:
#' \describe{
#'   \item{estimates_npt}{estimate for the non parametric term,}
#'   \item{estimates_ft}{estimate for the functional slope,}
#'   \item{value}{the minimum value of the loss,}
#'   \item{scale}{the scale estimate.}
#' }
#'
#' @examples
#'
#' ## Synthetic data
#' y <- rnorm(100)
#' x <- matrix(rnorm(5 * 100), ncol = 5)
#' u <- runif(100)
#' w <- 1
#'
#' ## Least-squares estimation
#' minimize(y, x, u, w, k_npt = 5, k_ft = 6, loss_fun = "ls", norder = 4)
#'
#' ## M-estimation
#' minimize(y, x, u, w, k_npt = 5, k_ft = 6, loss_fun = "huang", norder = 4)
#'
#' ## MM-estimation
#' minimize(y, x, u, w, k_npt = 5, k_ft = 6, loss_fun = "lmrob", norder = 4)
#'
#' ## MM-estimation with 95% of efficiency
#' minimize(y, x, u, w,
#'   k_npt = 5, k_ft = 6, loss_fun = "lmrob", norder = 4,
#'   rob.control = list(tuning.psi = 4.685061)
#' )
#' @import fda robustbase
#' @importFrom stats lm
minimize <- function(y, x_coeff, spl_u, intercept, loss_fun, ...) {
  
    ## Design matrix
    X <- cbind(x_coeff, spl_u)
    if (intercept) X <- cbind(1, X)
    
    loss <- loss_fun(y, X, ...)                             

    ## Estimated parameters
    k <- ncol(x_coeff)
    if (intercept) {
        intercept <- loss$coeff[1]
        slope_par <- loss$coeff[2:(k + 1)]
        spl_par <- loss$coeff[-(1:(k + 1))]
    } else {
        intercept <- 0
        slope_par <- loss$coeff[1:k]
        spl_par <- loss$coeff[-(1:k)]
    }

    return(
        list(
            est_npt = spl_par,
            est_ft = slope_par,
            est_intercept = intercept,
            value = loss$value,
            scale = loss$scale
        )
    )
}

##          if (!(loss_fun %in% c("ls", "huang", "lmrob"))) {
##              stop("Invalid loss function 'loss_fun'. Should be one of \'ls\' or \'huang\' or \'lmrob\'.")
##          }
