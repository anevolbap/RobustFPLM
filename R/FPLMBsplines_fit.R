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
#' @param w Numeric vector with varying coefficients for the non-parametric term
#'     (w = 1 for non varying coefficients).
#' @param k_ft basis size for the functional regression term.
#' @param k_npt basis size for the non-parametric term.
#' @param spl_order Integer for spline order (norder = 4 for cubic splines).
#' @param loss_fun String specifying the loss function. 'ls' for least squares,
#'     'huang' for Huber, 'lmrob' for MM-estimator.
#'
#' @return A list with components:
#' \describe{
#'   \item{estimates_npt}{estimate for the non parametric term}
#'   \item{estimates_ft}{estimate for the functional slope}
#'   \item{value}{the minimum value of the loss function}
#'   \item{scale}{the scale estimate}
#'   \item{slope_fun}{functional slope estimate}
#' }
#'
#' @examples
#'
#' # Synthetic data
#' n <- 100
#' m <- 50
#' u <- runif(n)
#' t <- runif(m)
#' w <- 1
#' b <- function(x) x^3
#' g <- function(x) sin(x)
#' x <- matrix(rnorm(n * m), nrow = n)
#' y <- x %*% b(t) * min(diff(t)) + w * g(u) + rnorm(n, sd = 0.1)
#'
#' # Robust FPLM fit for k_ft = 5 and k_npt = 7 with cubic splines
#' FPLM_one_fit <- FPLMBsplines_fit(y, x, u, t, w,
#'                                  k_ft = 5, k_npt = 7,
#'                                  spl_order = 4,
#'                                  loss_fun = loss_lmrob,
#'                                  model_selection = rbic)
#'                   
#' # Plot the estimates
#' par(mfrow = c(2, 1))
#' plot(t, FPLM_one_fit$est_ft_fun, pch = 16)
#' plot(u, FPLM_one_fit$est_npt_fun, pch = 16)
#' 
#' @import fda robustbase
#' 
#' @references Boente, Graciela & Salibian-Barrera, Matias & Vena,
#'     Pablo. (2020). Robust estimation for semi-functional linear regression
#'     models. Computational Statistics & Data
#'     Analysis. 152. 107041. 10.1016/j.csda.2020.107041.
#' @export
FPLMBsplines_fit <- function(y, x, u, t, w, k_ft, k_npt, spl_order, loss_fun, ...) {

    ## Functional slope decomposition
    func_term <- decompose_functional_term_splines(x, t, k_ft, spl_order)
    x_coeff <- func_term$x_coeff

    ## Non-parametric term decomposition
    non_parametric_term_spl_basis <- decompose_non_parametric_term_splines(u, w, k_npt, spl_order)
    spl_u <- non_parametric_term_spl_basis
        
    ## Check if an intercept is needed
    intercept <- !setequal(w, 1)

    ## Minimize loss function
    estimates <- minimize(y,
                          x_coeff,
                          spl_u,
                          intercept,
                          loss_fun, ...)
   
    ## Estimates
    estimates$est_npt_fun <- spl_u %*% estimates$est_npt
    estimates$est_ft_fun <- func_term$spl_basis_eval %*% estimates$est_ft +
        estimates$est_intercept
    estimates$y_est <- as.vector(x %*% estimates$est_ft_fun * min(diff(t)) +
                                 estimates$est_npt_fun + estimates$est_intercept)
    estimates$k_ft <- k_ft
    estimates$k_npt <- k_npt
    estimates$n <- length(y)
    
    return(estimates)
}
