#' Best FPLMBsplines_fit given by a model selection criterion.
#'
#' Fit a FPLM model for different spline basis sizes and picks the best one
#'     according to a specified model selection criterion.
#' 
#' @param y Numeric vector of scalar responses.
#' @param x Matrix of the functional covariates, where each row contains the
#'     functions evaluated on a (common) grid.
#' @param u Numeric vector with the values of the explanatory variable that
#'     enters the model non-parametrically.
#' @param t Numeric vector with the grid over which the functional covariates
#'     were evaluated.
#' @param w Numeric vector with varying coefficients for the non-parametric term
#'     (w = 1 for non varying coefficients).
#' @param range_nonparam_term a vector of B-spline basis sizes to explore for the
#'     functional regression coefficient.
#' @param range_func_term a vector of B-spline basis sizes to explore for the
#'     non-parametric component.
#' @param norder Integer for the order of the B-Splines.
#' @param loss_fun string specifying the loss function. 'ls' for least squares,
#'     'huang' for Huber, 'lmbrob' for MM-estimator.
#' @param criterion Function for the criterion for model selection.
#' @param verbose Logical argument indicating whether partial results are printed.
#'     
#' @return A list including the following components:
#' \describe{
#' \item{fit}{fitted parameters}
#' \item{best_k_npt}{chosen number of splines for the non-parametric term}
#' \item{best_k_ft}{chosen number of splines for the functional regression term}
#' }
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
#' #FIXME: add w
#' y <- x %*% b(t) * min(diff(t)) + w * g(u) + rnorm(n, sd = 0.1)
#'
#' # Best FPLM fit
#' FPLM_fit <- FPLMBsplines(y, x, u, t, w,
#'   range_nonparam_term = 4:13, range_func_term = 4:13,
#'   norder = 4, loss_fun = "ls", criterion = "bic1", verbose = FALSE
#' )
#'
#' # Plot the estimates
#' par(mfrow = c(2, 1))
#' plot(t, FPLM_fit$fit$slope_fun, pch = 16)
#' plot(u, FPLM_fit$fit$eta_est, pch = 16)
#'
#' @import fda robustbase
#' 
#' @export
FPLMBsplines <- function(y, x, u, t, w = 1,
                         range_nonparam_term = range_default,
                         range_func_term = range_default, norder = 4,
                         loss_fun = "lmrob", criterion = rbic,
                         verbose = FALSE) {

    ## Some Setup
    opt <- spl_opt <- freq_opt <- fit_opt <- Inf
    n <- length(y)
    range_default <- floor(max(n^(1 / 5), norder)):
        floor(2 * (norder + n^(1 / 5)))

    ## Double loop
    for (spl in range_func_term) {
        for (freq in range_nonparam_term) {
            fit <- FPLMBsplines_fit(y, x, u, t, w, freq, spl, norder, loss_fun)
            val <- fit$value
            scl <- fit$scale
            crt <- criterion(n, scl, val, spl, freq, criterion)
            if (crt < opt) {
                opt <- crt
                spl_opt <- spl
                freq_opt <- freq
                fit_opt <- fit
            }
            if (verbose) print(c("spl" = spl, "freq" = freq, "crit" = crt))
        }
    }

    ## Best fit
    kns <- seq(min(u), max(u), length = spl_opt - norder + 2)
    base <- create.bspline.basis(
        rangeval = range(u),
        norder = norder,
        breaks = kns
    )
    spl_uu <- getbasismatrix(u, base)
    fit_opt$eta_est <- (spl_uu * w) %*% fit_opt$spl
    dt <- min(diff(t))
    fit_opt$fitted <- as.vector(x %*% fit_opt$slope_fun * dt + fit_opt$eta_est)

    return(list(fit = fit_opt, best_k_npt = spl_opt, best_k_ft = freq_opt))
}
