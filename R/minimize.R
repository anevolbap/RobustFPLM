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
minimize <- function(y, x_coeff_ft, u, w, k_npt, k_ft, loss_fun, norder,
                     rob.control = lmrob.control(
                         trace.level = 0,
                         nResample = 5000,
                         tuning.psi = 3.443689, # 85% efficiency
                         subsampling = "simple",
                         rel.tol = 1e-5,
                         refine.tol = 1e-5,
                         k.max = 2e3,
                         maxit.scale = 2e3,
                         max.it = 2e3
                     )) {

    ## B-spline basis
    kns <- seq(min(u), max(u), length = k_npt - norder + 2)
    base <- create.bspline.basis(
        rangeval = range(u),
        norder = norder,
        breaks = kns
    )
    spl_u <- getbasismatrix(u, base)

    ## Design matrix
    INTERCEPT <- !setequal(w, 1)
    if (INTERCEPT) {
        X <- cbind(1, x_coeff_ft, spl_u * w)
    } else {
        X <- cbind(x_coeff_ft, spl_u)
    }

    if (!(loss_fun %in% c("ls", "huang", "lmrob"))) {
        stop("Invalid loss function 'loss_fun'. Should be one of \'ls\' or \'huang\' or \'lmrob\'.")
    }

    ## Minimization menu
    switch(loss_fun, # FIXME: replace switch statement (functions instead of strings)
           ls = {
               fit <- lm(y ~ X - 1)
               cf <- fit$coef
               vv <- sum(fit$res^2) # 1
               ss <- sqrt(mean(fit$res^2))
           },
           huang = {
               init <- lm(y ~ X - 1)$coef
               fit <- huber_estimates(X, y, init, 1.345, 1e-8)
               cf <- as.vector(fit$param)
               ss <- 1
               vv <- fit$value
           },
           lmrob = {
               fit <- lmrob(y ~ X - 1, control = rob.control)
               if (fit$init.S$converged) {
                   cf <- fit$coef
                   ss <- fit$scale
                   vv <- sum(Mpsi(fit$res / ss,
                                  cc = rob.control$tuning.psi,
                                  psi = rob.control$psi,
                                  deriv = -1
                                  ))
               } else {
                   stop("S-estimator did not converge.")
               }
           }
           )

    ## Estimated parameters
    if (INTERCEPT) {
        intercept <- cf[1]
        slope_par <- cf[2:(k_ft + 1)]
        spl_par <- cf[-(1:(k_ft + 1))]
    } else {
        intercept <- 0
        slope_par <- cf[1:k_ft]
        spl_par <- cf[-(1:k_ft)]
    }

    return(
        list(
            estimates_npt = spl_par,
            estimates_ft = slope_par,
            intercept = intercept,
            value = vv,
            scale = ss
        )
    )
}
