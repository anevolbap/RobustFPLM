#' Best FPLM fit given by a model selection criterion.
#'
#' Pick the best FPLM model across different basis sizes according to a
#' specified model selection criterion.
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
#' @param range_nonparam_term a vector of B-spline basis sizes to explore for
#'     the functional regression coefficient.
#' @param range_func_term a vector of B-spline basis sizes to explore for the
#'     non-parametric component.
#' @param spl_order Integer for the order of the B-Splines.
#' @param loss_fun string specifying the loss function. 'ls' for least squares,
#'     'huang' for Huber, 'lmbrob' for MM-estimator.
#' @param model_selection Function for the criterion for model selection.
#'
#' @return A list including the following components:
#' \describe{
#' \item{est_npt}{estimated spline coefficients for the non-parametric term,}
#' \item{est_ft}{estimated spline coefficients for the functional slope,}
#' \item{est_npt_fun}{estimated non-parametric term,}
#' \item{est_ft_fun}{estimated functional slope,}
#' \item{est_intercept}{estimated intercept,}
#' \item{k_ft}{chosen number of splines for the functional regression term,}
#' \item{k_npt}{chosen number of splines for the non-parametric term,}
#' \item{value}{value of the minimization,}
#' \item{scale}{scale,}
#' \item{y_est}{fitted responses,}
#' \item{n}{sample size.}
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
#' y <- x %*% b(t) * min(diff(t)) + w * g(u) + rnorm(n, sd = 0.1)
#'
#' # Best FPLM fit
#' FPLM_fit <- FPLMBsplines(y, x, u, t, w,
#'                          spl_order = 4,
#'                          loss_fun = ls,
#'                          model_selection = rbic)
#'
#' # Plot the estimates
#' par(mfrow = c(2, 1))
#' plot(t, FPLM_fit$est_ft_fun, pch = 16)
#' plot(u, FPLM_fit$est_npt_fun, pch = 16)
#' @import fda robustbase
#'
#' @export
FPLMBsplines <- function(y, x, u, t, w = 1,
                         range_func_term = range_default,
                         range_nonparam_term = range_default,
                         spl_order = 4,
                         loss_fun = loss_lmrob,
                         model_selection = rbic) {

    ## Set default range
    range_default <- define_default_range(sample_size = length(y),
                                          spl_order = spl_order)

    ## Progress bar
    total_fits <- length(range_nonparam_term) * length(range_func_term)
    message(paste("Searching the best fit among", total_fits, "models."))
    progress_bar <- txtProgressBar(max = total_fits, char = ".", style = 3)
    iter <- 0

    ## Double loop
    opt <- Inf
    for (iter_ft in range_func_term) {
        for (iter_npt in range_nonparam_term) {
            setTxtProgressBar(progress_bar, iter <- iter + 1)
            current_fit <- FPLMBsplines_fit(
                y,
                x,
                u,
                t,
                w,
                iter_ft,
                iter_npt,
                spl_order,
                loss_fun
            )
            goodness <- model_selection(current_fit)
            if (goodness < opt) best_fit <- current_fit
        }
    }
    return(best_fit)
}
