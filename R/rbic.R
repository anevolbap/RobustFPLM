#' Robust BIC-type criterion
#'
#' Computes a BIC-type information criterion for the current fit following the
#' proposal in Boente et al. (2020)
#' 
#' @param n Sample size.
#' @param scl Scale estimate (numeric).
#' @param val Minimization value (numeric).
#' @param k_npt Number of parameters for the non-parametric term.
#' @param k_ft Number of parameters for the functional term.
#'
#' @return Information criterion for the selected criterion.
#'
#' @references Boente, Graciela & Salibian-Barrera, Matias & Vena,
#'     Pablo. (2020). Robust estimation for semi-functional linear regression
#'     models. Computational Statistics & Data
#'     Analysis. 152. 107041. 10.1016/j.csda.2020.107041.
rbic <- function(fit) {
    n <- fit$n 
    log(fit$scale^2 * fit$value / n) + (fit$k_npt + fit$k_ft) * log(n) / (n)
}
