loss_ls = function(y, X) {
    fit <- lm(y ~ X - 1)
    coeff <- fit$coef
    value <- sum(fit$res^2)
    scale <- sqrt(mean(fit$res^2))
    return(list(coeff = coeff,
                value =value,
                scale = scale))
}

loss_huang = function(y, X) {
    init <- lm(y ~ X - 1)$coef
    fit <- huber_estimates(X, y, init, 1.345, 1e-8)
    coeff <- as.vector(fit$param)
    scale <- 1
    value <- fit$value
    return(list(coeff = coeff,
                value =value,
                scale = scale))

}

loss_lmrob = function(y, X, ...) {
    
    if (!("rob.control" %in% names(list(...)))) {      
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
        )
    }
    fit <- lmrob(y ~ X - 1, control = rob.control)
    if (!fit$init.S$converged) stop("S-estimator did not converge.")
    coeff <- fit$coef
    scale <- fit$scale
    value <- sum(Mpsi(fit$res / scale,
                      cc = rob.control$tuning.psi,
                      psi = rob.control$psi,
                      deriv = -1))
    return(list(coeff = coeff,
                value =value,
                scale = scale))
}
