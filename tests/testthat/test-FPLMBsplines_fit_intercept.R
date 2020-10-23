test_that("FPLMBplines_fit_tecator_lmrob", {
    
    ## Reality
    set.seed(124)
    data(tecator, package = "fda.usc")
    ab <- tecator$absorp.fdata
    ab2 <- fda.usc::fdata.deriv(ab, nderiv = 2) # segunda derivada
    train <- 1:155
    y <- tecator$y$Fat[train]
    u <- tecator$y$Protein[train]
    t <- ab$argvals
    w <- tecator$y$Water[train]
    ab2$data <- ab2$data[train, ]
    x <- ab2$data
    ret = FPLMBsplines_fit(y = y,
                           x = x,
                           u = u,
                           t = t,
                           w = w,
                           k_ft = 5,
                           k_npt = 7,
                           norder = 4,
                           loss_fun = loss_lmrob)
    output = c(ret$est_intercept, ret$est_ft, ret$est_npt, use.names = FALSE)

    ## Expected
    slope = c(575.921977, -0.186321, -90.393038, 257.141971,-129.948037)
    spl = c(-0.9981523, -0.9970536, -1.0035561, -1.0525582, -1.0649512, -1.0880852, -1.1045347)
    intercept = 86.9019 
    expected_output = c(intercept, slope, spl)

    ## Comparison (up to 4th decimal)
    expect_equal(output, expected_output, tolerance = 1e-4)
})



