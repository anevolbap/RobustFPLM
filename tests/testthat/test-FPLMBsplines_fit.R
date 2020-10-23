test_that("FPLMBplines_fit", {
    
    ## Reality
    set.seed(1)
    n = 300
    p = 100
    ret = FPLMBsplines_fit(y = rnorm(n),
                           x = matrix(rnorm(n * p), n, p),
                           u = sin(sort(seq(0, 1, length = n))),
                           t = sort(runif(p)),
                           w = 1,
                           k_ft = 4,
                           k_npt = 4,
                           spl_order = 4,
                           loss_fun = loss_ls)

    output = c(ret$est_ft, ret$est_npt, use.names = FALSE)

    ## Expected
    slope = c(317.9664, -422.1852, 617.6771, -430.8267)
    spl = c(-0.02783304, 0.03221540, 0.16196770, 0.18592464)
    expected_output = c(slope, spl)

    ## Comparison (up to 4th decimal)
    expect_equal(output, expected_output, tolerance = 1e-4)
})
