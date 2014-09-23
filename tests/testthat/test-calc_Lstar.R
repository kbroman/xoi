
context("calc_Lstar")

test_that("calc_Lstar works", {

    # interference + large chr : Lstar == L
    expect_equal(calc_Lstar(500, 10, 0), 500)

    # p==1 is same as m==0
    expect_equal(calc_Lstar(100, 5, 1), calc_Lstar(100, 0, 0))

    # no interference case
    L <- seq(55, 205, by=25)
    tol <- (.Machine$double.eps)^(1/3)
    for(i in L) {
        Lstar <- calc_Lstar(i, 0, 0)
        expect_equal(i*2, 2*Lstar/(1-dpois(0, Lstar/50)), tolerance=tol)
    }

})

