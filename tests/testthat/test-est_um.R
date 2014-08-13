
context("est_coi_um")

test_that("estimated intensity is correct in simple cases", {

    xoloc <- list(c(1, 2, 3), 2, c(1, 3))
    sclength <- c(4, 4, 4)
    centromere <- c(2, 2, 2)

    # estimated intensity with window = 0.05
    z05 <- est.coi.um(xoloc, sclength, centromere, intwindow=0.05, intloc=seq(0, 1, len=500))
    pos05 <- z05$intensity[,1]
    expected_intensity <- rep(0, length(pos05))
    expected_intensity[abs(pos05 - 0.25) <= 0.05/2 | abs(pos05 - 0.5) <= 0.05/2 | abs(pos05 - 0.75) <= 0.05/2] <- (2/3)/0.05
    expect_equal(z05$intensity[,2], expected_intensity, tolerance=1e-12)

    # estimated intensity with window = 0.10
    z10 <- est.coi.um(xoloc, sclength, centromere, intwindow=0.10, intloc=seq(0, 1, len=500))
    pos10 <- z10$intensity[,1]
    expected_intensity <- rep(0, length(pos10))
    expected_intensity[abs(pos10 - 0.25) <= 0.10/2 | abs(pos10 - 0.5) <= 0.10/2 | abs(pos10 - 0.75) <= 0.10/2] <- (2/3)/0.10
    expect_equal(z10$intensity[,2], expected_intensity, tolerance=1e-12)
})

