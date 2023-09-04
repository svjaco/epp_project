################################################################################
#              Tests for bw_CV_fun function (CV optimal bandwidth)             #
################################################################################

# We use the package "locpol" (published on CRAN) to validate our
# implementation of the cross-validated bandwidth selector.

test_that("warning messages appear", {

    expect_warning(bw_CV_fun(X = X, Y = Y, bw_grid = seq(0.075, 0.2, length.out = 100)),
                   "Selected bandwidth equals the smallest or largest grid value.")

})

test_that("CV bandwidth selection is correctly implemented", {

    # To search for the minimum, the package locpol uses a combination of golden
    # section search and successive parabolic interpolation.
    # In contrast, we perform explicit grid search.

    # Nadaraya-Watson (degree = 0)
    bw_CV_NW_own <- suppressWarnings(bw_CV_fun(X = X, Y = Y, degree = 0L, plot = FALSE,
                                               bw_grid = bw_grid))
    bw_CV_NW_locpol <- regCVBwSelC(x = X, y = Y, deg = 0, kernel = EpaK,
                                   interval = interval)

    expect_equal(round(bw_CV_NW_own, 3), round(bw_CV_NW_locpol, 3))

    # Local linear (degree = 1)
    bw_CV_LL_own <- bw_CV_fun(X = X, Y = Y, degree = 1L, plot = FALSE,
                              bw_grid = bw_grid)
    bw_CV_LL_locpol <- regCVBwSelC(x = X, y = Y, deg = 1, kernel = EpaK,
                                   interval = interval)

    expect_equal(round(bw_CV_LL_own, 3), round(bw_CV_LL_locpol, 3))

    # Higher degree (degree = 3)
    bw_CV_cubic_own <- suppressMessages(bw_CV_fun(X = X, Y = Y, degree = 3L, plot = FALSE,
                                                  bw_grid = bw_grid))
    bw_CV_cubic_locpol <- regCVBwSelC(x = X, y = Y, deg = 3, kernel = EpaK,
                                      interval = interval)

    expect_equal(round(bw_CV_cubic_own, 3), round(bw_CV_cubic_locpol, 3))

})