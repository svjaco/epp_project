################################################################################
#              Tests for LP function (Local polynomial estimator)              #
################################################################################

# We use the package "locpol" (published on CRAN) to validate our
# implementation of the local polynomial estimator.
# Due to one issue with this package (see below),
# also the package "np" is used.

test_that("error messages appear", {

    expect_error(LP(x = X, X = X[1:(n-1)], Y = Y, bw = bw),
                 "Length of X and Y has to be equal.")

    expect_error(LP(x = X, X = X, Y = Y, bw = -bw),
                 "Bandwidth has to be strictly positive.")

    expect_error(LP(x = X, X = X, Y = Y, kernel = stats::dnorm, bw = bw), # dnorm is the Gaussian kernel
                 "Only supported kernels can be chosen.")

    expect_error(LP(x = X, X = X, Y = Y, bw = bw, degree = -1L),
                 "degree has to be of type integer and nonnegative.")

    expect_error(LP(x = X, X = X, Y = Y, bw = bw, degree = 0.5),
                 "degree has to be of type integer and nonnegative.")

})

test_that("warning messages appear", {

    expect_warning(LP(x = X, X = X, Y = Y, bw = bw, degree = 0L),
                   message(strwrap("It is strongly recommended to use an
                                   odd order as odd order polynomial fits
                                   outperform even order ones.
                                   In particular, odd order fits are
                                   boundary adaptive.",
                                   prefix = " ", initial = "")))

})

test_that("general messages appear", {

    expect_message(LP(x = X, X = X, Y = Y, bw = bw, degree = 3L),
                   message(strwrap("It is recommended to use the lowest (odd) order,
                                   since the bandwidth is used to control the
                                   model complexity.",
                                   prefix = " ", initial = "")))

})

test_that("estimation of the regression function is correctly implemented", {

    # Nadaraya-Watson (degree = 0)

    #output_NW_LP <- LP(x = X, X = X, Y = Y, bw = bw, degree = 0L)
    #estimates_NW_LP <- output_NW_LP$estimates
    #output_NW_locpol <- locpol(Y ~ X, df, bw = bw, deg = 0) # locpol uses the Epanechnikov kernel by default
    #estimates_NW_locpol <- output_NW_locpol$lpFit$Y

    output_NW_LP <- suppressWarnings(LP(x = X, X = X, Y = Y, kernel = uniform, bw = bw, degree = 0L))
    estimates_NW_LP <- output_NW_LP$estimates
    output_NW_npreg <- npreg(Y ~ X, ckertype = "uniform", bws = bw, regtype = "lc")
    estimates_NW_npreg <- output_NW_npreg$mean

    expect_equal(estimates_NW_LP, estimates_NW_npreg)

    # The direct computation of the Nadaraya-Watson estimates via the locpol
    # package results in an error.
    # We informed the author, but he did not reply.
    # As an alternative, the np package with the npreg function is applied.

    # Local linear (degree = 1)
    output_LL_LP <- LP(x = X, X = X, Y = Y, bw = bw, degree = 1L)
    estimates_LL_LP <- output_LL_LP$estimates
    output_LL_locpol <- locpol(Y ~ X, df, bw = bw, deg = 1)
    estimates_LL_locpol <- output_LL_locpol$lpFit$Y

    expect_equal(estimates_LL_LP, estimates_LL_locpol)

    # Higher degree (degree = 3)
    output_cubic_LP <- LP(x = X, X = X, Y = Y, bw = bw, degree = 3L)
    estimates_cubic_LP <- output_cubic_LP$estimates
    output_cubic_locpol <- locpol(Y ~ X, df, bw = bw, deg = 3)
    estimates_cubic_locpol <- output_cubic_locpol$lpFit$Y

    expect_equal(estimates_cubic_LP, estimates_cubic_locpol)

})

test_that("estimation for the first derivative of the regression function is correctly implemented", {

    # Local linear (degree = 1)
    output_LL_LP <- LP(x = X, X = X, Y = Y, bw = bw, degree = 1L)
    slopes_LL_LP <- output_LL_LP$slopes
    output_LL_locpol <- locpol(Y ~ X, df, bw = bw, deg = 1)
    slopes_LL_locpol <- output_LL_locpol$lpFit$Y1

    expect_equal(slopes_LL_LP, slopes_LL_locpol)

    # Higher degree (degree = 3)
    output_cubic_LP <- LP(x = X, X = X, Y = Y, bw = bw, degree = 3L)
    slopes_cubic_LP <- output_cubic_LP$slopes
    output_cubic_locpol <- locpol(Y ~ X, df, bw = bw, deg = 3)
    slopes_cubic_locpol <- output_cubic_locpol$lpFit$Y1

    expect_equal(slopes_cubic_LP, slopes_cubic_locpol)

})

test_that("estimation for the second derivative of the regression function is correctly implemented", {

    # Higher degree (degree = 3)
    output_cubic_LP <- LP(x = X, X = X, Y = Y, bw = bw, degree = 3L)
    curvatures_cubic_LP <- output_cubic_LP$curvatures
    output_cubic_locpol <- locpol(Y ~ X, df, bw = bw, deg = 3)
    curvatures_cubic_locpol <- output_cubic_locpol$lpFit$Y2

    expect_equal(curvatures_cubic_LP, curvatures_cubic_locpol)

})