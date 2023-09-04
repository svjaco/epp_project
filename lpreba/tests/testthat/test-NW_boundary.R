################################################################################
#        Tests for NW_boundary function (Boundary-adjusted NW estimator)       #
################################################################################

test_that("error messages appear", {

    expect_error(NW_boundary(x = X, X = X[1:(n-1)], Y = Y, bw = bw,
                             boundary_left = 0, boundary_right = 1),
                 "Length of X and Y has to be equal.")

    expect_error(NW_boundary(x = X, X = X, Y = Y, bw = -bw,
                             boundary_left = 0, boundary_right = 1),
                 "Bandwidth has to be strictly positive.")

    expect_error(NW_boundary(x = X, X = X, Y = Y, bw = bw,
                             kernel_interior = dnorm, # dnorm is the Gaussian kernel
                             boundary_left = 0, boundary_right = 1),
                 "Only supported kernels can be chosen.")

    expect_error(NW_boundary(x = X, X = X, Y = Y, bw = bw,
                             kernel_left = dnorm,
                             boundary_left = 0, boundary_right = 1),
                 "Only supported boundary kernels can be chosen.")

    expect_error(NW_boundary(x = X, X = X, Y = Y, bw = bw,
                             boundary_left = 0.5, boundary_right = 1),
                 "The lower boundary cannot be larger than the smallest observed value of the regressor X.")

    expect_error(NW_boundary(x = X, X = X, Y = Y, bw = bw,
                             boundary_left = 0, boundary_right = 0.5),
                 "The upper boundary cannot be smaller than the largest observed value of the regressor X.")

})