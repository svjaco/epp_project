################################################################################
#                 Tests for CV_error_fun function (LOOCV error)                #
################################################################################

test_that("error messages appear", {

    expect_error(CV_error_fun(X = X, Y = Y, bw = bw, degree = 1L,
                              boundary_adjustment = TRUE),
                 "Explicit boundary adjustment can only be used for the Nadaraya-Watson estimator (degree = 0).",
                 fixed = TRUE)

})