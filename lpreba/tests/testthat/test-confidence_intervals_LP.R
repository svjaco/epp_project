################################################################################
# Tests for confidence_intervals_LP function (Asymptotic confidence intervals) #
################################################################################

test_that("error messages appear", {

    expect_error(confidence_intervals_LP(x = X, X = X, Y = Y, bw = bw, alpha = 0),
                 "The significance level has to be between zero and one.")

})

# We compute the standard errors for the construction of asymptotic confidence
# intervals for the local polynomial estimator according to Hansen (2022, Section 19.16).
# Currently, there is no R package available with such an implementation.
#
# Hansen, B. (2022). Econometrics. New Jersey: Princeton University Press.