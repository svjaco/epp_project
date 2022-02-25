# Tests to check the correct functioning of the functions

# Effective Programming Practices for Economists
# Prof. von Gaudecker
# Winter 2021/22, M.Sc. Economics, Bonn University
# Sven Jacobs

################################################################################

################################################################################
#                                 Preparation                                  #
################################################################################

set.seed(123) # Seed for reproducibility

# Packages
library(locpol) # Kernel local polynomial regression with weights

# We use the package "locpol" (published on CRAN) to validate our
# implementation of the local polynomial estimator and
# the cross-validated bandwidth selector.

################################################################################
#                           Generating example data                            #
################################################################################

m_fun <- function(x) {sin(2*pi*x)} # True regression function
n <- 100 # Sample size
X <- seq(0, 1, length.out = n) # Observed values for regressor
m_X <- m_fun(X) # True values of regression function
epsilon <- rnorm(n, sd = 0.25) # Error term
Y <- m_X + epsilon # Observed values for regressand
bw <- 0.2 # Bandwidth
df <- data.frame(Y, X) # data frame, required for locpol

################################################################################
#              Tests for LP function (Local polynomial estimator)              #
################################################################################

#########################################################
#    Checks for messages (errors, warnings, general)    #
#########################################################

# Error messages
LP_unequal_length <- LP(x = X, X = X[1:(n-1)], Y = Y, bw = bw)
LP_negative_bw <- LP(x = X, X = X, Y = Y, bw = -bw)
LP_non-supported_kernel <- LP(x = X, X = X, Y = Y, kernel = dnorm, bw = bw) # dnorm is the Gaussian kernel
LP_negative_degree <- LP(x = X, X = X, Y = Y, bw = bw, degree = -1L)
LP_non-integer_degree <- LP(x = X, X = X, Y = Y, bw = bw, degree = 0.5)

# Warning messages
LP_even_degree <- LP(x = X, X = X, Y = Y, bw = bw, degree = 0L)

# General messages
LP_higher_degree <- LP(x = X, X = X, Y = Y, bw = bw, degree = 3L)

##########################################################
#    Checks for estimates for the regression function    #
##########################################################

# Nadaraya-Watson (degree = 0)
output_NW_LP <- LP(x = X, X = X, Y = Y, bw = bw, degree = 0L)
estimates_NW_LP <- output_NW_LP$estimates
output_NW_locpol <- locpol(Y ~ X, df, bw = bw, deg = 0) # locpol uses the Epanechnikov kernel by default
estimates_NW_locpol <- output_NW_locpol$lpFit$Y
all.equal(estimates_NW_LP, estimates_NW_locpol) # Should yield TRUE

# Local linear (degree = 1)
output_LL_LP <- LP(x = X, X = X, Y = Y, bw = bw, degree = 1L)
estimates_LL_LP <- output_LL_LP$estimates
output_LL_locpol <- locpol(Y ~ X, df, bw = bw, deg = 1)
estimates_LL_locpol <- output_LL_locpol$lpFit$Y
all.equal(estimates_LL_LP, estimates_LL_locpol) # Should yield TRUE

# Higher degree (degree = 3)
output_cubic_LP <- LP(x = X, X = X, Y = Y, bw = bw, degree = 3L)
estimates_cubic_LP <- output_cubic_LP$estimates
output_cubic_locpol <- locpol(Y ~ X, df, bw = bw, deg = 3)
estimates_cubic_locpol <- output_cubic_locpol$lpFit$Y
all.equal(estimates_cubic_LP, estimates_cubic_locpol) # Should yield TRUE

##################################################################################
#    Checks for estimates for the first derivative of the regression function    #
##################################################################################

# Local linear (degree = 1)
slopes_LL_LP <- output_LL_LP$slopes
slopes_LL_locpol <- output_LL_locpol$lpFit$Y1
all.equal(slopes_LL_LP, slopes_LL_locpol) # Should yield TRUE

# Higher degree (degree = 3)
slopes_cubic_LP <- output_cubic_LP$slopes
slopes_cubic_locpol <- output_cubic_locpol$lpFit$Y1
all.equal(slopes_cubic_LP, slopes_cubic_locpol) # Should yield TRUE

###################################################################################
#    Checks for estimates for the second derivative of the regression function    #
###################################################################################

# Higher degree (degree = 3)
curvatures_cubic_LP <- output_cubic_LP$curvatures
curvatures_cubic_locpol <- output_cubic_locpol$lpFit$Y2
all.equal(curvatures_cubic_LP, curvatures_cubic_locpol) # Should yield TRUE

################################################################################
#        Tests for NW_boundary function (Boundary-adjusted NW estimator)       #
################################################################################

######################################
#    Checks for messages (errors)    #
######################################

NW_boundary_unequal_length <- NW_boundary(x = X, X = X[1:(n-1)], Y = Y, bw = bw,
                                          boundary_left = 0, boundary_right = 1)
NW_boundary_negative_bw <- NW_boundary(x = X, X = X, Y = Y, bw = -bw,
                                       boundary_left = 0, boundary_right = 1)
NW_boundary_non-supported_kernel <- NW_boundary(x = X, X = X, Y = Y, bw = bw,
                                                kernel_interior = dnorm,
                                                boundary_left = 0, boundary_right = 1)
NW_boundary_non-supported_boundary_kernel <- NW_boundary(x = X, X = X, Y = Y, bw = bw,
                                                         kernel_left = dnorm,
                                                         boundary_left = 0, boundary_right = 1)
NW_boundary_left_boundary_too_large <- NW_boundary(x = X, X = X, Y = Y, bw = bw,
                                                   boundary_left = 0.5, boundary_right = 1)
NW_boundary_right_boundary_too_small <- NW_boundary(x = X, X = X, Y = Y, bw = bw,
                                                    boundary_left = 0, boundary_right = 0.5)

################################################################################
#                 Tests for CV_error_fun function (LOOCV error)                #
################################################################################

######################################
#    Checks for messages (errors)    #
######################################

CV_error_disallowed_boundary_adjustment <- CV_error_fun(X = X, Y = Y, bw = bw, degree = 1L,
                                                        boundary_adjustment = TRUE)

################################################################################
#              Tests for bw_CV_fun function (CV optimal bandwidth)             #
################################################################################

########################################
#    Checks for messages (warnings)    #
########################################

# Selected bandwidth coincides with the smallest grid value (see also plot)
bw_CV_fun(X = X, Y = Y, bw_grid = seq(0.075, 0.2, length.out = 100))
# Expanding the grid reveals that global minimum was not reached (see also plot)
bw_CV_fun(X = X, Y = Y, bw_grid = seq(0.05, 0.2, length.out = 100))

#########################################
#    Checks for CV optimal bandwidth    #
#########################################

# Note:
#
# To search for the minimum, the package locpol uses a combination of golden
# section search and successive parabolic interpolation.
# In contrast, we perform explicit grid search.

bw_grid <- seq(0.05, 0.2, length.out = 1000)
interval <- c(min(bw_grid), max(bw_grid))

# Nadaraya-Watson (degree = 0)
bw_CV_NW_own <- suppressWarnings(bw_CV_fun(X = X, Y = Y, degree = 0L, plot = TRUE,
                                 bw_grid = bw_grid))
bw_CV_NW_locpol <- regCVBwSelC(x = X, y = Y, deg = 0, kernel = EpaK,
                               interval = interval)
all.equal(round(bw_CV_NW_own, 3), round(bw_CV_NW_locpol, 3)) # Should yield TRUE

# Local linear (degree = 1)
bw_CV_LL_own <- bw_CV_fun(X = X, Y = Y, degree = 1L, plot = TRUE,
                          bw_grid = bw_grid)
bw_CV_LL_locpol <- regCVBwSelC(x = X, y = Y, deg = 1, kernel = EpaK,
                               interval = interval)
all.equal(round(bw_CV_LL_own, 3), round(bw_CV_LL_locpol, 3)) # Should yield TRUE

# Higher degree (degree = 3)
bw_CV_cubic_own <- suppressMessages(bw_CV_fun(X = X, Y = Y, degree = 3L, plot = TRUE,
                                    bw_grid = bw_grid))
bw_CV_cubic_locpol <- regCVBwSelC(x = X, y = Y, deg = 3, kernel = EpaK,
                                  interval = interval)
all.equal(round(bw_CV_cubic_own, 3), round(bw_CV_cubic_locpol, 3)) # Should yield TRUE