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
# implementation of the local polynomial estimator.

################################################################################
#              Tests for LP function (Local polynomial estimator)              #
################################################################################

#################################
#    Generating example data    #
#################################

m_fun <- function(x) {sin(2*pi*x)} # True regression function
n <- 100 # Sample size
X <- seq(0, 1, length.out = n) # Observed values for regressor
m_X <- m_fun(X) # True values of regression function
epsilon <- rnorm(n, sd = 0.25) # Error term
Y <- m_X + epsilon # Observed values for regressand
h <- 0.2 # Bandwidth
df <- data.frame(Y, X) # data frame, required for locpol

#########################################################
#    Checks for messages (errors, warnings, general)    #
#########################################################

# Error messages
LP_unequal_length <- LP(x = X, X = X[1:(n-1)], Y = Y, bw = h)
LP_negative_bw <- LP(x = X, X = X, Y = Y, bw = -h)
LP_non-supported_kernel <- LP(x = X, X = X, Y = Y, kernel = dnorm, bw = h) # dnorm is the Gaussian kernel
LP_negative_degree <- LP(x = X, X = X, Y = Y, bw = h, degree = -1L)
LP_non-integer_degree <- LP(x = X, X = X, Y = Y, bw = h, degree = 0.5)

# Warning messages
LP_even_degree <- LP(x = X, X = X, Y = Y, bw = h, degree = 0L)

# General messages
LP_higher_degree <- LP(x = X, X = X, Y = Y, bw = h, degree = 3L)

##########################################################
#    Checks for estimates for the regression function    #
##########################################################

# Nadaraya-Watson (degree = 0)
output_NW_LP <- LP(x = X, X = X, Y = Y, bw = h, degree = 0L)
estimates_NW_LP <- output_NW_LP$estimates
output_NW_locpol <- locpol(Y ~ X, df, bw = h, deg = 0) # locpol uses the Epanechnikov kernel by default
estimates_NW_locpol <- output_NW_locpol$lpFit$Y
all.equal(estimates_NW_LP, estimates_NW_locpol) # Should yield TRUE

# Local linear (degree = 1)
output_LL_LP <- LP(x = X, X = X, Y = Y, bw = h, degree = 1L)
estimates_LL_LP <- output_LL_LP$estimates
output_LL_locpol <- locpol(Y ~ X, df, bw = h, deg = 1)
estimates_LL_locpol <- output_LL_locpol$lpFit$Y
all.equal(estimates_LL_LP, estimates_LL_locpol) # Should yield TRUE

# Higher degree (degree = 3)
output_cubic_LP <- LP(x = X, X = X, Y = Y, bw = h, degree = 3L)
estimates_cubic_LP <- output_cubic_LP$estimates
output_cubic_locpol <- locpol(Y ~ X, df, bw = h, deg = 3)
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