# Local polynomial regression with the option for explicit boundary adjustment

# Effective Programming Practices for Economists
# Prof. von Gaudecker
# Winter 2021/22, M.Sc. Economics, Bonn University
# Sven Jacobs

################################################################################

#############################
#    Auxiliary functions    #
#############################

# Creating standard basis vectors
make_basis <- function(length, position) {

    replace(numeric(length), position, 1)

}

#################
#    Kernels    #
#################

# Uniform
uniform <- function(u) {ifelse(abs(u) <= 1, 0.5, 0)}

# Triangular
triangular <- function(u) {ifelse(abs(u) <= 1, 1 - abs(u), 0)}

# Epanechnikov
epanechnikov <- function(u) {ifelse(abs(u) <= 1, 3/4 * (1 - u^2), 0)}

# Biweight
biweight <- function(u) {ifelse(abs(u) <= 1, 15/16 * (1 - u^2)^2, 0)}

# Triweight
triweight <- function(u) {ifelse(abs(u) <= 1, 35/32 * (1 - u^2)^3, 0)}

# Tricube
tricube <- function(u) {ifelse(abs(u) <= 1, 70/81 * (1 - abs(u)^3)^3, 0)}

# Cosine
cosine <- function(u) {ifelse(abs(u) <= 1, pi/4 * cos(pi/2 * u), 0)}

# List with all kernels
kernel_list <- list(uniform, triangular, epanechnikov, biweight,
                    triweight, tricube, cosine)

####################################
#    Local polynomial estimator    #
####################################

LP <- function(x, X, Y, kernel = epanechnikov, bw, degree = 1L) {

    # input: - x: evaluation points (vector)
    #        - X: data for the regressor (vector)
    #        - Y: data for the regressand (vector)
    #        - kernel: kernel (function), with default
    #        - bw: bandwidth (scalar)
    #        - degree: degree of the locally fitted polynomial (integer),
    #                  with default

    # output: - list containing
    #               - estimates: estimates for the regression function
    #                            at the evaluation points (vector)
    #               - effective_kernels: effective kernels
    #                                    at the evaluation points (matrix)
    #               - slopes: estimates for the first derivative (slope) of the
    #                         regression function at the evaluation points (vector)
    #               - curvatures: estimates for the second derivative (curvature) of the
    #                             regression function at the evaluation points (vector)

    # Error messages
    if (length(X) != length(Y)) {stop("Length of X and Y has to be equal.")}
    if (bw <= 0) {stop("Bandwidth has to be strictly positive.")}
    if (TRUE %in% lapply(kernel_list, all.equal, kernel) == FALSE) {
        stop("Only supported kernels can be chosen.")
    }
    if (degree < 0 | typeof(degree) != "integer") {
        stop("degree has to be of type integer and nonnegative.")
    }

    # Warning messages
    if ((degree %% 2) == 0) {warning(strwrap("It is strongly recommended to use an
                                             odd order as odd order polynomial fits
                                             outperform even order ones.
                                             In particular, odd order fits are
                                             boundary adaptive.",
                                             prefix = " ", initial = ""))}

    # General messages
    if (degree > 1) {message(strwrap("It is recommended to use the lowest (odd) order,
                                     since the bandwidth is used to control the
                                     model complexity.",
                                     prefix = " ", initial = ""))}

    X_matrix <- matrix(NA, length(X), degree + 1)

    estimates <- rep(NA, length(x))
    effective_kernels <- matrix(NA, length(x), length(Y))
    slopes <- rep(NA, length(x))
    curvatures <- rep(NA, length(x))

    for (i in 1:length(x)) {

        for (j in 1:(degree + 1)) {

            X_matrix[, j] <- (X - x[i])^(j - 1)

        }

        W_matrix <- diag(kernel((X - x[i])/bw))

        auxiliary_matrix <- solve(t(X_matrix)%*%W_matrix%*%X_matrix) %*%
                            t(X_matrix) %*% W_matrix

        effective_kernels[i, ] <- t(make_basis(degree + 1, 1)) %*% auxiliary_matrix
        estimates[i] <- effective_kernels[i, ] %*% Y

        if (degree >= 1) {
            slopes[i] <- t(make_basis(degree + 1, 2)) %*% auxiliary_matrix %*% Y
        }
        if (degree >= 2) {
            curvatures[i] <- 2 * t(make_basis(degree + 1, 3)) %*% auxiliary_matrix %*% Y
        }

    }

    list(estimates = estimates, effective_kernels = effective_kernels,
         slopes = slopes, curvatures = curvatures)

}

##########################
#    Boundary kernels    #
##########################

# Note:
#
# The following boundary kernels rely on Müller (1991).
# However, the left boundary kernels of Müller are right boundary kernels for us, and vice versa.
# This is because of the different definition of the kernel estimator concerning the argument of the kernel.
# Müller uses K(x - X), whereas we define K(X - x).

# Left Uniform kernels
uniform_left <- function(u, rho) {

    ifelse(u >= -rho & u <= 1, 1/(1+rho)*(1 + 3*((1-rho)/(1+rho))^2 - 6*((1-rho)/(1+rho)^2)*u), 0)

}

# Left Epanechnikov kernels
epanechnikov_left <- function(u, rho) {

    ifelse(u >= -rho & u <= 1, 6*(1-u)*(rho+u)*(1/(1+rho)^3)*(1 + 5*((1-rho)/(1+rho))^2 - 10*((1-rho)/(1+rho)^2)*u), 0)

}

# Left Biweight kernels
biweight_left <- function(u, rho) {

    ifelse(u >= -rho & u <= 1, 30*(1-u)^2*(rho+u)^2*(1/(1+rho)^5)*(1 + 7*((1-rho)/(1+rho))^2 - 14*((1-rho)/(1+rho)^2)*u), 0)

}

# Left Triweight kernels
triweight_left <- function(u, rho) {

    ifelse(u >= -rho & u <= 1, 140*(1-u)^3*(rho+u)^3*(1/(1+rho)^7)*(1 + 9*((1-rho)/(1+rho))^2 - 18*((1-rho)/(1+rho)^2)*u), 0)

}

# Note:
#
# The following boundary kernels rely on Müller and Wang (1994).
# These alternative kernels have greater asymptotic efficiency
# at the cost of discontinuity at their endpoints.

# Left alternative Uniform kernels
uniform_alt_left <- function(u, rho) {

    ifelse(u >= -rho & u <= 1, 2/(1+rho)^3*(-3*(1-rho)*u + 2*(1-rho+rho^2)), 0)

}

# Left alternative Epanechnikov kernels
epanechnikov_alt_left <- function(u, rho) {

    ifelse(u >= -rho & u <= 1, 12/(1+rho)^4*(1-u)*(-(1-2*rho)*u + 0.5*(3*rho^2-2*rho+1)), 0)

}

# Left alternative Biweight kernels
biweight_alt_left <- function(u, rho) {

    ifelse(u >= -rho & u <= 1, 15/(1+rho)^5*(1-u)^2*(rho+u)*(-2*u*(5*(1-rho)/(1+rho)-1) + 3*rho-1 + 5*(1-rho)^2/(1+rho)), 0)

}

# Left alternative Triweight kernels
triweight_alt_left <- function(u, rho) {

    ifelse(u >= -rho & u <= 1, 70/(1+rho)^7*(1-u)^3*(rho+u)^2*(-2*u*(7*(1-rho)/(1+rho)-1) + 3*rho-1 + 7*(1-rho)^2/(1+rho)), 0)

}

# List with all boundary kernels
boundary_kernel_list <- list(uniform_left, epanechnikov_left, biweight_left, triweight_left,
                             uniform_alt_left, epanechnikov_alt_left, biweight_alt_left, triweight_alt_left)

#####################################################
#    Boundary-adjusted Nadaraya-Watson estimator    #
#####################################################

NW_boundary <- function(x, X, Y, bw,
                        kernel_interior = epanechnikov, kernel_left = epanechnikov_left,
                        boundary_left, boundary_right) {

    # input: - x: evaluation points (vector)
    #        - X: data for the regressor (vector)
    #        - Y: data for the regressand (vector)
    #        - bw: bandwidth (scalar)
    #        - kernel_interior: kernel used in the interior (function), with default
    #        - kernel_left: left boundary kernels (function), with default
    #        - boundary_left: lower boundary of the support of X (scalar)
    #        - boundary_right: upper boundary of the support of X (scalar)

    # output: - list containing
    #               - estimates: estimates for the regression function
    #                            at the evaluation points (vector)
    #               - effective_kernels: effective kernels
    #                                    at the evaluation points (matrix)

    # Error messages
    if (length(X) != length(Y)) {stop("Length of X and Y has to be equal.")}
    if (bw <= 0) {stop("Bandwidth has to be strictly positive.")}
    if (TRUE %in% lapply(kernel_list, all.equal, kernel_interior) == FALSE) {
        stop("Only supported kernels can be chosen.")
    }
    if (TRUE %in% lapply(boundary_kernel_list, all.equal, kernel_left) == FALSE) {
        stop("Only supported boundary kernels can be chosen.")
    }
    if (boundary_left > min(X)) {
        stop("The lower boundary cannot be larger than the smallest observed value of the regressor X.")
    }
    if (boundary_right < max(X)) {
        stop("The upper boundary cannot be smaller than the largest observed value of the regressor X.")
    }

    absolute_weights <- matrix(NA, length(x), length(Y))

    for (i in 1:length(x)) {

        if (x[i] >= boundary_left & x[i] < boundary_left + bw) {

            absolute_weights[i, ] <- kernel_left(u = (X - x[i])/bw, rho = (x[i] - boundary_left)/bw)

        } else if (x[i] >= boundary_left + bw & x[i] <= boundary_right - bw) {

            absolute_weights[i, ] <- kernel_interior((X - x[i])/bw)

        } else if (x[i] > boundary_right - bw & x[i] <= boundary_right) {

            absolute_weights[i, ] <- kernel_left(u = -(X - x[i])/bw, rho = (boundary_right - x[i])/bw)

        }

    }

    effective_kernels <- absolute_weights/rowSums(absolute_weights)

    estimates <- as.vector(effective_kernels %*% Y)

    list(estimates = estimates, effective_kernels = effective_kernels)

}