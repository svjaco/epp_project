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
# The following boundary kernels rely on M?ller (1991).
# However, the left boundary kernels of M?ller are right boundary kernels for us, and vice versa.
# This is because of the different definition of the kernel estimator concerning the argument of the kernel.
# M?ller uses K(x - X), whereas we define K(X - x).

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
# The following boundary kernels rely on M?ller and Wang (1994).
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

######################################################
#    Leave-one-out cross-validation (LOOCV) error    #
######################################################

CV_error_fun <- function(X, Y, kernel = epanechnikov, bw, degree = 1L,
                         kernel_left = epanechnikov_left,
                         boundary_left = NA, boundary_right = NA,
                         boundary_adjustment = FALSE) {

    # input: - X: data for the regressor (vector)
    #        - Y: data for the regressand (vector)
    #        - kernel: kernel (function), with default
    #        - bw: bandwidth (scalar)
    #        - degree: degree of the locally fitted polynomial (integer), with default
    #        - kernel_left: left boundary kernels (function), with default
    #        - boundary_left: lower boundary of the support of X (scalar), with default
    #        - boundary_right: upper boundary of the support of X (scalar), with default
    #        - boundary_adjustment: explicit boundary adjustment (boolean), with default

    # output: - CV_error: LOOCV error for bandwidth bw (scalar)

    # Error messages
    if (degree != 0L & boundary_adjustment == TRUE) {
        stop("Explicit boundary adjustment can only be used for the Nadaraya-Watson estimator (degree = 0).")
    }

    if (boundary_adjustment == FALSE) {

        LP_output <- LP(x = X, X = X, Y = Y, kernel = kernel, bw = bw, degree = degree)
        estimates <- LP_output$estimates
        effective_kernels <- LP_output$effective_kernels

    } else if (boundary_adjustment == TRUE) {

        NW_boundary_output <- NW_boundary(x = X, X = X, Y = Y, bw = bw,
                                          kernel_interior = kernel, kernel_left = kernel_left,
                                          boundary_left = boundary_left, boundary_right = boundary_right)
        estimates <- NW_boundary_output$estimates
        effective_kernels <- NW_boundary_output$effective_kernels

    }

    CV_error <- 1/length(Y) * sum(((Y - estimates)/(1 - diag(effective_kernels)))^2)
    return(CV_error)

}

##############################
#    CV optimal bandwidth    #
##############################

bw_CV_fun <- function(X, Y, kernel = epanechnikov, bw_grid, degree = 1L,
                      kernel_left = epanechnikov_left,
                      boundary_left = NA, boundary_right = NA,
                      boundary_adjustment = FALSE,
                      plot = TRUE) {

    # input: - X: data for the regressor (vector)
    #        - Y: data for the regressand (vector)
    #        - kernel: kernel (function), with default
    #        - bw_grid: bandwidth grid to find global minimum of CV error (vector)
    #        - degree: degree of the locally fitted polynomial (integer), with default
    #        - kernel_left: left boundary kernels (function), with default
    #        - boundary_left: lower boundary of the support of X (scalar), with default
    #        - boundary_right: upper boundary of the support of X (scalar), with default
    #        - boundary_adjustment: explicit boundary adjustment (boolean), with default
    #        - plot: plot CV error over bandwidth grid (boolean), with default

    # output: - bw_CV: CV optimal bandwidth (scalar)
    #         - plot of CV errors if plot == TRUE (plot)

    CV_errors <- sapply(bw_grid, function(bw) {

        CV_error_fun(X = X, Y = Y, kernel = kernel, bw, degree = degree,
                     kernel_left = kernel_left,
                     boundary_left = boundary_left, boundary_right = boundary_right,
                     boundary_adjustment = boundary_adjustment)

    })

    bw_CV <- bw_grid[which.min(CV_errors)]

    if (bw_CV == bw_grid[1] | bw_CV == bw_grid[length(bw_grid)]) {
        warning("Selected bandwidth equals the smallest or largest grid value.
                You may wish to expand the grid after consulting the plot of the CV errors.")
    }

    if (plot == TRUE) {

        plot(bw_grid, CV_errors,
             xlab = "bw", ylab = "CV error", cex.lab = 1.25, cex.axis = 1.25, cex = 0.75)
        rug(bw_grid, ticksize = 0.015)
        abline(v = bw_CV, col = "red")
        legend("top", legend = bquote(paste("bw"["CV"]*" = ", .(round(bw_CV, 3)))), bty = "n", cex = 1.25)

    }

    return(bw_CV)

}

############################################################################
#    Asymptotic confidence intervals for the local polynomial estimator    #
############################################################################

confidence_intervals_LP <- function(x, X, Y, kernel = epanechnikov, bw, degree = 1L, alpha = 0.05) {

    # input: - x: evaluation points (vector)
    #        - X: data for the regressor (vector)
    #        - Y: data for the regressand (vector)
    #        - kernel: kernel (function), with default
    #        - bw: bandwidth (scalar)
    #        - degree: degree of the locally fitted polynomial (integer), with default
    #        - alpha: significance level (scalar), with default

    # output: - list containing
    #               - confidence_intervals_lower:
    #                 lower points for the (1 - alpha) asymptotic confidence intervals (vector)
    #               - confidence_intervals_upper:
    #                 upper points for the (1 - alpha) asymptotic confidence intervals (vector)

    # Error messages
    if (alpha <= 0 | alpha >= 1) {
        stop("The significance level has to be between zero and one.")
    }

    LP_output <- LP(x = x, X = X, Y = Y, kernel = kernel, bw = bw, degree = degree)
    estimates <- LP_output$estimates
    effective_kernels <- LP_output$effective_kernels

    squared_prediction_errors <- ((Y - estimates)/(1 - diag(effective_kernels)))^2
    error_matrix <- diag(squared_prediction_errors)

    X_matrix <- matrix(NA, length(X), degree + 1)
    estimates_variance <- rep(NA, length(x))

    for (i in 1:length(x)) {

        for (j in 1:(degree + 1)) {

            X_matrix[, j] <- (X - x[i])^(j - 1)

        }

        W_matrix <- diag(kernel((X - x[i])/bw))

        auxiliary_matrix <- solve(t(X_matrix)%*%W_matrix%*%X_matrix)
        estimate_covariance_matrix <- auxiliary_matrix %*%
                                      t(X_matrix)%*%W_matrix%*%error_matrix%*%W_matrix%*%X_matrix %*%
                                      auxiliary_matrix
        estimates_variance[i] <- estimate_covariance_matrix[1, 1]

    }

    confidence_intervals_LP_lower <- estimates - qnorm(p = 1 - alpha/2) * sqrt(estimates_variance)
    confidence_intervals_LP_upper <- estimates + qnorm(p = 1 - alpha/2) * sqrt(estimates_variance)

    list(confidence_intervals_lower = confidence_intervals_LP_lower,
         confidence_intervals_upper = confidence_intervals_LP_upper)

}
