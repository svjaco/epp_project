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
    if(length(X) != length(Y)) {stop("Length of X and Y has to be equal.")}
    if(bw <= 0) {stop("Bandwidth has to be strictly positive.")}
    if(TRUE %in% lapply(kernel_list, all.equal, kernel) == FALSE) {
        stop("Only supported kernels can be chosen.")
    }
    if(degree < 0 | typeof(degree) != "integer") {
        stop("degree has to be of type integer and nonnegative.")
    }

    # Warning messages
    if((degree %% 2) == 0) {warning(strwrap("It is strongly recommended to use an
                                            odd order as odd order polynomial fits
                                            outperform even order ones.
                                            In particular, odd order fits are
                                            boundary adaptive.",
                                            prefix = " ", initial = ""))}

    # General messages
    if(degree > 1) {message(strwrap("It is recommended to use the lowest (odd) order,
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

        if(degree >= 1) {
            slopes[i] <- t(make_basis(degree + 1, 2)) %*% auxiliary_matrix %*% Y
        }
        if(degree >= 2) {
            curvatures[i] <- 2 * t(make_basis(degree + 1, 3)) %*% auxiliary_matrix %*% Y
        }

    }

    list(estimates = estimates, effective_kernels = effective_kernels,
         slopes = slopes, curvatures = curvatures)

}